# ----------------------------------------------------------------------------
# Copyright (c) 2016-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import json
import os
import pkg_resources
import shutil
from urllib.parse import quote
import functools
import matplotlib.pyplot as plt

import scipy
import numpy as np
import pandas as pd
import qiime2
from statsmodels.sandbox.stats.multicomp import multipletests
import q2templates
import biom
import itertools
from qiime2 import Artifact

from kneed import KneeLocator
#from qiime2.plugins.diversity.methods import alpha
#from . import METRICS
from q2_types.tree import NewickFormat


#helper functions
def subsample_feature_table(feature_table, depths):
    """
    Perform rarefaction on the feature table for multiple rarefaction depths.
    
    Parameters:
    feature_table (pd.DataFrame): The feature table with samples as rows and species as columns.
                                  Values represent the abundance of each species in each sample.
    depths (list): A list of rarefaction depths to subsample to.
    
    Returns:
    pd.DataFrame: A DataFrame where each entry is the rarefied abundance data for each depth.
                  The resulting DataFrame will have the same number of samples and species.
    """
    rarefied_data = []

    for depth in depths:
        
        for sample in feature_table.index:
            sample_counts = feature_table.loc[sample].values  # Get species counts for this sample
            
            if sample_counts.sum() >= depth:
                rarefied_sample = rarefy(sample_counts, depth)
                rarefied_data.append((depth, sample, rarefied_sample))
            else:
                # If a sample has fewer reads than the rarefaction depth, set it to NaN 
                rarefied_data.append((depth, sample, np.nan * np.ones(feature_table.shape[1])))

        
    index = pd.MultiIndex.from_tuples([(depth, sample) for depth, sample, _ in rarefied_data])
    print("after index")
    
    data_array = np.array([data for _, _, data in rarefied_data])
    rarefied_df = pd.DataFrame(data_array, index=index, columns=feature_table.columns)
    #rarefied_df = pd.DataFrame((data for _, _, data in rarefied_data), index=index, columns=feature_table.columns)
    #rarefied_df = pd.DataFrame([data for _, _, data in rarefied_data], index=index, columns=feature_table.columns)
    print("after rarefied_df")
    return rarefied_df
    

def rarefy(counts, depth):
    """
    Rarefy a single sample to a specified depth by subsampling without replacement.
    
    Parameters:
    counts (array-like): A 1D array of counts for each species/feature in a sample.
    depth (int): The number of reads to subsample (rarefaction depth).
    
    Returns:
    array-like: The rarefied counts for each species.
    """
    if sum(counts) < depth:
        raise ValueError(f"Sample has fewer reads ({sum(counts)}) than the rarefaction depth ({depth}).")
    
    # Generate the rarefied sample by randomly subsampling without replacement
    counts = counts.astype(int)
    reads = np.repeat(np.arange(len(counts)), counts)  # Create a list of read indices based on counts
    subsampled_reads = np.random.choice(reads, size=depth, replace=False)  # Subsample without replacement
    rarefied_counts = np.bincount(subsampled_reads, minlength=len(counts))  # Count how many times each species was selected
    
    return rarefied_counts



#my automated rarefaction depth function
def automated_rarefaction_depth(outpur_dir: str, table: biom.Table, phylogeny: NewickFormat = None, metrics: set = None,
                                metadata: qiime2.Metadata = None, iterations: int = 10, p_samples: float = 0.8) -> None:
    
    min_depth = 1
    steps = 20
    #calculate the max reads that were used
    table_df = table.view(pd.DataFrame)
    num_samples = len(table_df)
    #only take the first half of the samples bc df is too big
    #table_df =  table_df.iloc[:(num_samples // 8)]
    reads_per_sample = table_df.sum(axis=1) 
    max_depth =  int(reads_per_sample.max())
    print("shape of table_df:", table_df.shape)

    sorted_depths = reads_per_sample.sort_values()

    # calculate the index for the p_sample loss (as an integer)
    sample_loss_index = int(np.ceil((1-p_samples) * len(sorted_depths))) - 1 #len(sorted_depths)=number of samples

    # get the sequencing depth at the sample loss index (p_samples of all samples will have fewer reads)
    # 1 - depth_threshold is the range of possible & valid values for the knee
    depth_threshold = int(sorted_depths.iloc[sample_loss_index])
    print("depth threshold:", depth_threshold)
    
    #the points that I will look at & calculate everything for
    depth_range = np.linspace(min_depth, depth_threshold, num=steps, dtype=int) 
    max_range = np.linspace(min_depth, max_depth, num=steps, dtype=int)
    #max_range = np.linspace(min_depth, 0.6*max_depth, num=steps, dtype=int)
    print("depth range:", depth_range)
    print("max range:", max_range)
    #depth_range = max_range
    spacing = depth_range[1] - depth_range[0]
    # make a new filtered df to use for finding the knee, including subsampling at the specified depths
    rarefied_df = subsample_feature_table(table_df, depth_range)

    print("rarefied_df:")
    print(rarefied_df)
    print("rarefied_df_index:")
    print(rarefied_df.index)
    print(rarefied_df.shape)  
    print(rarefied_df.columns)  

    #only for dev -> delete later
    plt.figure(figsize=(8, 6))

    #finding the knee for each sample (curve)
    located_points = [None] * num_samples
    i = 0
    d_range = depth_range
    #using the gradient method
    for sample in table_df.index:
        depth_range = d_range
        sample_rarefied_data = rarefied_df.xs(sample, level=1, axis=0)
        #print("sample_rarefied_data: ", sample_rarefied_data)
        
        # calculate how many different species have been seen at each depth
        curr_array = ((sample_rarefied_data != 0) & ~np.isnan(sample_rarefied_data)).sum(axis=1).values
        #print("curr_array:", curr_array)
        #deleting all zero values to ensure the gradient calculation works
        non_zero_array = curr_array[curr_array != 0]
        #print("non_zero_array:", non_zero_array)
        
        # Calculate the difference in lengths, if different, truncate the longer array so lengths match
        length_difference = len(curr_array) - len(non_zero_array)
        if length_difference > 0:
            depth_range = depth_range[:-length_difference]  # Remove the last 'length_difference' elements
        curr_array = non_zero_array

        plt.plot(depth_range, non_zero_array, marker='o', linestyle='-', label=sample)
        
        
        first_derivative = np.gradient(curr_array, spacing)
        second_derivative = np.gradient(first_derivative, spacing)
        max_index = np.argmax(second_derivative)
        corresponding_depth = depth_range[max_index]
        #print("Maximum second derivative:", second_derivative[max_index])
        #print("Corresponding depth value:", corresponding_depth)

        #p = curr_array[np.argmax(second_derivative)]
        located_points[i] = corresponding_depth
        i += 1
        if(i % 100 == 0):
            print("Processed", i, "samples.")

    print("located depths:", located_points)
    plt.xlabel('Sequencing Depth')
    plt.ylabel('Observed Features')
    plt.title('Rarefaction Curve for Multiple Samples')
    plt.legend()
    plt.savefig('rarefaction_curves_gradient.png')

    #taking the mean of all located_points to get the average knee point
    knee_point = np.mean(located_points)
    print("Knee_point:", knee_point)

    #finding +-5% points
    #finding out what percentile my knee point is at
    #print("sorted_depths:", sorted_depths)
    index = np.searchsorted(sorted_depths, knee_point)
    percentile = (index / len(sorted_depths)) * 100

    #calculate the +-5% percentiles & ensuring it won't be out of bounds
    lower_percentile = max(percentile - 5, 0) 
    upper_percentile = min(percentile + 5, 100) 

    #convert these percentiles to indices for the sorted_depths array
    lower_index = int((lower_percentile / 100) * len(sorted_depths))
    upper_index = min(int((upper_percentile / 100) * len(sorted_depths)), len(sorted_depths) - 1)

    #get the values corresponding to these indices
    lower_value = sorted_depths.iloc[lower_index]
    upper_value = sorted_depths.iloc[upper_index]

    print(f"The target value {knee_point} is at the {percentile:.2f}% percentile.")
    print(f"The value at {percentile-5:.2f}% is approximately {lower_value}.")
    print(f"The value at {percentile+5:.2f}% is approximately {upper_value}.")

    
    #plotting the rarefaction curve including the determined depth etc
    f_table = _compute_rarefaction_data(table, min_depth, max_depth, steps, iterations, phylogeny, metrics)


#to test & get outputs -> delete in the end
#feature_table_path = "../../table.qza"
#other feature table -> a lot bigger
feature_table_path = "../../feature-table.qza"
ft_artifact = qiime2.Artifact.load(feature_table_path)
automated_rarefaction_depth("../../", ft_artifact)







#example visualizer copied from q2-diversity/q2_diversity/_alpha/_visualizer.py
def _compute_rarefaction_data(feature_table, min_depth, max_depth, steps,
                              iterations, phylogeny, metrics):
    depth_range = np.linspace(min_depth, max_depth, num=steps, dtype=int)
    iter_range = range(1, iterations + 1)

    rows = feature_table.ids(axis='sample')
    cols = pd.MultiIndex.from_product(
        [list(depth_range), list(iter_range)],
        names=['_alpha_rarefaction_depth_column_', 'iter'])
    data = {k: pd.DataFrame(np.NaN, index=rows, columns=cols)
            for k in metrics}

    with qiime2.sdk.Context() as scope:
        feature_table = scope.ctx.make_artifact(
                'FeatureTable[Frequency]', feature_table)

        if phylogeny:
            phylogeny = scope.ctx.make_artifact('Phylogeny[Rooted]', phylogeny)

        for depth, i in itertools.product(depth_range, iter_range):
            rarefy_method = scope.ctx.get_action('feature_table', 'rarefy')
            rt, = rarefy_method(feature_table, depth)

            for metric in metrics:
                if metric in (METRICS['PHYLO']['IMPL'] |
                              METRICS['PHYLO']['UNIMPL']):
                    alpha_phylo = scope.ctx.get_action('diversity',
                                                       'alpha_phylogenetic')
                    vector, = alpha_phylo(table=rt, metric=metric,
                                          phylogeny=phylogeny)
                else:
                    alpha = scope.ctx.get_action('diversity', 'alpha')
                    vector, = alpha(table=rt, metric=metric)

                vector = vector.view(pd.Series)
                data[metric][(depth, i)] = vector
        return data


def alpha_rarefaction(output_dir: str, table: biom.Table, max_depth: int,
                      phylogeny: NewickFormat = None, metrics: set = None,
                      metadata: qiime2.Metadata = None, min_depth: int = 1,
                      steps: int = 10, iterations: int = 10) -> None:

    if metrics is None:
        metrics = {'observed_features', 'shannon'}
        if phylogeny is not None:
            metrics.add('faith_pd')
    elif not metrics:
        raise ValueError('`metrics` was given an empty set.')
    else:
        phylo_overlap = ((METRICS['PHYLO']['IMPL'] |
                          METRICS['PHYLO']['UNIMPL']) &
                         metrics)
        if phylo_overlap and phylogeny is None:
            raise ValueError('Phylogenetic metric %s was requested but '
                             'phylogeny was not provided.' % phylo_overlap)

    if max_depth <= min_depth:
        raise ValueError('Provided max_depth of %d must be greater than '
                         'provided min_depth of %d.' % (max_depth, min_depth))
    possible_steps = max_depth - min_depth
    if possible_steps < steps:
        raise ValueError('Provided number of steps (%d) is greater than the '
                         'steps possible between min_depth and '
                         'max_depth (%d).' % (steps, possible_steps))
    if table.is_empty():
        raise ValueError('Provided table is empty.')
    max_frequency = max(table.sum(axis='sample'))
    if max_frequency < max_depth:
        raise ValueError('Provided max_depth of %d is greater than '
                         'the maximum sample total frequency of the '
                         'feature_table (%d).' % (max_depth, max_frequency))

    if metadata is None:
        columns, filtered_columns = set(), set()
    else:
        # Filter metadata to only include sample IDs present in the feature
        # table. Also ensures every feature table sample ID is present in the
        # metadata.
        metadata = metadata.filter_ids(table.ids(axis='sample'))

        # Drop metadata columns that aren't categorical, or consist solely of
        # missing values.
        pre_filtered_cols = set(metadata.columns)
        metadata = metadata.filter_columns(column_type='categorical',
                                           drop_all_missing=True)
        filtered_columns = pre_filtered_cols - set(metadata.columns)

        metadata_df = metadata.to_dataframe()
        if metadata_df.empty or len(metadata.columns) == 0:
            raise ValueError("All metadata filtered after dropping columns "
                             "that contained non-categorical data.")
        metadata_df.columns = pd.MultiIndex.from_tuples(
            [(c, '') for c in metadata_df.columns],
            names=('_alpha_rarefaction_depth_column_', 'iter'))
        columns = metadata_df.columns.get_level_values(0)
    data = _compute_rarefaction_data(table, min_depth, max_depth,
                                     steps, iterations, phylogeny, metrics)

    filenames = []
    for m, data in data.items():
        metric_name = quote(m)
        filename = '%s.csv' % metric_name

        if metadata is None:
            n_df = _compute_summary(data, 'sample-id')
            jsonp_filename = '%s.jsonp' % metric_name
            _alpha_rarefaction_jsonp(output_dir, jsonp_filename, metric_name,
                                     n_df, '')
            filenames.append(jsonp_filename)
        else:
            merged = data.join(metadata_df, how='left')
            for column in columns:
                column_name = quote(column)
                reindexed_df, counts = _reindex_with_metadata(column,
                                                              columns,
                                                              merged)
                c_df = _compute_summary(reindexed_df, column, counts=counts)
                jsonp_filename = "%s-%s.jsonp" % (metric_name, column_name)
                _alpha_rarefaction_jsonp(output_dir, jsonp_filename,
                                         metric_name, c_df, column)
                filenames.append(jsonp_filename)

        with open(os.path.join(output_dir, filename), 'w') as fh:
            data.columns = [
                'depth-%d_iter-%d' % (t[0], t[1])
                for t in data.columns.values]
            if metadata is not None:
                data = data.join(metadata.to_dataframe(), how='left')
            data.to_csv(fh, index_label=['sample-id'])

    index = os.path.join(TEMPLATES, 'alpha_rarefaction_assets', 'index.html')
    q2templates.render(index, output_dir,
                       context={'metrics': list(metrics),
                                'filenames': [quote(f) for f in filenames],
                                'columns': list(columns),
                                'steps': steps,
                                'filtered_columns': sorted(filtered_columns)})

    shutil.copytree(os.path.join(TEMPLATES, 'alpha_rarefaction_assets',
                                 'dist'),
                    os.path.join(output_dir, 'dist'))


alpha_rarefaction_unsupported_metrics = {'osd', 'lladser_ci', 'strong',
                                         'esty_ci', 'kempton_taylor_q',
                                         'chao1_ci'}