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
                # If a sample has fewer reads than the rarefaction depth, set it to NaN or 0
                rarefied_data.append((depth, sample, np.nan * np.ones(feature_table.shape[1])))

        
    index = pd.MultiIndex.from_tuples([(depth, sample) for depth, sample, _ in rarefied_data])
    rarefied_df = pd.DataFrame([data for _, _, data in rarefied_data], index=index, columns=feature_table.columns)

    return rarefied_df
    # Combine results into a single DataFrame
    # Using list comprehension
    """return pd.concat(
        {depth: pd.DataFrame(rarefied_tables[depth], index=feature_table.index, columns=feature_table.columns)
         for depth in depths}, axis=1)"""

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
    reads_per_sample = table_df.sum(axis=1) 
    #print(table_df)
    #print(reads_per_sample)
    #print(reads_per_sample.max())
    max_depth =  int(reads_per_sample.max())

    sorted_depths = reads_per_sample.sort_values()

    # calculate the index for the p_sample loss (as an integer)
    sample_loss_index = int(np.ceil((1-p_samples) * len(sorted_depths))) - 1 #len(sorted_depths)=number of samples

    # get the sequencing depth at the sample loss index (p_samples of all samples will have fewer reads)
    # 1 - depth_threshold is the range of possible & valid values for the knee
    depth_threshold = int(sorted_depths.iloc[sample_loss_index])
    print("depth threshold:", depth_threshold)
    
    #the points that I will look at & calculate everything for
    depth_range = np.linspace(min_depth, depth_threshold, num=steps, dtype=int) 
    print("depth range:", depth_range)

    # make a new filtered df to use for finding the knee, including subsampling at the specified depths
    rarefied_df = subsample_feature_table(table_df, depth_range)
    print("rarefied_df:")
    print(rarefied_df)
    print("rarefied_df_index:")
    print(rarefied_df.index)
    print(rarefied_df.shape)  # Number of rows and columns
    print(rarefied_df.columns)  # Column names

    """# filter samples by depth_threshold, only keep the ones <= 
    filtered_samples = reads_per_sample[reads_per_sample <= depth_threshold].index
    print("DataFrame Columns:", table_df.columns)
    print("Filtered Samples:", filtered_samples.columns)
    filtered_table = table_df[filtered_samples]

    # function to subsample a sample down to the depth_threshold, incl several iterations for statistical validity 
    def subsample_sample(sample_series, depth):
        subsample_counts = np.zeros(sample_series.shape)  
        for i in range(iterations):
            subsample = sample_series.sample(n=depth, replace=False)
            subsample_counts += subsample.reindex(sample_series.index, fill_value=0)
        return (subsample_counts / iterations)

    # Subsample samples with more than depth_threshold reads
    # each row is a sample, each column is a different feature/species
    subsam_t = filtered_table.subsample(depth_threshold, axis=1, by_id=False,
                            with_replacement=False)"""
    
    num_samples = len(table_df)

    #calculating the knee point
    knee_points = {}
    #only for dev -> delete later
    plt.figure(figsize=(8, 6))

    for sample in table_df.index:
        # Extract rarefied counts for this sample across all depths
        sample_rarefied_data = rarefied_df.xs(sample, level=1, axis=0)
        print("sample_rarefied_data: ", sample_rarefied_data)
        
        # calculate how many different species have been seen at each depth
        num_different_features = (sample_rarefied_data != 0).sum(axis=1).values
        print("num_different_features", num_different_features)

        #just for visualizations of a single sample during development -> delete later!!
        #plt.figure(figsize=(8, 6))
        plt.plot(depth_range, num_different_features, marker='o', linestyle='-', label=sample)
        """plt.xlabel('Sequencing Depth')
        plt.ylabel('Observed Features')
        plt.title('Rarefaction Curve')
        plt.legend()
        plt.savefig('example_curve.png')"""
        
        # Use KneeLocator to find the knee point (depth where total abundance starts leveling off)
        kneedle = KneeLocator(depth_range, num_different_features, curve="concave", direction="increasing")
        
        # Store the knee point for this sample
        knee_points[sample] = kneedle.knee

    #just for dev ->delete later    
    plt.xlabel('Sequencing Depth')
    plt.ylabel('Observed Features')
    plt.title('Rarefaction Curve for Multiple Samples')
    plt.legend()
    plt.savefig('rarefaction_curves_overlay.png')

    print("Knee_points:", knee_points)
    #using kneedle algorithm/ Knee_locator function
    #finding the knee for each sample (curve)
    """ for i in range(num_samples):
        #using the kneeLocator object to find the knee of the curve
        #-> change x & y values!!! supposed to be arrays
        x = 1    # number of reads (probably just index of the y-value)
        y =  1   # sum of all different features so far
        knee_locator = KneeLocator(x, y, curve='concave', direction='increasing')
        located_points[i] = knee_locator.elbow #gives the x-value of the proposed knee point"""

    """#using the gradient method
    for i in range(num_samples):
        curr_array = subsampled_table.iloc[:, i].to_numpy() 
        first_derivative = np.gradient(curr_array)
        second_derivative = np.gradient(first_derivative)
        located_points[i] = subsampled_table[np.argmax(second_derivative)]"""

    #taking the mean of all located_points to get the average knee point
    knee_point = sum(knee_points.values()) / len(knee_points.values()) 
    print("knee_point (avg)", knee_point)

    #plot the different knee points for each sample, just for dev -> delete later!!
    keys = list(knee_points.keys())
    values = list(knee_points.values())
    plt.figure(figsize=(8, 6))
    plt.subplot(1, 2, 2)  # 1 row, 2 columns, 2nd subplot
    plt.bar(keys, values, color='orange')
    plt.title('Bar Plot of Dictionary Values')
    plt.xlabel('Keys')
    plt.ylabel('Values')
    plt.grid(axis='y')
    plt.tight_layout()
    plt.savefig('knee_points.png')

    
    #plotting the rarefaction curve including the determined depth etc
    f_table = _compute_rarefaction_data(table, min_depth, max_depth, steps, iterations, phylogeny, metrics)


#to test & get outputs -> delete in the end
feature_table_path = "../../table.qza"
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