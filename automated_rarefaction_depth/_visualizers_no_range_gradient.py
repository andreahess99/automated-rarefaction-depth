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
import time
import psutil
import tracemalloc
#from qiime2.plugins.diversity.methods import alpha
#from . import METRICS
from q2_types.tree import NewickFormat


#helper functions
def rarefy(counts, depth, iter, seed):
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
    #set a seed for reproducibility
    np.random.seed(iter+seed)
    counts = counts.astype(int)
    reads = np.repeat(np.arange(len(counts)), counts)  # Create a list of read indices based on counts
    subsampled_reads = np.random.choice(reads, size=depth, replace=False)  # Subsample without replacement
    rarefied_counts = np.bincount(subsampled_reads, minlength=len(counts))  # Count how many times each species was selected
    
    return rarefied_counts



#my automated rarefaction depth function
def automated_rarefaction_depth(outpur_dir: str, table: biom.Table, phylogeny: NewickFormat = None, seed: int = 42,
                                metadata: qiime2.Metadata = None, iterations: int = 10, table_size: int = None,
                                steps: int = 20, p_samples: float = 0.8, algorithm: str = 'kneedle') -> None:

    # Measure runtime & memory usage
    start_time = time.time()
    tracemalloc.start() 

    min_depth = 1
    #calculate the max reads that were used
    table_df = table.view(pd.DataFrame)

    #make a new, smaller df where all the max depths are similar & save it as a new qza file
    #somehow doesn't work -> no clue why
    """alt_table = table_df.sum(axis=1).sort_values().iloc[5300:5500]
    filtered_df = table_df.loc[alt_table.index]
    filtered_df = filtered_df.loc[:, ~(filtered_df.isna() | (filtered_df == 0)).all(axis=0)]
    filtered_df = filtered_df.loc[~(filtered_df.sum(axis=1) == 0)]  
    filtered_df.columns = filtered_df.columns.astype(str).str.strip()
    if filtered_df.columns.duplicated().any():
        print("Duplicate observation IDs found:", filtered_df.columns[filtered_df.columns.duplicated()])
    filtered_df = filtered_df.T.groupby(filtered_df.columns).sum()
    ind = filtered_df.index.tolist()
    col = filtered_df.columns.tolist()
    biom_table = biom.Table(filtered_df.values, 
                        observation_ids=[str(x) for x in col], 
                        sample_ids=[str(x) for x in ind])
    feature_table = qiime2.Artifact.import_data('FeatureTable[Frequency]', biom_table)
    feature_table.save('filtered_cancer_microbiome_table.qza')
    print("done")"""

    #for filtering the cancer microbiome dataset to a smaller size
    """alt_table = table_df.sum(axis=1).sort_values().iloc[6950:7000]
    filtered_df = table_df.loc[alt_table.index]
    table_df = filtered_df
"""
    #for filtering the moving pictures tuotrial dataset to a smaller size
    """alt_table = table_df.sum(axis=1).sort_values().iloc[22:33]
    filtered_df = table_df.loc[alt_table.index]
    table_df = filtered_df"""
    

    if table_df.empty:
        raise ValueError("The feature table is empty.")
    table_size = 200
    #adjusting table size if it's too big -> keep table_size rows
    if (table_size is not None and len(table_df) > table_size):
        table_df = table_df.sample(n=table_size, random_state=seed)  #randomly select rows, set random_state for reproducibility
        table_df = table_df.loc[:, ~(table_df.isna() | (table_df == 0)).all(axis=0)] #remove columns where all values are either 0 or NaN

    """if(len(table_df) > 100):
        #delete the top 5 and bottom 5 samples
        table_df['row_sum'] = table_df.sum(axis=1)
        df_sorted = table_df.sort_values(by='row_sum')
        df_filtered = df_sorted.iloc[5:-5]
        df_filtered = df_filtered.drop(columns='row_sum')
        table_df = df_filtered"""

    num_samples = len(table_df)
    reads_per_sample = table_df.sum(axis=1) 
    max_depth =  int(reads_per_sample.max())
    sorted_depths = reads_per_sample.sort_values()
    print("sorted depths:", sorted_depths)

    # calculate the index for the p_sample loss (as an integer)
    sample_loss_index = int(np.ceil((1-p_samples) * len(sorted_depths))) - 1 #len(sorted_depths)=number of samples

    # get the sequencing depth at the sample loss index (p_samples of all samples will have fewer reads)
    # 1 - depth_threshold is the range of possible & valid values for the knee
    depth_threshold = int(sorted_depths.iloc[sample_loss_index])
    print("depth threshold:", depth_threshold)
  

    #go through every sample and calculate the observed features for each depth
    knee_points = [None] * num_samples
    s = 0
    plt.figure(figsize=(8, 6))
    

    for sample in table_df.index:
        max_range = np.linspace(min_depth, reads_per_sample.loc[sample], num=steps, dtype=int)
        #array to store the observed features for each depth, rows: iterations, columns: depths
        #the different iterations are used to calculate the average of the observed features for each depth & have statistical validity
        array_sample = np.empty((iterations, steps))

        for i in range(steps):
            for j in range(iterations):
                rarefied_sample = rarefy(table_df.loc[sample].values, max_range[i], j, seed) 
                c = np.count_nonzero(rarefied_sample)
                array_sample[j, i] = c
                
        
        array_sample_avg = np.mean(array_sample, axis=0)
        #array_sample_median = np.median(array_sample, axis=0)

        
        if (s < 50 and s > 4): #plot only the first 5 samples
            plt.plot(max_range, array_sample_avg, marker='o', linestyle='-', label=sample)
            #print("array sample:", array_sample)
            #print("array sample avg:", array_sample_avg)
        
        #using the gradient method
        curr_array = array_sample_avg
        first_derivative = np.gradient(curr_array, max_range)
        second_derivative = np.gradient(first_derivative, max_range)
        max_index = np.argmax(second_derivative)
        knee_points[s] = max_range[max_index]

        #ratio of the slopes method
        """slopes = np.gradient(curr_array, max_range)
        ratios = slopes[:-1] / slopes[1:]
        knee_point_index = np.argmin(ratios)  # or where the ratio falls below a chosen threshold
        knee_points[s] = max_range[knee_point_index]"""
       
        if(s % 200 == 0):
            print("Processed", s, "samples.")
        s += 1


    #calculate the average of all knee points
    knee_point_avg = round(np.mean(knee_points))
    knee_point_median = round(np.median(knee_points))
    print("knee point avg:", knee_point_avg)
    print("knee point median:", knee_point_median)

    #just for visualizations of a single sample during development -> delete later!!
    plt.axvspan(0, depth_threshold, color='peachpuff', alpha=0.75, label='Acceptable Range')
    plt.axvline(x=knee_point_avg, color='red', linestyle='--', label='Vertical Line at Knee Point (avg)')
    plt.axvline(x=knee_point_median, color='black', linestyle='--', label='Vertical Line at Knee Point (median)')
    plt.xlabel('Sequencing Depth')
    plt.ylabel('Observed Features')
    plt.title('Rarefaction Curve')
    plt.legend()
    plt.savefig('example_curve_no_range_gradient.png')


    knee_point = knee_point_avg
    perc_avg = (np.searchsorted(sorted_depths, knee_point_avg) / len(sorted_depths)) * 100
    #checking if the knee point is in the range of acceptable values
    if (knee_point_avg > depth_threshold):
        print("The knee point is above the depth threshold.")
        knee_point = depth_threshold
        print("actual knee point:", knee_point_avg)
        print("corrected knee point:", knee_point)


    #finding +-5% points
    #finding out what percentile my knee point is at
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

    print(f"The target value {knee_point_avg} is at the {perc_avg:.2f}% percentile.")
    print(f"The target value {knee_point} is at the {percentile:.2f}% percentile.")
    print(f"The value at {max(percentile-5, 0):.2f}% is approximately {lower_value}.")
    print(f"The value at {percentile+5:.2f}% is approximately {upper_value}.")

    #end measuring runtime & memory usage
    current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    end_time = time.time()
    print(f"Runtime: {end_time - start_time:.4f} seconds")
    print(f"Current memory usage: {current / 10**6:.4f} MB")
    print(f"Peak memory usage: {peak / 10**6:.4f} MB")

    

#to test & get outputs -> delete in the end
#feature_table_path = "../../table.qza"
#other feature tables
#feature_table_path = "../../atacama_soil_table.qza"
#feature_table_path = "../../parkinson_mouse_dada2_table.qza"
#very big one, cancer microbiome
feature_table_path = "../../feature-table.qza"
ft_artifact = qiime2.Artifact.load(feature_table_path)
automated_rarefaction_depth("../../", ft_artifact)


