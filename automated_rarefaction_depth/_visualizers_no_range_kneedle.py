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
    # Set a seed for reproducibility
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
    
    #check which algo is selected and run appropriate one -> to be implemented!!

    # Measure runtime & memory usage
    start_time = time.time()
    tracemalloc.start() 
    # Get the current process & start cpu measurement
    #process = psutil.Process(os.getpid())
    
    min_depth = 1
    table_df = table.view(pd.DataFrame)

    #for filtering the cancer microbiome dataset to a smaller size
    """alt_table = table_df.sum(axis=1).sort_values().iloc[1000:8000]
    filtered_df = table_df.loc[alt_table.index]
    table_df = filtered_df"""

    #for filtering the moving pictures tuotrial dataset to a smaller size
    """alt_table = table_df.sum(axis=1).sort_values().iloc[22:33]
    filtered_df = table_df.loc[alt_table.index]
    table_df = filtered_df"""

    if table_df.empty:
        raise ValueError("The feature table is empty.")
    table_size = 100
    #adjusting table size if it's too big -> keep table_size rows
    if (table_size is not None and len(table_df) > table_size):
        table_df = table_df.sample(n=table_size, random_state=seed) #randomly select rows, set random_state for reproducibility
        table_df = table_df.loc[:, ~(table_df.isna() | (table_df == 0)).all(axis=0)] #remove columns where all values are either 0 or NaN

    """if (len(table_df) > 100):
        #delete the top 5 and bottom 5 samples
        table_df['row_sum'] = table_df.sum(axis=1)
        df_sorted = table_df.sort_values(by='row_sum')
        df_filtered = df_sorted.iloc[5:-5]
        df_filtered = df_filtered.drop(columns='row_sum')
        table_df = df_filtered
        print("truncated")"""

    num_samples = len(table_df)
    #calculate the max reads that were used
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

        if (s < 50 and s > 10): #plot only the first 5 samples
            plt.plot(max_range, array_sample_avg, marker='o', linestyle='-', label=sample)
        
        # Use KneeLocator to find the knee point (depth where total abundance starts leveling off)
        kneedle = KneeLocator(max_range, array_sample_avg, curve="concave", direction="increasing")
        # Store the knee point 
        knee_points[s] = kneedle.knee
        s += 1
        if(s % 100 == 0):
            print("Processed", s, "samples.")


    #calculate the average of all knee points
    #print("knee points:", knee_points)
    # Filter out None values
    knee_points_filtered = [point for point in knee_points if point is not None]
    knee_point_avg = round(np.mean(knee_points_filtered))
    knee_point_median = round(np.median(knee_points_filtered))

    #print("array sample: ", array_sample)
    #just for visualizations of a single sample during development -> delete later!!
    plt.axvspan(0, depth_threshold, color='peachpuff', alpha=0.75, label='Acceptable Range')
    plt.axvline(x=knee_point_avg, color='red', linestyle='--', label='Vertical Line at Knee Point (avg)')
    plt.axvline(x=knee_point_median, color='black', linestyle='--', label='Vertical Line at Knee Point (median)')
    plt.xlabel('Sequencing Depth')
    plt.ylabel('Observed Features')
    plt.title('Rarefaction Curve')
    plt.legend()
    plt.savefig('example_curve_no_range.png')
    
    print("knee point avg:", knee_point_avg)
    print("knee point median:", knee_point_median)
    

    knee_point = knee_point_avg
    #checking if the knee point is in the range of acceptable values
    if (knee_point_avg > depth_threshold):
        print("The knee point is above the depth threshold.")
        knee_point = depth_threshold


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

    perc_avg = (np.searchsorted(sorted_depths, knee_point_avg) / len(sorted_depths)) * 100
    print(f"The target value {knee_point_avg} is at the {perc_avg:.2f}% percentile.")
    print(f"The target value {knee_point} is at the {percentile:.2f}% percentile.")
    print(f"The value at {max(percentile-5, 0):.2f}% is approximately {lower_value}.")
    print(f"The value at {percentile+5:.2f}% is approximately {upper_value}.")

    #end measuring runtime & memory usage
    current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    end_time = time.time()
    # Measure CPU usage and runtime
    #cpu_usage = process.cpu_percent(interval=3)  # CPU usage since the last call

    # Report CPU usage
    #print(f"CPU Usage: {cpu_usage}%")
    print(f"Runtime: {end_time - start_time:.4f} seconds")
    print(f"Current memory usage: {current / 10**6:.4f} MB")
    print(f"Peak memory usage: {peak / 10**6:.4f} MB")

    

#to test & get outputs -> delete in the end
#feature_table_path = "../../table.qza"
#other feature tables
#feature_table_path = "../../atacama_soil_table.qza"
#feature_table_path = "../../parkinson_mouse_dada2_table.qza"
#very big one
feature_table_path = "../../feature-table.qza"
ft_artifact = qiime2.Artifact.load(feature_table_path)
automated_rarefaction_depth("../../", ft_artifact)


