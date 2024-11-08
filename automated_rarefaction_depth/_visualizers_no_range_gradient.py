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
    rarefied_df = pd.DataFrame([data for _, _, data in rarefied_data], index=index, columns=feature_table.columns)

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
    steps = 10
    #calculate the max reads that were used
    table_df = table.view(pd.DataFrame)
    num_samples = len(table_df)
    reads_per_sample = table_df.sum(axis=1) 
    print(table_df)
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
                rarefied_sample = rarefy(table_df.loc[sample].values, i)
                c = np.count_nonzero(rarefied_sample)
                array_sample[j, i] = c
                
        
        array_sample_avg = np.mean(array_sample, axis=0)

        if (s < 44 and s > 34): #plot only the first 5 samples
            plt.plot(max_range, array_sample_avg, marker='o', linestyle='-', label=sample)
            #print("array sample:", array_sample)
            #print("array sample avg:", array_sample_avg)
        
        #using the gradient method
        curr_array = array_sample_avg

        first_derivative = np.gradient(curr_array, max_range)
        second_derivative = np.gradient(first_derivative, max_range)
        max_index = np.argmax(second_derivative)
        corresponding_depth = max_range[max_index]
        knee_points[s] = corresponding_depth
       
        if(s % 200 == 0):
            print("Processed", s, "samples.")
        s += 1

    #just for visualizations of a single sample during development -> delete later!!
    plt.xlabel('Sequencing Depth')
    plt.ylabel('Observed Features')
    plt.title('Rarefaction Curve')
    plt.legend()
    plt.savefig('example_curve_no_range_gradient.png')

    #calculate the average of all knee points
    knee_point_avg = round(np.mean(knee_points))
    print("knee points:", knee_points)

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
    print(f"The value at {percentile-5:.2f}% is approximately {lower_value}.")
    print(f"The value at {percentile+5:.2f}% is approximately {upper_value}.")

    

#to test & get outputs -> delete in the end
#feature_table_path = "../../table.qza"
#other feature tables
feature_table_path = "../../atacama_soil_table.qza"
#feature_table_path = "../../parkinson_mouse_dada2_table.qza"
#very big one
#feature_table_path = "../../feature-table.qza"
ft_artifact = qiime2.Artifact.load(feature_table_path)
automated_rarefaction_depth("../../", ft_artifact)


