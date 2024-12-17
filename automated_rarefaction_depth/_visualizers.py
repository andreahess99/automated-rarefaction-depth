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
import q2templates
import tempfile
import biom
import itertools
from qiime2 import Artifact


from kneed import KneeLocator
import altair as alt
import time
import psutil
import tracemalloc



#helper function
def rarefy(counts, depth, iteration, seed):
    if sum(counts) < depth:
        raise ValueError(f"Sample has fewer reads ({sum(counts)}) than the rarefaction depth ({depth}).")
    
    # Generate the rarefied sample by randomly subsampling without replacement
    np.random.seed(iteration+seed)
    counts = counts.astype(int)
    reads = np.repeat(np.arange(len(counts)), counts)  # Create a list of read indices based on counts
    subsampled_reads = np.random.choice(reads, size=depth, replace=False)  # Subsample without replacement
    rarefied_counts = np.bincount(subsampled_reads, minlength=len(counts))  # Count how many times each species was selected
    
    return rarefied_counts



#my automated rarefaction depth function
def rarefaction_depth(outpur_dir: str, table: biom.Table, seed: int = 42,
                                iterations: int = 10, table_size: int = None, steps: int = 20,
                                percent_samples: float = 0.95, algorithm: str = 'kneedle') -> None:
    
    # Measure runtime & memory usage
    start_time = time.time()
    tracemalloc.start()
    
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

    #checks to make sure everything has the format it should have
    if table_df.empty:
        raise ValueError("The feature table is empty.")
    if not np.issubdtype(table_df.values.dtype, np.number):
        raise ValueError("The feature table contains non-numerical values.")
    
    
    

    #table_size = 200
    #adjusting table size if it's too big -> keep table_size rows
    if (table_size is not None and len(table_df) > table_size):
        table_df = table_df.sample(n=table_size, random_state=seed) #randomly select rows, set random_state for reproducibility
        table_df = table_df.loc[:, ~(table_df.isna() | (table_df == 0)).all(axis=0)] #remove columns where all values are either 0 or NaN

    num_samples = len(table_df)
    reads_per_sample = table_df.sum(axis=1) 
    max_reads = reads_per_sample.max()
    sorted_depths = reads_per_sample.sort_values()

    # calculate the index for the p_sample loss (as an integer)
    sample_loss_index = int(np.ceil((1-percent_samples) * len(sorted_depths))) - 1 #len(sorted_depths)=number of samples

    # get the sequencing depth at the sample loss index (p_samples of all samples will have fewer reads)
    # 1 - depth_threshold is the range of possible & valid values for the knee
    depth_threshold = int(sorted_depths.iloc[sample_loss_index])

    #go through every sample and calculate the observed features for each depth
    knee_points = [None] * num_samples
    s = 0
    #plt.figure(figsize=(8, 6))
    df_list = []
    
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
        
        sample_df = pd.DataFrame({'depth': max_range, 'observed_features': array_sample_avg, 'sample': sample})
        df_list.append(sample_df)
        
        if(algorithm.lower().strip() == 'kneedle'):
            # Use KneeLocator to find the knee point (depth where total abundance starts leveling off)
            kneedle = KneeLocator(max_range, array_sample_avg, curve="concave", direction="increasing")
            knee_points[s] = kneedle.knee
        else:
            #using the gradient method
            curr_array = array_sample_avg
            first_derivative = np.gradient(curr_array, max_range)
            second_derivative = np.gradient(first_derivative, max_range)
            max_index = np.argmax(second_derivative)
            knee_points[s] = max_range[max_index]

        s += 1
        

    combined_df = pd.concat(df_list, ignore_index=True)

    # Filter out None values
    knee_points_filtered = [point for point in knee_points if point is not None]
    #calculate the average of all knee points
    knee_point_avg = round(np.mean(knee_points_filtered))
    knee_point_median = round(np.median(knee_points_filtered))
    
    #just for visualizations during development -> delete later!!
    """plt.axvspan(0, depth_threshold, color='peachpuff', alpha=0.75, label='Acceptable Range')
    plt.axvline(x=knee_point_avg, color='red', linestyle='--', label='Vertical Line at Knee Point (avg)')
    plt.axvline(x=knee_point_median, color='black', linestyle='--', label='Vertical Line at Knee Point (median)')
    plt.xlabel('Sequencing Depth')
    plt.ylabel('Observed Features')
    plt.title('Rarefaction Curve')
    plt.legend()
    plt.savefig('example_curve_visualizations.png')"""


    knee_point = knee_point_avg

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


    #plotting with altair
    #plotting the rarefaction curves including the shaded area
    depth_lines = pd.DataFrame({
        'x': [lower_value, knee_point, upper_value],  # x-values for vertical lines
        'label': ['-5%', 'Knee point', '+5%']  
    })

    reads_per_sample_df = reads_per_sample.reset_index()
    reads_per_sample_df.columns = ['sample', 'reads_per_sample']

    zoom = alt.selection_interval(bind='scales')

    param_checkbox = alt.param(bind=alt.binding_checkbox(name='Show depth specified by slider as a line on the plot: '))

    s = alt.param(
        name='position', bind=alt.binding_range(min=0, max=max_reads, step=50, name='Rarefaction Depth Line'), value=knee_point
    )

    chart = alt.Chart(combined_df).mark_line(point=True).encode(
        x=alt.X('depth:Q', title='Read Depth'),
        y=alt.Y('observed_features:Q', title='# Observed Features'),
        color=alt.Color('sample:N', legend=None).scale(scheme='category10') #remove .scale(..) if I want the og colors
    ).add_params(
        zoom,
        param_checkbox
    ).properties(
        title='Rarefaction Curves'
    )
    shaded_area = alt.Chart(pd.DataFrame({
        'x_min': [0],
        'x_max': [depth_threshold]
    })).mark_rect(opacity=0.7, color='peachpuff').encode(
        x='x_min:Q',
        x2='x_max:Q'
    )
    my_colors = ["#377eb8", "#f781bf", "#984ea3"]  # possibly change purple to red?, 984ea3 purple, 17becf cyan, 377eb8 blue, f781bf pink, e41a1c red

    vertical_lines = alt.Chart(depth_lines).mark_rule(strokeWidth=2).encode(
        x='x:Q',
        color=alt.Color('label:N', legend=alt.Legend(title="Thresholds"), scale=alt.Scale(range=my_colors))
    )

    final_chart = alt.layer(
        shaded_area, chart, vertical_lines
    ).resolve_scale(
        x='shared',
        color='independent'
    ).properties(
        title='Rarefaction Curves',
        width=500,
        height=350
    )
    
    vertical_line = alt.Chart(pd.DataFrame({'position': [0]})).mark_rule(color='black', strokeWidth=2).encode(
        x='position:Q'
    ).add_params(
        s
    ).transform_calculate(
        position='position'
    ).transform_filter(
        param_checkbox
    ).properties(
        title='Interactive Vertical Line',
        width=500, #600
        height=350 #400
    )
    
    final_with_line = alt.layer(final_chart, vertical_line).resolve_scale( 
        x='shared',
        y='shared'
    )

    text = [f'Knee point: {knee_point} (at the {percentile:.2f} percentile)',
            f'Knee point -5%: {lower_value:.0f} (at the {lower_percentile:.2f} percentile)',
            f'Knee point +5%: {upper_value:.0f} (at the {upper_percentile:.2f} percentile)',
            f' ']
          
    if (knee_point > depth_threshold):
        text.append(f"The knee point is above the specified, acceptable depth range!")
        text.append(f"The upper bound of the acceptable range is {depth_threshold}.")
        text.append(f" ")
        text.append(f"If the calculated knee point is used, {percentile:.2f}% of the samples will")
        text.append(f" be excluded because they have too little reads.")
        
    text_l = pd.DataFrame({'text': ['\n'.join(text)]})
    text_lines = alt.Chart(text_l).mark_text(size=14, align='left', baseline='top', lineBreak="\n", dx=-50).encode(text='text:N').properties(width=100, height=50)#(width=100, height=300)
    upper_chart = alt.vconcat(final_with_line, text_lines).properties(spacing=0)
    #final_chart.save('rarefaction_curves.html')


    #barplot of reads_per_sample
    predicate = alt.datum.reads_per_sample >= s 
  
    barplot = alt.Chart(reads_per_sample_df).mark_bar().encode(
        x=alt.X('reads_per_sample:Q', bin=alt.Bin(maxbins=50), title='Reads per Sample'),
        y=alt.Y('count()', title='# Samples'),
        color=alt.value('steelblue')  
    ).add_params(
        s,
        zoom
    ).transform_filter(
        predicate
    ).transform_calculate(
        position='position'
    ).properties(
        title='Histogram of Reads per Sample',
        width=500,
        height=350
    )

    background = alt.Chart(reads_per_sample_df).mark_bar().encode(
        x=alt.X('reads_per_sample:Q', bin=alt.Bin(maxbins=50), title='Reads per Sample'),  
        y=alt.Y('count()', title='# Samples'),
        color=alt.value('lightgrey')
    ).add_params(
        zoom
    ).properties(
        title='Histogram of Reads per Sample',
        width=500,
        height=350
    )
    barplot_combined = alt.layer(background, barplot).resolve_scale(
        x='shared',
        y='shared'
    )
        
    combined_chart = alt.hconcat(upper_chart, barplot_combined).properties(spacing=60)
    #combined_chart.save('combined_chart.html')
    combined_chart.save('combined_chart.html')



    #end measuring runtime & memory usage
    current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    end_time = time.time()
    
    print(f"Runtime: {end_time - start_time:.4f} seconds")
    print(f"Current memory usage: {current / 10**6:.4f} MB")
    print(f"Peak memory usage: {peak / 10**6:.4f} MB")

    # copied from https://github.com/bokulich-lab/q2-moshpit/blob/ea8cb818462a098651169f2c884b9005f76d75fe/q2_moshpit/busco/busco.py
    #do I need this??
    # Render
    tabbed_context = {}
    #vega_json = json.dumps(combined_chart.to_dict())
    vega_json = combined_chart.to_json()
    temp_dir = tempfile.TemporaryDirectory()
    output_dir = os.path.join(temp_dir.name, "rarefaction-depth")
    tabbed_context.update({
        "tabs":[{"title": "Rarefaction Depth Overview", "url": "index.html"}],
        "vega_json": vega_json,
        "page_size": 100
    })
    templates = os.path.join(os.path.dirname(__file__), 'assets')
    q2templates.render(templates, output_dir, context=tabbed_context)

    

#to test & get outputs -> delete in the end
feature_table_path = "../../table.qza"
#other feature tables
#feature_table_path = "../../atacama_soil_table.qza"
#feature_table_path = "../../parkinson_mouse_dada2_table.qza"
#very big one
#feature_table_path = "../../feature-table.qza"
ft_artifact = qiime2.Artifact.load(feature_table_path)
rarefaction_depth("../../", ft_artifact)
