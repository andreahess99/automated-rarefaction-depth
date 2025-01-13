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
from shutil import copytree
from bs4 import BeautifulSoup



def rarefy(counts, depth, iteration, seed):
    if sum(counts) < depth:
        raise ValueError(f"Sample has fewer reads ({sum(counts)}) than the rarefaction depth ({depth}).")
    
    # Generate the rarefied sample by randomly subsampling without replacement
    np.random.seed(iteration+seed)
    counts = counts.astype(int)
    reads = np.repeat(np.arange(len(counts)), counts)  
    subsampled_reads = np.random.choice(reads, size=depth, replace=False)  
    rarefied_counts = np.bincount(subsampled_reads, minlength=len(counts)) 
    
    return rarefied_counts

def change_html_file(file_path: str) -> None:
    with open(file_path, 'r') as file:
        soup = BeautifulSoup(file, 'html.parser')
        link_tag = soup.new_tag('link', rel='stylesheet', href='./css/styles.css')#href='./automated_rarefaction_depth/assets/css/styles.css'
        soup.head.append(link_tag)

    with open(file_path, 'w') as file:
        file.write(str(soup))


def rarefaction_depth(output_dir: str, table: pd.DataFrame, seed: int = 42,
                                iterations: int = 10, table_size: int = None, steps: int = 20,
                                percent_samples: float = 0.8, algorithm: str = 'kneedle') -> None:
    
    # Measure runtime & memory usage
    start_time = time.time()
    tracemalloc.start()
    
    min_depth = 1
    #table_df = table.view(pd.DataFrame)
    table_df = table

    #checks to make sure everything has the format it should have
    if table_df.empty:
        raise ValueError("The feature table is empty.")
    if not np.issubdtype(table_df.values.dtype, np.number):
        raise ValueError("The feature table contains non-numerical values.")
    
    #table_size = 500
    #adjusting table size if it's too big -> keep table_size rows
    if (table_size is not None and len(table_df) > table_size):
        table_df = table_df.sample(n=table_size, random_state=seed) #randomly select rows, set random_state for reproducibility
        table_df = table_df.loc[:, ~(table_df.isna() | (table_df == 0)).all(axis=0)] #remove columns where all values are either 0 or NaN


    num_samples = len(table_df)
    print(f"Table size: {num_samples} samples; algorithm: {algorithm}")
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
    df_list = []
    
    for sample in table_df.index:
        max_range = np.linspace(min_depth, reads_per_sample.loc[sample], num=steps, dtype=int)
        #array to store the observed features for each depth, rows: iterations, columns: depths
        #the different iterations are used to calculate the average of the observed features for each depth & have statistical validity
        array_sample = np.empty((iterations, steps))

        sample_values = table_df.loc[sample].values	

        for i in range(steps):
            temp = max_range[i]
            for j in range(iterations):
                rarefied_sample = rarefy(sample_values, temp, j, seed)
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
    print(f"df_list size: {combined_df.shape}")

    # Filter out None values
    knee_points_filtered = [point for point in knee_points if point is not None]
    #calculate the average of all knee points
    knee_point = round(np.mean(knee_points_filtered))
    print(f"Knee point: {knee_point}")


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

    after_algo_time = time.time()

    #plotting with altair
    #registering the theme
    def boots():
        font = "system-ui, -apple-system, 'Segoe UI', Roboto, 'Helvetica Neue', Arial, 'Noto Sans', 'Liberation Sans', sans-serif"
        return {
        "config" : {
             "title": {'font': font},
             "axis": {
                  "labelFont": font,
                  "titleFont": font,
                  "labelFontSize": 14,
                  "titleFontSize": 13
             },
             "header": {
                  "labelFont": font,
                  "titleFont": font,
                  "labelFontSize": 14,
                  "titleFontSize": 13 
             },
             "legend": {
                  "labelFont": font,
                  "titleFont": font,
                  #"labelfontSize": 26,
                  "titleFontSize": 14
             },
             "text": {
                  "font": font,
                  "fontSize": 16
             }
        }
    } 

    # Register and enable the theme
    alt.themes.register('boots', boots)
    alt.themes.enable('boots')
    alt.data_transformers.disable_max_rows()

    #plotting the rarefaction curves including the shaded area
    depth_lines = pd.DataFrame({
        'x': [lower_value, knee_point, upper_value],  # x-values for vertical lines
        'label': ['-5%', 'Knee point', '+5%']  
    })

    reads_per_sample_df = reads_per_sample.reset_index()
    reads_per_sample_df.columns = ['sample', 'reads_per_sample']
    print(f"Rows reads per sample df: {len(reads_per_sample_df)}")
    print(f"Rows combined df: {len(combined_df)}")

    zoom = alt.selection_interval(bind='scales')

    param_checkbox = alt.param(
        bind=alt.binding_checkbox(name='Show read depth specified by slider as a line on the plot: ')
    )

    s = alt.param(
        name='position', bind=alt.binding_range(min=0, max=max_reads, step=20, name='Rarefaction Depth Line'), value=knee_point
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
        width=450,
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
        width=450, #600
        height=350 #400
    )
    

    final_with_line = alt.layer(final_chart, vertical_line).resolve_scale( 
        x='shared',
        y='shared'
    )

    text = [f'Knee point: {knee_point} (at the {percentile:.2f} percentile)',
            f'Knee point -5%: {lower_value:.0f} (at the {lower_percentile:.2f} percentile)',
            f'Knee point +5%: {upper_value:.0f} (at the {upper_percentile:.2f} percentile)',
            f'The shaded area is where at least {percent_samples * 100:.2f}% of the samples ',
            f'are kept, if used as the rarefaction depth.',
            f' ']
          
    if (knee_point > depth_threshold):
        ps = percent_samples * 100
        text.append(f"The knee point is above the specified, acceptable depth range!")
        text.append(f"The upper bound of the acceptable range (keeping at least {ps:.2f}% of the samples) is {depth_threshold}.")
        text.append(f" ")
        text.append(f"If the calculated knee point is used, {percentile:.2f}% of the samples will")
        text.append(f" be excluded because they have too little reads.")
        
    text_l = pd.DataFrame({'text': ['\n'.join(text)]})
    text_lines = alt.Chart(text_l).mark_text(fontSize=16, size=14, align='left', baseline='top', lineBreak="\n", dx=-95).encode(text='text:N').properties(width=100, height=50)#(width=100, height=300)
    upper_chart = alt.vconcat(final_with_line, text_lines).properties(spacing=0)


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
        width=450,
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
        width=450,
        height=350
    )
    barplot_combined = alt.layer(background, barplot).resolve_scale(
        x='shared',
        y='shared'
    )
        
    combined_chart = alt.hconcat(upper_chart, barplot_combined).properties(spacing=60).configure_legend(
        labelFontSize=14,  # Font size of the legend labels
        titleFontSize=14   # Font size of the legend title
    )
    #combined_chart.save('new_chart.png')
    
    #change css file path if I want to save it this way
    #combined_chart.save('new_chart.html')
    #change_html_file('new_chart.html')

    # Check if the file already exists and delete it if it does
    new_chart_path = os.path.join(output_dir, 'new_chart.html')
    if os.path.exists(new_chart_path):
        os.remove(new_chart_path)
    combined_chart.save(new_chart_path, inline=True)
    change_html_file(new_chart_path)
    
    #end measuring runtime & memory usage
    end_time = time.time()

    #copied from q2-moshpit/busco
    tabbed_context = {}
    vega_json = combined_chart.to_json()
    TEMPLATES = os.path.join(
        os.path.dirname(__file__),
        "assets"
    )
    
    tabbed_context.update({
        "vega_json": vega_json
    })
    templates = os.path.join(TEMPLATES, 'index.html')
    copytree(
        src=TEMPLATES,
        dst=output_dir,
        dirs_exist_ok=True
    )

    q2templates.render(templates, output_dir, context=tabbed_context)


    #end measuring runtime & memory usage
    e_time = time.time()
    current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    
    print(f"Current memory usage: {current / 10**6:.4f} MB")
    print(f"Peak memory usage: {peak / 10**6:.4f} MB")

    total_runtime = e_time - start_time
    render_time = e_time - end_time
    algorithm_time = after_algo_time - start_time
    plot_time = end_time - after_algo_time
    print(f"Total runtime: {total_runtime:.4f} seconds")
    print(f"Algorithm time: {algorithm_time:.4f} seconds")
    print(f"Plotting time: {plot_time:.4f} seconds")
    print(f"Rendering time: {render_time:.4f} seconds")

    """ #runtime csv
    # Check if the file exists
    if os.path.exists("runtime_.csv"):
        # Load existing data
        existing_data = np.loadtxt("runtime_.csv", delimiter=",", ndmin=1)
        print(existing_data)
        # Append the new value to the existing data
        updated_data = np.append(existing_data, [total_runtime, algorithm_time, plot_time, render_time])
    else:
        # If the file does not exist, create an array with the new value
        updated_data = np.array([total_runtime, algorithm_time, plot_time, render_time])
    
    # Save the updated data back to the file
    np.savetxt("runtime_.csv", updated_data, delimiter=",", fmt="%f")

    #memory csv in MB
    if os.path.exists("memory_.csv"):
        # Load existing data
        existing_data = np.loadtxt("memory_.csv", delimiter=",", ndmin=1)
        # Append the new value to the existing data
        updated_data = np.append(existing_data, (peak / 10**6))
        
    else:
        # If the file does not exist, create an array with the new value
        updated_data = np.array([peak / 10**6])
    
    # Save the updated data back to the file
    np.savetxt("memory_.csv", updated_data, delimiter=",", fmt="%f")"""

    """existing_data = np.loadtxt("memory_kneedle.csv", delimiter=",", ndmin=1)
    
    # Prepare the data for each output file
    total_time = existing_data[::4]  # 1st, 5th, 9th, etc. (index 0, 4, 8, ...)
    algo_time = existing_data[1::4]   # 2nd, 6th, 10th, etc. (index 1, 5, 9, ...)
    plot_time = existing_data[2::4]   # 3rd, 7th, 11th, etc. (index 2, 6, 10, ...)
    render_time = existing_data[3::4] # 4th, 8th, 12th, etc. (index 3, 7, 11, ...)

    data = [total_time, algo_time, plot_time, render_time]

    plt.boxplot(data)
    time_categories = ['Total', 'Algorithm', 'Plotting', 'Rendering']
    plt.xticks([1, 2, 3, 4], time_categories)
    plt.ylabel("Time in seconds")
    plt.savefig("box_plot_runtime_kneedle.png", format='png')"""

     



#to test & get outputs -> delete in the end
#feature_table_path = "../../table.qza"
#other feature tables
#feature_table_path = "../../atacama_soil_table.qza"
#feature_table_path = "../../parkinson_mouse_dada2_table.qza"
#very big one
#feature_table_path = "../../feature-table.qza"
#ft_artifact = qiime2.Artifact.load(feature_table_path)
#rarefaction_depth("../../", ft_artifact)
