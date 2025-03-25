# ----------------------------------------------------------------------------
# Copyright (c) 2016-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------


import os
from urllib.parse import quote
import numpy as np
import pandas as pd
import qiime2
import q2templates
from qiime2 import Artifact
from kneed import KneeLocator
import altair as alt
import time
import tracemalloc
from shutil import copytree
import shutil
from bs4 import BeautifulSoup
from qiime2 import sdk
from tempfile import TemporaryDirectory
from qiime2.plugin import Visualization
#import warnings


# change so that q2-boots is used and it's a pipeline, not a visualizer anymore
# semester project 


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
    #add css stylesheet to the html file
    with open(file_path, 'r') as file:
        soup = BeautifulSoup(file, 'html.parser')
        link_tag = soup.new_tag('link', rel='stylesheet', href='./css/styles.css')
        soup.head.append(link_tag)

    with open(file_path, 'w') as file:
        file.write(str(soup))



_pipe_defaults = {
    'seed': 42,
    'iterations': 2, #10
    'table_size': None,
    'steps': 5, #20
    'percent_samples': 0.8,
    'algorithm': 'gradient'
}

#doesn't work for an inexplicable reason, delete later
def rf_depth_pipe(ctx, table, seed=_pipe_defaults['seed'], iterations=_pipe_defaults['iterations'], table_size=_pipe_defaults['table_size'],
                   steps=_pipe_defaults['steps'], percent_samples=_pipe_defaults['percent_samples'], algorithm=_pipe_defaults['algorithm']):#output_dir: str,
    print("started!")
    alpha_action = ctx.get_action('boots', 'alpha')

    #warnings.simplefilter(action='ignore', category=FutureWarning)

    min_depth = 1
    table_df = table.view(pd.DataFrame)

    if table_df.empty:
        raise ValueError("The feature table is empty.")
    if not np.issubdtype(table_df.values.dtype, np.number):
        raise ValueError("The feature table contains non-numerical values.")
    

    #adjusting table size if it's too big -> keep table_size rows
    if (table_size is not None and len(table_df) > table_size):
        table_df = table_df.sample(n=table_size, random_state=seed) 
        table_df = table_df.loc[:, ~(table_df.isna() | (table_df == 0)).all(axis=0)] 

    print("after table size adjustment")
    num_samples = len(table_df)
    reads_per_sample = table_df.sum(axis=1) 
    max_reads = reads_per_sample.max()
    sorted_depths = reads_per_sample.sort_values()

    sample_loss_index = int(np.ceil((1-percent_samples) * len(sorted_depths))) - 1 
    depth_threshold = int(sorted_depths.iloc[sample_loss_index])

    knee_points = [None] * num_samples
    s = 0
    df_list = []

    print("before loop to calculate rf")
    #go through every sample and calculate the observed features for each depth
    for sample in table_df.index:
        max_range = np.linspace(min_depth, reads_per_sample.loc[sample], num=steps, dtype=int)
        array_sample = np.empty(steps)

        sample_values = table_df.loc[sample].values	

        for i in range(steps):
            print("before alpha_action")
            table_artifact = Artifact.import_data('FeatureTable[Frequency]', table_df)
            result, = alpha_action(table=table_artifact, sampling_depth=int(max_range[i]), metric='observed_features', n=iterations, replacement=False, average_method='mean')
            metadata = result.view(qiime2.Metadata)
            rarefied_sample = metadata.to_dataframe()
            #rarefied_sample = rarefy(sample_values, temp, j, seed)
            array_sample[i] = rarefied_sample.loc[sample, "observed_features"]#rarefied_sample[sample]
            print(array_sample[i])  
        
        sample_df = pd.DataFrame({'depth': max_range, 'observed_features': array_sample, 'sample': sample})
        df_list.append(sample_df)
        if(algorithm.lower().strip() == 'kneedle'):
            #using KneeLocator to find the knee point 
            kneedle = KneeLocator(max_range, array_sample, curve="concave", direction="increasing")
            knee_points[s] = kneedle.knee
        else:
            #using the gradient method
            curr_array = array_sample
            first_derivative = np.gradient(curr_array, max_range)
            second_derivative = np.gradient(first_derivative, max_range)
            max_index = np.argmax(second_derivative)
            knee_points[s] = max_range[max_index]

        s += 1

    combined_df = pd.concat(df_list, ignore_index=True)
    print("after rf loop")
    
    rf_depth(output_dir=ctx.output_dir, percent_samples=percent_samples,
              reads_per_sample=reads_per_sample, knee_points=knee_points, sorted_depths=sorted_depths, max_reads=max_reads, combined_df=combined_df, depth_threshold=depth_threshold)
    print("after calling rf_depth")
    
    #visualization = rf_depth(table=result, output_dir=output_dir, seed=seed, iterations=iterations, table_size=table_size, steps=steps, percent_samples=percent_samples, algorithm=algorithm)
    visualization = qiime2.Visualization.load(os.path.join(output_dir, "output.qzv"))
    #rf_depth(table=result, output_dir='/home/andrea/automated-rarefaction-depth/', seed=seed, iterations=iterations, table_size=table_size, steps=steps, percent_samples=percent_samples, algorithm=algorithm)
    print("calculated rf_depth")
    return visualization
    #return qiime2.Visualization.load('/home/andrea/automated-rarefaction-depth/rarefaction-depth.qzv')
    


def pipeline_test_new(ctx, table, seed=_pipe_defaults['seed'], iterations=_pipe_defaults['iterations'], table_size=_pipe_defaults['table_size'],
                   steps=_pipe_defaults['steps'], percent_samples=_pipe_defaults['percent_samples'], algorithm=_pipe_defaults['algorithm']):#output_dir: str,
    print("started!")
    alpha_action = ctx.get_action('boots', 'alpha')


    min_depth = 1
    table_df = table.view(pd.DataFrame)

    if table_df.empty:
        raise ValueError("The feature table is empty.")
    if not np.issubdtype(table_df.values.dtype, np.number):
        raise ValueError("The feature table contains non-numerical values.")
    

    #adjusting table size if it's too big -> keep table_size rows
    #delete this later:
    table_size = 10
    if (table_size is not None and len(table_df) > table_size):
        table_df = table_df.sample(n=table_size, random_state=seed) 
        table_df = table_df.loc[:, ~(table_df.isna() | (table_df == 0)).all(axis=0)] 

    print("after table size adjustment")
    #num_samples = len(table_df)
    reads_per_sample = table_df.sum(axis=1) 
    max_reads = reads_per_sample.max()
    sorted_depths = reads_per_sample.sort_values() 

    #pass as list of integers 
    sorted_depths_pass = sorted_depths.tolist()
    sorted_depths_pass = [int(depth) for depth in sorted_depths_pass]
    reads_per_sample_pass = reads_per_sample.tolist()
    reads_per_sample_pass = [int(read) for read in reads_per_sample_pass]

    sample_loss_index = int(np.ceil((1-percent_samples) * len(sorted_depths))) - 1 
    depth_threshold = int(sorted_depths.iloc[sample_loss_index])

    #knee_points = [None] * num_samples
    #s = 0
    #df_list = []
    artifacts_list = []
    sample_list = table_df.index.tolist()
    depths_list = []

    print("before loop to calculate rf")
    #go through every sample and calculate the observed features for each depth
    for sample in table_df.index:
        max_range = np.linspace(min_depth, reads_per_sample.loc[sample], num=steps, dtype=int)
        depths_list.append(max_range.tolist())

        #array_sample = np.empty(steps)

        #sample_values = table_df.loc[sample].values	
        # Create a new DataFrame with only the current sample
        current_sample_df = table_df.loc[[sample]]

        for i in range(steps):
            print("before alpha_action")
            #subsample so there's only the current sample in the table
            table_artifact = Artifact.import_data('FeatureTable[Frequency]', current_sample_df)
            result, = alpha_action(table=table_artifact, sampling_depth=int(max_range[i]), metric='observed_features', n=iterations, replacement=False, average_method='mean')
            artifacts_list.append(result)

            """metadata = result.view(qiime2.Metadata)
            rarefied_sample = metadata.to_dataframe()"""
            #rarefied_sample = rarefy(sample_values, temp, j, seed)
            
            #array_sample[i] = rarefied_sample.loc[sample, "observed_features"]#rarefied_sample[sample] 
            #pass list of artifacts to visualizer
        
        """sample_df = pd.DataFrame({'depth': max_range, 'observed_features': array_sample, 'sample': sample})
        df_list.append(sample_df)
        if(algorithm.lower().strip() == 'kneedle'):
            #using KneeLocator to find the knee point 
            kneedle = KneeLocator(max_range, array_sample, curve="concave", direction="increasing")
            knee_points[s] = kneedle.knee
        else:
            #using the gradient method
            curr_array = array_sample
            first_derivative = np.gradient(curr_array, max_range)
            second_derivative = np.gradient(first_derivative, max_range)
            max_index = np.argmax(second_derivative)
            knee_points[s] = max_range[max_index]

        s += 1"""

    #combined_df = pd.concat(df_list, ignore_index=True)
    print("after rf loop")
    
    _rf_visualizer(output_dir='/home/andrea/automated-rarefaction-depth/', percent_samples=percent_samples, reads_per_sample=reads_per_sample_pass, artifacts_list=artifacts_list, 
                   sorted_depths=sorted_depths_pass, max_reads=max_reads, depth_threshold=depth_threshold, sample_list=sample_list, depths_list=depths_list, steps=steps, algorithm=algorithm)
    print("after calling rf_visualizer")
    

    #rf_depth(table=result, output_dir='/home/andrea/automated-rarefaction-depth/', seed=seed, iterations=iterations, table_size=table_size, steps=steps, percent_samples=percent_samples, algorithm=algorithm)
    #print("calculated rf_depth")
    #change to real return value!! just a placeholder for right now to not get errors!
    return qiime2.Visualization.load('/home/andrea/automated-rarefaction-depth/rarefaction-depth.qzv')
    
    



def _rf_visualizer(output_dir: str, percent_samples: float, reads_per_sample: set[int], artifacts_list: set[pd.DataFrame], sorted_depths: set[int], max_reads: int, depth_threshold: int, sample_list: set[str], depths_list: set[int], steps: int, algorithm: str)-> None:
    
    print("in rf_visualizer")
    sorted_depths = pd.Series(sorted_depths)
    if isinstance(reads_per_sample, list):
        print("reads_per_sample is a list")
        reads_per_sample = pd.DataFrame(reads_per_sample)

    counter = 0
    knee_points = [None] * len(sample_list)
    df_list = []
    pd_list = [artifact.view(qiime2.Metadata).to_dataframe() for artifact in artifacts_list]

    for sample in sample_list:
        max_range = np.array(depths_list[counter])
        array_sample = np.array(pd_list[counter * steps : (counter + 1) * steps])

        print(f"array_sample shape: {(array_sample)}")
        array_sample = array_sample.flatten()  
        sample = np.full(len(max_range), sample)  # Create an array of repeated values


        print(f"max_range shape: {np.shape(max_range)}")
        print(f"array_sample shape: {np.shape(array_sample)}")
        print(f"sample shape: {np.shape(sample)}")

        sample_df = pd.DataFrame({'depth': max_range, 'observed_features': array_sample, 'sample': sample})
        df_list.append(sample_df)
        if(algorithm.lower().strip() == 'kneedle'):
            #using KneeLocator to find the knee point 
            kneedle = KneeLocator(max_range, array_sample, curve="concave", direction="increasing")
            knee_points[counter] = kneedle.knee
        else:
            #using the gradient method
            curr_array = array_sample
            print("curr_array:")
            print(curr_array)
            print(type(curr_array))
            print(curr_array.shape)
            print("max_range:")
            print(max_range)
            print(type(max_range))
            print(max_range.shape)
            first_derivative = np.gradient(curr_array, max_range)
            second_derivative = np.gradient(first_derivative, max_range)
            max_index = np.argmax(second_derivative)
            knee_points[counter] = max_range[max_index]

        counter += 1

    print("after knee points calculations")
    combined_df = pd.concat(df_list, ignore_index=True)
    
    knee_points_filtered = [point for point in knee_points if point is not None]
    knee_point = round(np.mean(knee_points_filtered))
    print("knee_point:")
    print(knee_point)
    print("knee_points:")
    print(knee_points)
   

    #finding +-5% points & at what percentile knee point is
    index = np.searchsorted(sorted_depths, knee_point)
    percentile = (index / len(sorted_depths)) * 100

    lower_percentile = max(percentile - 5, 0) 
    upper_percentile = min(percentile + 5, 100) 

    lower_index = int((lower_percentile / 100) * len(sorted_depths))
    upper_index = min(int((upper_percentile / 100) * len(sorted_depths)), len(sorted_depths) - 1)

    lower_value = sorted_depths.iloc[lower_index]
    upper_value = sorted_depths.iloc[upper_index]

    #plotting with altair
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
                  "titleFontSize": 14
             },
             "text": {
                  "font": font,
                  "fontSize": 16
             }
        }
    } 

    #register and enable the defined theme
    alt.themes.register('boots', boots)
    alt.themes.enable('boots')
    alt.data_transformers.disable_max_rows()

    #plotting the rarefaction curves including the shaded area
    depth_lines = pd.DataFrame({
        'x': [lower_value, knee_point, upper_value],  
        'label': ['-5%', 'Knee point', '+5%']  
    })

    reads_per_sample_df = reads_per_sample.reset_index()
    reads_per_sample_df.columns = ['sample', 'reads_per_sample']

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
        color=alt.Color('sample:N', legend=None).scale(scheme='category10')  
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
    my_colors = ["#377eb8", "#f781bf", "#984ea3"]

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
        width=450, 
        height=350 
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
    text_lines = alt.Chart(text_l).mark_text(fontSize=16, size=14, align='left', baseline='top', lineBreak="\n", dx=-95).encode(text='text:N').properties(width=100, height=50)
    upper_chart = alt.vconcat(final_with_line, text_lines).properties(spacing=0)
    empty_lines = alt.Chart(pd.DataFrame({'text': ['\n\n']})).mark_text(fontSize=12, size=10, align='left', baseline='top', lineBreak="\n", dx=-95).encode(text='text:N').properties(width=100, height=50)
    upper_chart = alt.vconcat(empty_lines, upper_chart).properties(spacing=0)


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
    
    empty_lines = alt.Chart(pd.DataFrame({'text': ['\n\n']})).mark_text(fontSize=12, size=10, align='left', baseline='top', lineBreak="\n", dx=-95, dy=-5).encode(text='text:N').properties(width=100, height=50)
    barplot_combined = alt.vconcat(empty_lines, barplot_combined).properties(spacing=0)
        
    combined_chart = alt.hconcat(upper_chart, barplot_combined).properties(spacing=60).configure_legend(
        labelFontSize=14,  
        titleFontSize=14   
    )
   
    # Check if the file already exists and delete it if it does
    new_chart_path = os.path.join(output_dir, 'new_chart.html')
    if os.path.exists(new_chart_path):
        os.remove(new_chart_path)
    combined_chart.save(new_chart_path, inline=True)
    change_html_file(new_chart_path)
    

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

    output_dir_n = os.path.join(output_dir, 'q2templateassets')
    if os.path.exists(output_dir_n):
        shutil.rmtree(output_dir_n)

    copytree(
        src=TEMPLATES,
        dst=output_dir_n,#output_assets,
        dirs_exist_ok=True
    )
    print("output_dir:")
    print(output_dir_n)
    
    q2templates.render(templates, output_dir_n, context=tabbed_context)
    """output_html_path = os.path.join(output_dir, 'pipe_output.qzv')
    provenance = qiime2.sdk.parse_context().provenance
    visualization = qiime2.Visualization._from_data_dir(output_dir, provenance_capture=provenance)
    visualization.save(output_html_path)"""
    #return qiime2.Visualization.load(output_html_path)
    

