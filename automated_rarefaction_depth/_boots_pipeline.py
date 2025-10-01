# ----------------------------------------------------------------------------
# Copyright (c) 2016-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------


import json
import os
import numpy as np
import pandas as pd
import qiime2
import q2templates
from qiime2 import Artifact
from kneed import KneeLocator
import altair as alt
from shutil import copytree
import warnings


def altair_theme():
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


_pipe_defaults = {
    'iterations': 10,
    'table_size': None,
    'steps': 20,
    'percent_samples': 0.8,
    'algorithm': 'kneedle', 
    'seed': 42,
    'kmer_size': 16,
    'tfidf': False,
    'max_df': 1.0,
    'min_df': 1,
    'max_features': None,
    'norm': 'None',
    'metric': 'observed_features'
}

#ignore warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=RuntimeWarning)
warnings.simplefilter(action='ignore', category=UserWarning)
 

def pipeline_boots(ctx, table, sequence=None, iterations=_pipe_defaults['iterations'], table_size=_pipe_defaults['table_size'], metric=_pipe_defaults['metric'],
                   steps=_pipe_defaults['steps'], percent_samples=_pipe_defaults['percent_samples'], algorithm=_pipe_defaults['algorithm'],
                   seed = _pipe_defaults['seed'], kmer_size=_pipe_defaults['kmer_size'], tfidf=_pipe_defaults['tfidf'], max_df=_pipe_defaults['max_df'],
                   min_df=_pipe_defaults['min_df'], max_features=_pipe_defaults['max_features'], norm=_pipe_defaults['norm']):
    
    alpha_action = ctx.get_action('boots', 'alpha')
    kmer_action = ctx.get_action('kmerizer', 'seqs_to_kmers')
    viz_action = ctx.get_action('rarefaction-depth', '_rf_visualizer_boots')

    table_df = table.view(pd.DataFrame)

    #run seqs_to_kmers if sequence is provided
    kmer_run = False
    if sequence is not None:
        print("Sequences were provided as input.")
        print("Therefore, the kmerizer is run to generate a new feature table with kmer sequences.")
        print("This new kmer feature table will be used for the analysis.")
        table, = kmer_action(table=table, sequences=sequence, kmer_size=kmer_size, tfidf=tfidf, max_df=max_df, min_df=min_df, max_features=max_features, norm=norm)
        table_df = table.view(pd.DataFrame)
        kmer_run = True
    else:
        print("No sequences were provided as input, and therefore, the kmerizer is not run.")
        print("The feature table given as input will be used for the analysis.")
   
    #adjusting table size if it's too big -> keep table_size rows
    if (table_size is not None and len(table_df) > table_size):
        table_df = table_df.sample(n=table_size, random_state=seed)
        table_df = table_df.loc[:, ~(table_df.isna() | (table_df == 0)).all(axis=0)] 

    table_size = len(table_df)
    reads_per_sample = table_df.sum(axis=1) 
    percentile = 90
    #adjust max_reads if kneedle and shannon -> maybe adjust the percentile later
    if (metric == 'shannon' and kmer_run == True):
        percentile = 30
    elif (metric == 'shannon' and kmer_run == False):
        percentile = 60
    max_reads = int(np.percentile(reads_per_sample, percentile))
    
    sorted_depths = reads_per_sample.sort_values()
    sorted_depths_pass = [int(depth) for depth in sorted_depths.tolist()] # converting float to int
    reads_per_sample_pass = [int(read) for read in reads_per_sample.tolist()]

    sample_loss_index = int(np.ceil((1-percent_samples) * len(sorted_depths))) - 1 
    depth_threshold = int(sorted_depths.iloc[sample_loss_index])

    artifacts_list = []
    sample_list = table_df.index.tolist()

    max_range = np.linspace(1, max_reads, num=steps, dtype=int)
    
    for i in range(steps):
        print(f"step {i+1}: {max_range[i]}")
        result, = alpha_action(table=table, sampling_depth=int(max_range[i]), metric=metric, n=iterations, replacement=False, average_method='mean')
        artifacts_list.append(result)

    pd_new = pd.DataFrame(
        np.nan,  
        index=table_df.index, 
        columns=[f"Step_{j+1}" for j in range(steps)]  
    )
    
    dfs = []
    for i, artifact in enumerate(artifacts_list):
        df = artifact.view(pd.Series).to_frame(name='observed_features')
        dfs.append(df)
        for j, sample in enumerate(sample_list):
            if sample in df.index:
                pd_new.iloc[j, i] = df.loc[sample].values[0].round().astype(int)
            else:
                pd_new.iloc[j, i] = np.nan

    
    combined_df, knee_point = _rf_knee_locator(artifacts_list=pd_new, sample_list=sample_list,
                                               steps=steps, algorithm=algorithm, max_reads=max_reads)
    
    sample_names = combined_df.iloc[:, -1].tolist()
    combined_df = combined_df.drop(combined_df.columns[-1], axis=1)

    combined_df.index = combined_df.index.astype(str)
    combined_artifact = qiime2.Artifact.import_data("FeatureTable[Frequency]", combined_df)

    visualization, = viz_action(sample_names=sample_names, metric=metric, kmer_run=kmer_run, percent_samples=percent_samples, reads_per_sample=reads_per_sample_pass, combined_df=combined_artifact,
                    sorted_depths=sorted_depths_pass, knee_point=knee_point, max_reads=int(max_reads), depth_threshold=int(depth_threshold), max_read_percentile=percentile)

    return visualization


    
def _rf_knee_locator(artifacts_list: pd.DataFrame, sample_list: list[str], steps: int, algorithm: str, max_reads: int) -> tuple[pd.DataFrame, int]:

    knee_points = [None] * len(sample_list)
    df_list = []
    max_range = np.linspace(1, max_reads, num=steps, dtype=int)

    for i, sample in enumerate(sample_list):
        array_sample = np.array(artifacts_list.iloc[i]).flatten() 
        sample = np.full(len(max_range), sample)  

        sample_df = pd.DataFrame({'depth': max_range, 'observed_features': array_sample, 'sample': sample})
        df_list.append(sample_df)
    
        if(algorithm == 'kneedle'):
            #using KneeLocator to find the knee point
            kneedle = KneeLocator(max_range, array_sample, curve="concave", direction="increasing", S=3)
            knee_points[i] = kneedle.knee
        else:
            #using the gradient method
            first_derivative = np.gradient(array_sample, max_range)
            second_derivative = np.gradient(first_derivative, max_range)
            max_index = np.argmax(second_derivative)
            knee_points[i] = max_range[max_index]

    combined_df = pd.concat(df_list, ignore_index=True)
    
    knee_points_filtered = [point for point in knee_points if point is not None] 
    knee_point = round(np.mean(knee_points_filtered))
    print("calculated rarefaction depth:")
    print(knee_point)

    return combined_df, knee_point

    


def _rf_visualizer_boots(output_dir: str, sample_names: list[str], percent_samples: float, reads_per_sample: list[int], sorted_depths: list[int],
                         metric: str, max_reads: int, depth_threshold: int, knee_point: int, kmer_run: bool, combined_df: pd.DataFrame,
                         max_read_percentile: int)-> None: 
    
    combined_df['sample'] = sample_names
    sorted_depths = pd.Series(sorted_depths)
    reads_per_sample = pd.DataFrame(reads_per_sample)

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
    #register and enable the defined theme
    alt.themes.register('altair_theme', altair_theme)
    alt.themes.enable('altair_theme')
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

    if kmer_run:
        title_x = 'Total Kmer Count'
        title_y = '# Observed Distinct Kmers'
        graph_data = "number of observed kmers"
        graph_name = "Kmers"
        barplot_title = 'Kmers per Sample'
        property_title = 'Histogram of Kmers per Sample'
    else:
        title_x = 'Read Depth'
        title_y = '# Observed Features'
        graph_data = "number of observed features"
        graph_name = "Reads"
        barplot_title = 'Reads per Sample'
        property_title = 'Histogram of Reads per Sample'

    if metric == 'shannon':
        title_y = 'Shannon Index'
        graph_data = "Shannon Index"

    chart = alt.Chart(combined_df).mark_line(point=True).encode( 
        x=alt.X('depth:Q', title=title_x),
        y=alt.Y('observed_features:Q', title=title_y),
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

    text_lines = alt.Chart(pd.DataFrame({'text': ['']})).mark_text(fontSize=16, size=6, align='left', baseline='top', lineBreak="\n", dx=-95).encode(text='text:N').properties(width=100, height=50)
    upper_chart = alt.vconcat(final_with_line, text_lines).properties(spacing=0)
    empty_lines = alt.Chart(pd.DataFrame({'text': ['\n\n']})).mark_text(fontSize=12, size=6, align='left', baseline='top', lineBreak="\n", dx=-95).encode(text='text:N').properties(width=100, height=50)
    upper_chart = alt.vconcat(empty_lines, upper_chart).properties(spacing=0)

    #barplot of reads_per_sample
    predicate = alt.datum.reads_per_sample >= s 

    barplot = alt.Chart(reads_per_sample_df).mark_bar().encode(
        x=alt.X('reads_per_sample:Q', bin=alt.Bin(maxbins=50), title=barplot_title),
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
        title=property_title,
        width=450,
        height=350
    )

    background = alt.Chart(reads_per_sample_df).mark_bar().encode(
        x=alt.X('reads_per_sample:Q', bin=alt.Bin(maxbins=50), title=barplot_title),  
        y=alt.Y('count()', title='# Samples'),
        color=alt.value('lightgrey')
    ).add_params(
        zoom
    ).properties(
        title=property_title,
        width=450,
        height=350
    )

    barplot_combined = alt.layer(background, barplot).resolve_scale(
        x='shared',
        y='shared'
    )
    
    empty_lines = alt.Chart(pd.DataFrame({'text': ['\n\n']})).mark_text(fontSize=12, size=6, align='left', baseline='top', lineBreak="\n", dx=-95, dy=-5).encode(text='text:N').properties(width=100, height=50)
    barplot_combined = alt.vconcat(empty_lines, barplot_combined).properties(spacing=0)

    combined_chart = alt.hconcat(upper_chart, barplot_combined).properties(spacing=60).configure_legend(
        labelFontSize=14,  
        titleFontSize=14   
    )

    combined_chart["padding"] = {
        "top": 0,
        "left": 0,
        "right": 0,
        "bottom": -49  
    }

    vega_json = combined_chart.to_json()

    TEMPLATES = os.path.join(
        os.path.dirname(__file__),
        "assets"
    )

    add_text = False
    if (knee_point > depth_threshold):
        add_text = True

    percent_samples_100 = round(percent_samples * 100, 2)
    percentile = round(percentile, 2)
    lower_percentile = round(lower_percentile, 2)
    upper_percentile = round(upper_percentile, 2)
    lower_value = round(lower_value, 0)
    upper_value = round(upper_value, 0)

    tabbed_context = {
        "vega_json": vega_json,
        "knee_point": json.dumps(int(knee_point)),
        "percent_samples_100": json.dumps(float(percent_samples_100)),
        "depth_threshold": json.dumps(int(depth_threshold)),
        "add_text": bool(add_text),
        "percentile": json.dumps(float(percentile)),
        "lower_percentile": json.dumps(float(lower_percentile)),
        "upper_percentile": json.dumps(float(upper_percentile)),
        "lower_value": json.dumps(int(lower_value)),
        "upper_value": json.dumps(int(upper_value)),
        "graph_data": graph_data,
        "graph_name": graph_name,
        "max_read_percentile": json.dumps(int(max_read_percentile))
    }
    
    templates = os.path.join(TEMPLATES, 'index.html')

    copytree(
        src=TEMPLATES,
        dst=output_dir,
        dirs_exist_ok=True 
    )
   
    q2templates.render(templates, output_dir, context=tabbed_context)
    