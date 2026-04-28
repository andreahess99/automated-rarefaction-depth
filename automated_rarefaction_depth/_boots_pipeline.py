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
#from qiime2 import Artifact
from kneed import KneeLocator
import altair as alt
from shutil import copytree
import warnings
from skbio import DistanceMatrix
#import time
#import matplotlib.pyplot as plt


"""def altair_theme():
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
    }"""


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
    'metrics': {'observed_features'}
}


#ignore warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=RuntimeWarning)
warnings.simplefilter(action='ignore', category=UserWarning)
 

def pipeline_boots(ctx, table, meta_data, sequence=None, iterations=_pipe_defaults['iterations'], table_size=_pipe_defaults['table_size'], metrics=_pipe_defaults['metrics'],
                   steps=_pipe_defaults['steps'], percent_samples=_pipe_defaults['percent_samples'], algorithm=_pipe_defaults['algorithm'],
                   seed = _pipe_defaults['seed'], kmer_size=_pipe_defaults['kmer_size'], tfidf=_pipe_defaults['tfidf'], max_df=_pipe_defaults['max_df'],
                   min_df=_pipe_defaults['min_df'], max_features=_pipe_defaults['max_features'], norm=_pipe_defaults['norm']):
    
    #alpha_action = ctx.get_action('boots', 'alpha')
    beta_action = ctx.get_action('boots', 'beta')
    kmer_action = ctx.get_action('kmerizer', 'seqs_to_kmers')
    viz_combined_action = ctx.get_action('rarefaction-depth', '_combined_viz')
    alpha_collection_action = ctx.get_action("boots", "alpha_collection")

    table_df = table.view(pd.DataFrame)
    alpha = False
    beta = False
    #observed_features and braycurtis are always included
    metrics.add('observed_features')
    metrics.add('braycurtis')

    if any(m in ['braycurtis', 'jaccard', 'hamming', 'dice', 'jensenshannon', 'matching', 'minkowski', 'rogerstanimoto', 'russellrao',
                        'seuclidean', 'sokalmichener', 'sokalsneath', 'sqeuclidean', 'yule', 'canberra_adkins', 'chebyshev', 'cityblock', 'correlation', 'cosine',
                        'euclidean', 'aitchison',  'canberra'] for m in metrics):
        beta = True
        print("Beta", beta)
        metrics_beta = [m for m in metrics if m in ['braycurtis', 'jaccard', 'hamming', 'dice', 'jensenshannon', 'matching', 'minkowski', 'rogerstanimoto', 'russellrao',
                        'seuclidean', 'sokalmichener', 'sokalsneath', 'sqeuclidean', 'yule', 'canberra_adkins', 'chebyshev', 'cityblock', 'correlation', 'cosine',
                        'euclidean', 'aitchison',  'canberra']]
        print("metrics beta:", metrics_beta)
        
    if any(m not in metrics_beta for m in metrics):
        alpha = True 
        metrics_alpha = [m for m in metrics if m not in metrics_beta] #adjust this list later
        print("metrics alpha:", metrics_alpha)

    meta = meta_data.to_dataframe()
    meta.index.name = "sample"
    metadata_columns = ["sample"] + meta.columns.tolist()
    
    print(metrics)
    metric = list(metrics)[0] #for now just take the first metric in the set, so I don't get errors from my beta code
    
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
    """if (metric == 'shannon' and kmer_run == True):
        percentile = 30
    elif (metric == 'shannon' and kmer_run == False):
        percentile = 60"""
    max_reads = int(np.percentile(reads_per_sample, percentile))
    
    sorted_depths = reads_per_sample.sort_values()
    sorted_depths_pass = [int(depth) for depth in sorted_depths.tolist()] 
    reads_per_sample_pass = [int(read) for read in reads_per_sample.tolist()]

    meta.fillna(0, inplace=True)
    reads_per_sample_merged = pd.DataFrame({"sample": table_df.index.tolist(), "reads": reads_per_sample_pass})
    reads_per_sample_merged = reads_per_sample_merged.merge(meta, left_on="sample", right_index=True, how="left")     
    reads_per_sample_merged.index.name = "sample"
    reads_per_sample_merged = (reads_per_sample_merged.rename(columns={"sample": "sample-id"}).set_index("sample-id"))
   

    sample_loss_index = int(np.ceil((1-percent_samples) * len(sorted_depths))) - 1 
    depth_threshold = int(sorted_depths.iloc[sample_loss_index])

    sample_list = table_df.index.tolist()

    max_range = np.linspace(1, max_reads, num=steps, dtype=int)

    # beta metric specific code
    if beta:
        knee_points_beta = []
        data_beta = []
        for k, metric in enumerate(metrics_beta):
            print("metric:", metric)
            num_samples_left = [None] * (steps)
            avg_difference = [None] * (steps-1)
            median_difference = [None] * (steps-1)
            std_difference = [None] * (steps-1)
            p75_25_difference = [None] * (steps-1)
            avg_range = [None] * (steps-1)
            #workaround as call fails for sampling_depth=1
            max_range[0] = 2
            for i in range(steps):
                print(f"step {i+1}: {max_range[i]}")
                #beta_result, = beta_action(table=table, sampling_depth=1, metric=metric, n=iterations, replacement=False)
                beta_result, = beta_action(table=table, sampling_depth=(int(max_range[i])), metric=metric, n=iterations, replacement=False)
                if k == 0:
                    num_samples_left[i] = beta_result.view(DistanceMatrix).shape[0] 

                if i > 0:
                    # reduce old_beta_result to only the ids present in the current beta_result
                    old_dm_full = old_beta_result.view(DistanceMatrix)
                    new_ids = list(beta_result.view(DistanceMatrix).ids)
                    common_ids = [sid for sid in old_dm_full.ids if sid in new_ids]

                    if common_ids:
                        idxs = [old_dm_full.ids.index(sid) for sid in common_ids]
                        subdata = old_dm_full.data[np.ix_(idxs, idxs)]
                        old_dm = DistanceMatrix(subdata, ids=common_ids)
                
                    new_dm = beta_result.view(DistanceMatrix)
                    difference = np.abs(new_dm.data - old_dm.data)
                    
                    avg_difference[i-1] = np.mean(difference[np.triu_indices_from(difference, k=1)])
                    median_difference[i-1] = np.median(difference[np.triu_indices_from(difference, k=1)])
                    std_difference[i-1] = np.std(difference[np.triu_indices_from(difference, k=1)])
                    p75_25_difference[i-1] = np.percentile(difference[np.triu_indices_from(difference, k=1)], 75) - np.percentile(difference[np.triu_indices_from(difference, k=1)], 25)

                    avg_range[i-1] = ((max_range[i] - max_range[i-1]) / 2) + max_range[i-1]
                
                old_beta_result = beta_result

            clean_max_range = [float(x) for x in max_range]
            clean_avg_diff = [float(x) if x is not None else 0.0 for x in avg_difference]
            data_beta.append(pd.DataFrame({'metric': metric, 'depth': avg_range[1:], 'observed': clean_avg_diff[1:]}))
            if k==0:
                #clean_samples_left = [int(x) for x in num_samples_left]
                df_bars = pd.DataFrame({'depth': avg_range[1:], 'num_samples_left': num_samples_left[2:]})
            if metric in ['braycurtis', 'jaccard', 'sokalsneath', 'matching', 'cosine', 'yule', 'canberra_adkins', 'jensenshannon', 'hamming',
                          'dice', 'correlation', 'aitchison', 'canberra', 'rogerstanimoto', 'sokalmichener', 'russellrao']:
                curve_type = "convex"
                direction = "decreasing"
            kpb = knee_point_locator(avg_range[1:], clean_avg_diff[1:], algorithm, "convex", "decreasing") #curve_type, direction
            print("knee point for metric", metric, ":", kpb)
            kpb = round(float(kpb)) if kpb is not None else 0
            knee_points_beta.append(pd.DataFrame({'knee': kpb, 'metric': metric}, index=[0]))
            print("knee_points_beta:", knee_points_beta)


        data_beta = pd.concat(data_beta, ignore_index=True)
        data_beta.columns = ['metric', 'depth', 'observed']
        data_beta.insert(0, 'id', [f"row{i}" for i in range(len(data_beta))])
        
        data_beta = data_beta.set_index('id')
        data_beta = qiime2.Metadata(data_beta)

        kp_beta = pd.concat(knee_points_beta, ignore_index=True)
        kp_beta.columns = ['knee', 'metric']
        kp_beta.insert(0, 'id', [f"row{i}" for i in range(len(kp_beta))])
        kp_beta = kp_beta.set_index('id')
        kp_beta.index = kp_beta.index.astype(str)
        kp_list_beta = qiime2.Metadata(kp_beta)
        df_bars.insert(0, 'id', [f"row{i}" for i in range(len(df_bars))])
        df_bars = df_bars.set_index('id')
        num_samples = qiime2.Metadata(df_bars)
        #visualization, = viz_combined_action(max_range=clean_max_range, kmer_run=kmer_run, metric=metric, algorithm=algorithm, num_samples=num_samples, data_beta=data_beta, kp_list_beta=kp_list_beta, beta_metrics=metrics_beta)
        
    
    if alpha:
        #if alpha metric was chosen
        dfs = []
        combined_dfs = []
        knee_point_list = []

        for metric in metrics_alpha:
            print("metric:", metric)
            for i in range(steps):
                print(f"step {i+1}: {max_range[i]}")
                """result, = alpha_action(table=table, sampling_depth=int(max_range[i]), metric=metric, n=iterations, replacement=False, average_method='mean')
                artifacts_list.append(result)"""

                #function to use if we want to get the data for each iteration and sample instead of just the average
                #and then combine it into a dataframe to use for the visualization and knee point calculation
                sample_data, = alpha_collection_action(table=table, sampling_depth=int(max_range[i]), metric=metric, n=iterations, replacement=False)
                for key, artifact in sample_data.items():
                    series = artifact.view(pd.Series)
                    df = series.reset_index()
                    df.columns = ['sample', 'observed']
                    df['iteration'] = key
                    df['read_depth'] = int(max_range[i])  
                    df['metric'] = metric
                    dfs.append(df)
                    
                combined = pd.concat(dfs, ignore_index=True)
                #potentially do it after the loop so you only merge once 
                meta.fillna(0, inplace=True)
                combined = combined.merge(meta, left_on="sample", right_index=True, how="left")

                # Calculate average after collecting all iterations for this depth
                #making the mean_df which I would need for the knee point calculation
                mean_df = (combined.groupby(['sample', 'read_depth'], as_index=False).agg(mean_observed=('observed', 'mean')))
            
            combined_dfs.append(combined)
            #calculatin knee point here as data formats are a bit different
            if metric in ['observed_features', 'shannon', 'simpson', 'brillouin_d', 'chao1', 'enspie', 'goods_coverage', 'michaelis_menten_fit']:
                curve_type = "concave"
                direction = "increasing"
            elif metric in ['dominance', 'robbins', 'simpson_e', 'mcintosh_e', 'berger_parker_d', 'jaccard', 'braycurtis']:
                curve_type = "convex"
                direction = "decreasing"
            knee_points = [None] * len(sample_list)
            for i, sample in enumerate(sample_list):
                array_sample = mean_df[mean_df['sample'] == sample]['mean_observed'].values
                max_range_for_sample = max_range[:len(array_sample)]
                knee_points[i] = knee_point_locator(max_range_for_sample, array_sample, algorithm, curve_type, direction)
                
        
            knee_points_filtered = [point for point in knee_points if point is not None] 
            knee_point = round(np.mean(knee_points_filtered))
            print("calculated rarefaction depth:")
            print(knee_point)
            knee_point_list.append((knee_point, metric))

    
        combined_df = mean_df.pivot(index='sample', columns='read_depth', values='mean_observed').reset_index()
        combined_df = combined_df.drop(combined_df.columns[-1], axis=1)
        combined_df.index = combined_df.index.astype(str)
    
        percent_samples_100 = round(percent_samples * 100, 2)
 
        combined.insert(0, 'id', [f"row{i}" for i in range(len(combined))])
        combined = combined.set_index('id')
       
        combined = qiime2.Metadata(combined)
        
        kp_df = pd.DataFrame(knee_point_list, columns=['knee', 'metric'])
        kp_df.index.name = 'id'
        kp_df.index = kp_df.index.astype(str)
        knee_point_list = qiime2.Metadata(kp_df)
        metrics = list(metrics)
        
        visualization, = viz_combined_action(metric=metric, kmer_run=kmer_run, percent_samples_100=percent_samples_100, steps=steps, algorithm=algorithm,
                        sorted_depths=sorted_depths_pass, knee_point=knee_point, max_reads=int(max_reads), depth_threshold=int(depth_threshold), max_read_percentile=percentile, 
                        combined=combined, metadata_columns=metadata_columns, rps=qiime2.Metadata(reads_per_sample_merged), kp_list=knee_point_list,
                        kp_list_beta=kp_list_beta, data_beta=data_beta, alpha_metrics=metrics_alpha, beta_metrics=metrics_beta, num_samples=num_samples, max_range=clean_max_range) 

    return visualization

#calculates the knee point based on the chosen algorithm & metric
def knee_point_locator(range: list[float], samples: list[float], algorithm: str, curve_type:str, direction:str) -> float:
    if algorithm == 'kneedle':
        kneedle = KneeLocator(range, samples, curve=curve_type, direction=direction, S=3)
        knee_point = kneedle.knee
    else:
        first_derivative = np.gradient(samples, range)
        second_derivative = np.gradient(first_derivative, range)
        knee_point = np.argmax(second_derivative)   
    return knee_point


#calculates the knee point and makes the visualization for beta rarefaction
"""def _beta_viz(output_dir: str, max_range: list[float], kmer_run: bool, calc_array: list[float] , metric: str, algorithm: str, num_samples_left: list[int])->None:
 
    avg_range = [None] * (len(max_range)-1)
    for i in range(1, len(max_range)):
        avg_range[i-1] = ((max_range[i] - max_range[i-1]) / 2) + max_range[i-1]
    
    #calculate knee  point
    knee_point = knee_point_locator(avg_range[1:], calc_array[1:], algorithm, "beta", "increasing")
    if knee_point is None:
        knee_point = 0
    else:
        knee_point = round(float(knee_point))

    #plotting with altair
    alt.themes.register('altair_theme', altair_theme)
    alt.themes.enable('altair_theme')

    #make 2 plots again 1 scatter, 1 barplot with same functionality as alpha_viz with zoom etc
    zoom = alt.selection_interval(bind='scales')
    param_checkbox = alt.param(
        bind=alt.binding_checkbox(name='Show read depth specified by slider as a line on the plot: ')
    ) 
    s = alt.param(
        name='position', bind=alt.binding_range(min=0, max=avg_range[-1], step=20, name='Rarefaction Depth Line'), value=knee_point
    )

    df = pd.DataFrame({'max_range': avg_range[1:], 'calc_array': calc_array[1:]})
    base = alt.Chart(df).mark_line(point=True).encode(
            x=alt.X('max_range:Q', title='Read Depth'),
            y=alt.Y('calc_array:Q', title='Calculated Value'), #adjust to what it is 
        ).properties(
            width=400,
            height=300,
            title='Beta Rarefaction Curve'
        ).add_params(zoom, param_checkbox) 
    
    static_line = alt.Chart(pd.DataFrame({'position': [knee_point]})).mark_rule(
        color='red', strokeWidth=2).encode(x='position:Q').properties()
    
    moving_line = alt.Chart(pd.DataFrame({'position': [0]})).mark_rule(
        color='black', strokeWidth=2, strokeDash=[4, 4]).encode(
            x='position:Q'
        ).add_params(s).transform_calculate(position='position').transform_filter(param_checkbox)
    

    beta_rf_plot = alt.layer(base, static_line, moving_line).resolve_scale(x='shared', y='shared')

    #barplot
    df_bars = pd.DataFrame({'max_range': max_range, 'num_samples_left': num_samples_left})
    bar_chart = alt.Chart(df_bars).mark_bar(color='steelblue').encode(
            x=alt.X('max_range:Q', title='Read Depth'),
            y=alt.Y('num_samples_left:Q', title='Number of Samples Left')
        ).properties(
            width=400,
            height=300,
            title='Samples Remaining per Rarefaction Depth')

    # trying to fix the layout
    text_lines = alt.Chart(pd.DataFrame({'text': ['']})).mark_text(fontSize=16, size=6, align='left', baseline='top', lineBreak="\n", dx=-95).encode(text='text:N').properties(width=100, height=50)
    upper_chart = alt.vconcat(beta_rf_plot, text_lines).properties(spacing=0)
    empty_lines = alt.Chart(pd.DataFrame({'text': ['\n\n']})).mark_text(fontSize=12, size=6, align='left', baseline='top', lineBreak="\n", dx=-95).encode(text='text:N').properties(width=100, height=50)
    upper_chart = alt.vconcat(empty_lines, upper_chart).properties(spacing=0)

    bar_space = alt.vconcat(bar_chart, text_lines).properties(spacing=0)
    bar_space = alt.vconcat(empty_lines, bar_space).properties(spacing=0)

    combined_chart = alt.hconcat(upper_chart, bar_space).properties(spacing=60).configure_legend(
        labelFontSize=14,  
        titleFontSize=14   
    )
    
    combined_chart["padding"] = {
        "top": 0,
        "left": 0,
        "right": 0,
        "bottom":  -49#0  
    }
    #define and make all necessary files & definitions
    vega_json = combined_chart.to_json() 

    TEMPLATES = os.path.join(
        os.path.dirname(__file__),
        "assets"
    )
    
    #add some text somewhere if kmer was run?
    beta_context = { 
        "vega_json": vega_json,
        "beta_metric": (str(metric)),
        "knee_point": (knee_point),
        "beta": (True),
        "algorithm": (str(algorithm)),
        "percent_samples_100": (0),
        "depth_threshold": (0),
        "add_text": (False),
        "percentile": (0),
        "lower_percentile": (0),
        "upper_percentile": (0),
        "lower_value": (0),
        "upper_value": (0),
        "graph_data": (str(metric)),
        "graph_name": (str(metric)),
        "max_read_percentile": (0)
    }
    
    templates = os.path.join(TEMPLATES, 'index.html')

    copytree(
        src=TEMPLATES,
        dst=output_dir,
        dirs_exist_ok=True 
    )
   
    q2templates.render(templates, output_dir, context=beta_context)"""


# combined visualization function for alpha and beta metrics
# to do: add some text somewhere if kmer was run?
def _combined_viz(output_dir: str, metric: str, kmer_run: bool, max_range: list[float] = None, algorithm: str = "kneedle", num_samples: qiime2.Metadata = None, steps: int = None, 
                  percent_samples_100: float = 0, sorted_depths: list[int] = None, max_reads: int = 1, depth_threshold: int = 1, knee_point: int = 0, max_read_percentile: int = 1,
                  metadata_columns: list[str] = None, combined: qiime2.Metadata = None, rps:qiime2.Metadata = None, kp_list: qiime2.Metadata = None,
                  kp_list_beta: qiime2.Metadata = None, data_beta: qiime2.Metadata = None, alpha_metrics: list[str] = None, beta_metrics: list[str] = None)->None:  
    
    #default values for the tabbed_context
    add_text = False
    percentile = 0
    lower_percentile = 0
    upper_percentile = 0
    lower_value = 0
    upper_value = 0
    graph_data = str(metric)
    graph_name = str(metric)
    beta = False
    alpha = False


    # beta metric specific code
    if beta_metrics is not None and len(beta_metrics) > 0:
        beta = True
        #line_plot_title = 'Beta Rarefaction Curve'

        line_chart_df = data_beta.to_dataframe().reset_index()
        line_chart_df = line_chart_df.drop('id', axis=1)
        kp_list_beta = kp_list_beta.to_dataframe().reset_index()
        kp_list_beta = kp_list_beta.drop('id', axis=1)
        kp_list_beta = kp_list_beta.to_dict(orient='records')
        df_bars = num_samples.to_dataframe().reset_index()
        df_bars = df_bars.drop('id', axis=1)


    #alpha metric specific code
    if alpha_metrics is not None and len(alpha_metrics) > 0:
        alpha = True
        kp_list = kp_list.to_dataframe().reset_index()
        kp_list = kp_list.drop('id', axis=1)
        kp_list = kp_list.to_dict(orient='records')
        sorted_depths = pd.Series(sorted_depths)
        #reads_per_sample = pd.DataFrame(reads_per_sample)

        #finding +-5% points & at what percentile knee point is
        """index = np.searchsorted(sorted_depths, knee_point)
        percentile = round((index / len(sorted_depths)) * 100, 2)
        lower_percentile = max(percentile - 5, 0) 
        upper_percentile = min(percentile + 5, 100) 
        lower_index = int((lower_percentile / 100) * len(sorted_depths))
        upper_index = min(int((upper_percentile / 100) * len(sorted_depths)), len(sorted_depths) - 1)
        lower_value = round(sorted_depths.iloc[lower_index])
        upper_value = round(sorted_depths.iloc[upper_index])"""

        #is this still needed?
        """reads_per_sample_df = reads_per_sample.reset_index()
        reads_per_sample_df.columns = ['sample', 'reads_per_sample']
        reads_per_sample_df["sample"] = reads_per_sample_df["sample"].astype(str)"""
        
        rps = rps.to_dataframe().reset_index()
   
    # specify names and titles according to what was run
    if kmer_run:
        graph_data = "number of observed kmers"
        graph_name = "Kmers"
        x_title = 'Sequencing Depth [Kmers]'
    else:
        graph_data = "number of observed features"
        graph_name = "Reads"
        x_title = 'Sequencing Depth [Reads]'

    if beta:
        graph_data = "Distance"

    TEMPLATES = os.path.join(
        os.path.dirname(__file__),
        "assets"
    ) 

    #dynamically populating the vega plot
    if beta:
        #new mm plot
        with open(os.path.join(TEMPLATES, "mm-beta.json")) as f:
            spec_beta = json.load(f)
        for d in spec_beta["data"]:
            if d["name"] == "raw":
                d["values"] = line_chart_df.to_dict(orient='records')
            if d["name"] == "samples":
                d["values"] = df_bars.to_dict(orient='records') 
            if d["name"] == "knee_points":
                d["values"] = kp_list_beta
        
        for signal in spec_beta["signals"]:
            if signal["name"] == "metricField":
                signal["bind"]["options"] = beta_metrics

    if alpha:
        with open(os.path.join(TEMPLATES, "mm_alpha_div.json")) as f:
            spec = json.load(f)

        for signal in spec["signals"]:
            if signal["name"] == "groupField":
                signal["bind"]["options"] = metadata_columns
            if signal["name"] == "metricField":
                signal["bind"]["options"] = alpha_metrics
            if signal["name"] == "x_axis_title":
                signal["value"] = x_title

        for d in spec["data"]:
            if d["name"] == "raw":
                combined = combined.to_dataframe().reset_index()
                combined = combined.drop('id', axis=1)
                d["values"] = combined.to_dict(orient='records')
            if d["name"] == "samples":
                max_range = np.linspace(1, max_reads, num=steps, dtype=int)
                depths_list = [int(d) for d in max_range]
                rps.rename(columns={"sample-id": "sample"}, inplace=True)
                rps = rps.set_index("sample").reset_index()
                samples_records = rps.to_dict(orient='records')
                for s in samples_records:
                    s["all_depths"] = depths_list
                d["values"] = samples_records
                """rps.rename(columns={"sample-id": "sample"}, inplace=True)
                rps = rps.set_index("sample").reset_index()
                d["values"] = rps.to_dict(orient='records')"""
            if d["name"] == "knee_points":
                d["values"] = kp_list

    
    vega_json = json.dumps(spec)
    vega_json2 = json.dumps(spec_beta)

    tabbed_context = {
        "tabs": [
            {"title": "Alpha Diversity", "url": "index.html"},
            {"title": "Beta Diversity", "url": "stats.html"},
        ],
        "vega_json": vega_json,
        "vega_json2": vega_json2,
        "beta_metric": str(metric),
        "algorithm": str(algorithm),
        "knee_point": (knee_point),
        "beta": beta,
        "percent_samples_100": json.dumps(float(percent_samples_100)),
        "depth_threshold": json.dumps(int(depth_threshold)),
        "percentile": json.dumps(float(percentile)),
        "lower_percentile": json.dumps(float(lower_percentile)),
        "upper_percentile": json.dumps(float(upper_percentile)),
        "lower_value": json.dumps(float(lower_value)),
        "upper_value": json.dumps(float(upper_value)),
        "graph_data": graph_data,
        "graph_name": graph_name,
        "max_read_percentile": json.dumps(int(max_read_percentile))
    }
    
    #templates = os.path.join(TEMPLATES, 'index.html')
    # new version with tabs
    templates = [
        os.path.join(TEMPLATES, 'index.html'),
        os.path.join(TEMPLATES, 'stats.html')
    ]

    copytree(
        src=TEMPLATES,
        dst=output_dir,
        dirs_exist_ok=True 
    )
   
    q2templates.render(templates, output_dir, context=tabbed_context)


    
def _rf_knee_locator(artifacts_list: pd.DataFrame, sample_list: list[str], steps: int, algorithm: str, max_reads: int) -> tuple[pd.DataFrame, int]:

    knee_points = [None] * len(sample_list)
    df_list = []
    max_range = np.linspace(1, max_reads, num=steps, dtype=int)

    for i, sample in enumerate(sample_list):
        array_sample = np.array(artifacts_list.iloc[i]).flatten() 
        sample = np.full(len(max_range), sample)  

        sample_df = pd.DataFrame({'depth': max_range, 'observed_features': array_sample, 'sample': sample})
        df_list.append(sample_df)

        knee_points[i] = knee_point_locator(max_range, array_sample, algorithm, "alpha")

    combined_df = pd.concat(df_list, ignore_index=True)
    
    knee_points_filtered = [point for point in knee_points if point is not None] 
    knee_point = round(np.mean(knee_points_filtered))
    print("calculated rarefaction depth:")
    print(knee_point)

    return combined_df, knee_point

    

def _rf_visualizer_boots(output_dir: str, sample_names: list[str], percent_samples: float, reads_per_sample: list[int], sorted_depths: list[int],
                         metric: str, max_reads: int, depth_threshold: int, knee_point: int, kmer_run: bool, combined_df: pd.DataFrame,
                         max_read_percentile: int, algorithm: str)-> None: 
    
    combined_df['sample'] = sample_names
    sorted_depths = pd.Series(sorted_depths)
    reads_per_sample = pd.DataFrame(reads_per_sample)

    #finding +-5% points & at what percentile knee point is
    index = np.searchsorted(sorted_depths, knee_point)
    percentile = round((index / len(sorted_depths)) * 100, 2)

    lower_percentile = max(percentile - 5, 0) 
    upper_percentile = min(percentile + 5, 100) 

    lower_index = int((lower_percentile / 100) * len(sorted_depths))
    upper_index = min(int((upper_percentile / 100) * len(sorted_depths)), len(sorted_depths) - 1)

    lower_value = round(sorted_depths.iloc[lower_index])
    upper_value = round(sorted_depths.iloc[upper_index])

    percent_samples_100 = round(percent_samples * 100, 2)

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
        bind=alt.binding_checkbox(name='Show read depth specified by slider as a line on the plot: '))

    s = alt.param(
        name='position', bind=alt.binding_range(min=0, max=max_reads, step=20, name='Rarefaction Depth Line'), value=knee_point)

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
    ).add_params(zoom, param_checkbox).properties(title='Rarefaction Curves')

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
    ).add_params(s).transform_calculate(position='position').transform_filter(param_checkbox
    ).properties(
        title='Interactive Vertical Line',
        width=450, 
        height=350)
    
    final_with_line = alt.layer(final_chart, vertical_line).resolve_scale(x='shared', y='shared')

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

    barplot_combined = alt.layer(background, barplot).resolve_scale(x='shared', y='shared')
    
    barplot_combined = alt.vconcat(barplot_combined, text_lines).properties(spacing=0)
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
        "max_read_percentile": json.dumps(int(max_read_percentile)),
        "beta": False,
        "beta_metric": metric,
        "algorithm": str(algorithm)
    }
    
    templates = os.path.join(TEMPLATES, 'index.html')

    copytree(
        src=TEMPLATES,
        dst=output_dir,
        dirs_exist_ok=True 
    )
   
    q2templates.render(templates, output_dir, context=tabbed_context)
    