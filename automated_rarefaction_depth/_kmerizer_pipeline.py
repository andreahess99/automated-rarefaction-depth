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
#import q2templates
from qiime2 import Artifact, sdk
#from kneed import KneeLocator
#import altair as alt
import time
#import tracemalloc
#from shutil import copytree
import shutil
from bs4 import BeautifulSoup
#from tempfile import TemporaryDirectory
#from qiime2.plugin import Visualization
#from q2_types.feature_table import FeatureTable, Frequency
#from q2_types.sample_data import AlphaDiversity, SampleData
import warnings
from biom import Table



def df_to_feature_table(df: pd.DataFrame) -> qiime2.Artifact:
    # Convert DataFrame to BIOM format
    biom_table = Table(df.values, observation_ids=df.index, sample_ids=df.columns)
    
    # Convert BIOM table to QIIME 2 FeatureTable[Frequency] Artifact
    feature_table_artifact = qiime2.Artifact.import_data("FeatureTable[Frequency]", biom_table)
    
    return feature_table_artifact  

def change_html_file(file_path: str) -> None:
    #add css stylesheet to the html file
    with open(file_path, 'r') as file:
        soup = BeautifulSoup(file, 'html.parser')
        link_tag = soup.new_tag('link', rel='stylesheet', href='./css/styles.css')
        soup.head.append(link_tag)

    with open(file_path, 'w') as file:
        file.write(str(soup))


_pipe_defaults = {
    'iterations': 1, #10
    'table_size': None,
    'steps': 20, #20
    'percent_samples': 0.8,
    'algorithm': 'kneedle'
}

warnings.simplefilter(action='ignore', category=FutureWarning)
 
#boots and diversity pipeline
#change so that the same depths are used for all samples -> only run boots alpha once
#linearly space according to the highest number of reads a sample has
def pipeline_kmerizer(ctx, table, sequence, metadata, iterations=_pipe_defaults['iterations'], table_size=_pipe_defaults['table_size'],
                   steps=_pipe_defaults['steps'], percent_samples=_pipe_defaults['percent_samples'], algorithm=_pipe_defaults['algorithm']):
    start_time = time.time()
    #alpha_action = ctx.get_action('boots', 'alpha')
    kmerizer_action = ctx.get_action('kmerizer', 'core_metrics')
    #change this to this files visualizer??
    viz_action = ctx.get_action('rarefaction-depth', '_rf_visualizer_boots')

    table_df = table.view(pd.DataFrame)
    print(len(table_df))

    if table_df.empty:
        raise ValueError("The feature table is empty.")
    if not np.issubdtype(table_df.values.dtype, np.number):
        raise ValueError("The feature table contains non-numerical values.")
    

    #adjusting table size if it's too big -> keep table_size rows
    if (table_size is not None and len(table_df) > table_size):
        table_df = table_df.sample(n=table_size, random_state=42) 
        table_df = table_df.loc[:, ~(table_df.isna() | (table_df == 0)).all(axis=0)] 
    table_size = len(table_df)

    reads_per_sample = table_df.sum(axis=1) 
    #max_reads = reads_per_sample.max()
    percentile_90 = int(np.percentile(reads_per_sample, 90))
    print(f"90th percentile of reads per sample: {percentile_90}")
    max_reads = percentile_90
    

    sorted_depths = reads_per_sample.sort_values() 
    print("sorted_depths:")
    print(sorted_depths)

    sorted_depths_pass = sorted_depths.tolist()
    sorted_depths_pass = [int(depth) for depth in sorted_depths_pass]
    reads_per_sample_pass = reads_per_sample.tolist()
    reads_per_sample_pass = [int(read) for read in reads_per_sample_pass]

    sample_loss_index = int(np.ceil((1-percent_samples) * len(sorted_depths))) - 1 
    depth_threshold = int(sorted_depths.iloc[sample_loss_index])

    artifacts_list = []
    sample_list = table_df.index.tolist()

    max_range = np.linspace(1, max_reads, num=steps, dtype=int)
    table_artifact = Artifact.import_data('FeatureTable[Frequency]', table_df)

    for i in range(steps):
        print(f"step: {max_range[i]}")   
        result = kmerizer_action(table=table_artifact, sampling_depth=int(max_range[i]), sequences=sequence, with_replacement=False, metadata=metadata)[2]
        artifacts_list.append(result)

    pd_new = pd.DataFrame(
        np.nan,  
        index=[sample_list[i] for i in range(table_size)],  
        columns=[f"Step_{j+1}" for j in range(steps)]  
    )
    
    #trying new approach because nothing else works
    for i, artifact in enumerate(artifacts_list):
        artifact.export_data(f'output_directory_{i}')  # Export each artifact
    dfs = []
    for i in range(len(artifacts_list)):
        df = pd.read_csv(f'output_directory_{i}/alpha-diversity.tsv', sep='\t', index_col=0)
        dfs.append(df)
        for j, sample in enumerate(sample_list):
            if sample in df.index:
                pd_new.iloc[j, i] = df.loc[sample].values[0].round().astype(int)
            else:
                pd_new.iloc[j, i] = np.nan

    final_df = pd.concat(dfs, axis=0)
    final_df = final_df.reset_index(drop=True)
  
    artifact_ft = df_to_feature_table(pd_new)

    visualization, = viz_action( percent_samples=percent_samples, reads_per_sample=reads_per_sample_pass, artifacts_list=artifact_ft, steps=int(steps),
                    sorted_depths=sorted_depths_pass, max_reads=int(max_reads), depth_threshold=int(depth_threshold), sample_list=sample_list, algorithm=algorithm)
  
    # Clean up the exported directories and files
    for i in range(len(artifacts_list)):
        dir_path = f'output_directory_{i}'
        file_path = os.path.join(dir_path, 'alpha-diversity.tsv')
        
        # Remove files & directories if they exist
        if os.path.exists(file_path):
            os.remove(file_path)
        if os.path.exists(dir_path):
            shutil.rmtree(dir_path)

    #measuring time
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Pipeline execution time: {elapsed_time:.2f} seconds")

    return visualization
    
