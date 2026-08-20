# ----------------------------------------------------------------------------
# Copyright (c) 2016-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import json
import os
import base64 
#import shutil
from shutil import copytree, copy
from xxlimited import Str
import tempfile, zipfile
import numpy as np
import pandas as pd
import qiime2
import q2templates
from kneed import KneeLocator

import warnings
from skbio import DistanceMatrix


_pipe_defaults = {
    'iterations': 10,
    'max_samples': None,
    'steps': 20,
    'algorithm': 'kneedle', 
    'seed': 42,
    'kmer_size': 16,
    'tfidf': False,
    'max_df': 1.0,
    'min_df': 1,
    'max_features': None,
    'norm': 'None',
    'metrics': {'observed_features', 'shannon'},
    'max_depth': None
}

_DEFAULT_MAX_DEPTH_PERCENTILE = 90
_alpha_metrics = ['observed_features', 'shannon', 'simpson', 'brillouin_d', 'chao1', 'enspie', 'goods_coverage', 'michaelis_menten_fit',
                  'dominance', 'robbins', 'simpson_e', 'mcintosh_e', 'berger_parker_d', 'jaccard', 'braycurtis']
_beta_metrics = ['braycurtis', 'jaccard', 'hamming', 'dice', 'jensenshannon', 'matching', 'rogerstanimoto', 'russellrao', 'canberra_adkins',
                 'sokalmichener', 'sokalsneath', 'yule', 'correlation', 'cosine', 'aitchison',  'canberra']


#ignore future warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=RuntimeWarning)
warnings.simplefilter(action='ignore', category=UserWarning)

def pipeline_boots(ctx, table, metadata, sequence=None, iterations=_pipe_defaults['iterations'], max_samples=_pipe_defaults['max_samples'], metrics=_pipe_defaults['metrics'],
                   steps=_pipe_defaults['steps'], algorithm=_pipe_defaults['algorithm'], max_depth=_pipe_defaults['max_depth'],
                   seed = _pipe_defaults['seed'], kmer_size=_pipe_defaults['kmer_size'], tfidf=_pipe_defaults['tfidf'], max_df=_pipe_defaults['max_df'],
                   min_df=_pipe_defaults['min_df'], max_features=_pipe_defaults['max_features'], norm=_pipe_defaults['norm']):
    
    beta_action = ctx.get_action('boots', 'beta')
    kmer_action = ctx.get_action('kmerizer', 'seqs_to_kmers')
    viz_combined_action = ctx.get_action('rarefaction-depth', '_combined_viz')
    alpha_collection_action = ctx.get_action("boots", "alpha_collection")

    table_df = table.view(pd.DataFrame)
    alpha, metrics_alpha = False, []
    beta, metrics_beta = False, []

    if any(m in _beta_metrics for m in metrics):
        beta = True
        metrics_beta = [m for m in metrics if m in _beta_metrics]
    else:
        kp_list_beta = None
        df_bars = None
        num_samples = None
        data_beta = None
        
    if any(m in _alpha_metrics for m in metrics):
        alpha = True 
        metrics_alpha = [m for m in metrics if m not in metrics_beta]

    # 2. instead of converting into DF and filtering, you could use the "filter_columns" method
    # of the Metadata object directly to only keep the categorical columns, see here:
    # https://develop.qiime2.org/en/stable/plugins/references/metadata-api.html#the-qiime-metadata-class
    meta = metadata.to_dataframe()
    meta.index.name = "sample"
    # find all numeric metadata columns
    numeric_columns = meta.select_dtypes(include=[np.number]).columns.tolist()
    meta = meta.drop(columns=numeric_columns)
    metadata_columns = ["sample"] + meta.columns.tolist()
    
    #run seqs_to_kmers if sequence is provided
    kmer_run = False
    if sequence is not None:
        print("Sequences were provided as input.")
        print("Therefore, the kmerizer is run to generate a new feature table with kmer sequences.")
        print("This new kmer feature table will be used for the analysis.")
        table, = kmer_action(table=table, sequences=sequence, kmer_size=kmer_size, tfidf=tfidf, max_df=max_df, min_df=min_df, max_features=max_features, norm=norm)
        table_df = table.view(pd.DataFrame)
        kmer_run = True

   
    #adjusting table size if it's too big -> keep max_samples rows
    if (max_samples is not None and len(table_df) > max_samples):
        table_df = table_df.loc[:, ~(table_df.isna() | (table_df == 0)).all(axis=0)] 
        table_df = table_df.sample(n=max_samples, random_state=seed)
        print(f"Table size is larger than {max_samples}. Randomly subsampling to {max_samples} samples.")


    reads_per_sample = table_df.sum(axis=1)

    #round the max_depth to a nice number
    if max_depth is None:
        max_reads = int(np.percentile(reads_per_sample, _DEFAULT_MAX_DEPTH_PERCENTILE))
        if max_reads < 2000:
            max_reads = (int(max_reads / 100) * 100)
        else:
            max_reads = (int(max_reads / 500) * 500)
    else:
        max_reads = max_depth
   
    reads_per_sample_merged = pd.DataFrame(
        {"reads": reads_per_sample.astype(int)},
        index=table_df.index.rename("sample-id"),
    ).join(meta, how="left")
   
    sample_list = table_df.index.tolist()

    # I would move this to right under the max_range calculation
    # also, you should jsut do these two steps in one go and use only one variable
    max_range = np.linspace(1, max_reads, num=steps, dtype=int)
    clean_max_range = [float(x) for x in max_range]

    # this should not be necessary: you should not be zipping these artifacts and passing
    # them to the visualziation action; instead you should pass the list of artifacts directly
    temp_zip_path = None 

    # beta metric specific code
    # this entire section should move to a separate function
    with tempfile.TemporaryDirectory() as temp_dir:
        if beta:
            knee_points_beta = []
            data_beta = []
            # this whole loop could be simplified to somehting like this:
            for k, metric in enumerate(metrics_beta):
                print("metric:", metric)
                avg_diffs, avg_ranges, num_samples_left = [], [], []
                old_dm = None

                for i, depth in enumerate(max_range):
                    depth_int = int(depth)
                    print(f"step {i + 1}: {depth_int}")

                    beta_result, = beta_action(
                        table=table,
                        sampling_depth=depth_int,
                        metric=metric,
                        n=iterations,
                        replacement=False,
                    )

                    dm = beta_result.view(DistanceMatrix)
                    if k == 0:
                        num_samples_left.append(dm.shape[0])

                    if old_dm is not None:
                        # Filter and align old distance matrix to current IDs
                        aligned_old_dm = old_dm.filter(dm.ids)
                        diff = np.abs(dm.data - aligned_old_dm.data)

                        upper_vals = diff[np.triu_indices_from(diff, k=1)]
                        avg_diffs.append(
                            float(np.mean(upper_vals)) if upper_vals.size else 0.0
                        )
                        avg_ranges.append((max_range[i - 1] + depth) / 2)

                    old_dm = dm

                # Drop the first delta (step 1 vs 0) to match your original [1:] slicing
                eval_ranges = avg_ranges[1:]
                eval_diffs = avg_diffs[1:]

                data_beta.append(
                    pd.DataFrame(
                        {"metric": metric, "depth": eval_ranges, "observed": eval_diffs}
                    )
                )
                if k == 0:
                    df_bars = pd.DataFrame(
                        {"depth": eval_ranges, "num_samples_left": num_samples_left[2:]}
                    )

                kpb = knee_point_locator(
                    eval_ranges, eval_diffs, algorithm, "convex", "decreasing"
                )
                kpb = round(float(kpb)) if kpb is not None else 0
                knee_points_beta.append(pd.DataFrame({"knee": [kpb], "metric": [metric]}))
                print("calculated rarefaction depth:", kpb)

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

            # make the beta artifacts zip file
            # temp_zip_path = os.path.join(temp_dir, 'beta_matrices.zip')
            # with zipfile.ZipFile(temp_zip_path, 'w', zipfile.ZIP_DEFLATED) as zipf:
            #     for file_path in beta_artifact_paths:
            #         zipf.write(file_path, arcname=os.path.basename(file_path))
        else:
            metrics_beta = None

        # this entire section should move to a separate function
        if alpha:
            #if alpha metric was chosen
            combined_dfs = []
            knee_point_list = []

            for metric in metrics_alpha:
                print("metric:", metric)
                dfs = []
                for i in range(steps):
                    print(f"step {i+1}: {max_range[i]}")

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

                    # Calculate average after collecting all iterations for this depth
                    mean_df = (combined.groupby(['sample', 'read_depth'], as_index=False).agg(mean_observed=('observed', 'mean')))
                
                num_cols = meta.select_dtypes(include=[np.number]).columns 
                meta[num_cols] = meta[num_cols].fillna(0)

                combined = combined.merge(meta, left_on="sample", right_index=True, how="left")
                combined_dfs.append(combined)

                #calculating knee point 
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
                if len(knee_points_filtered) > 0:
                    knee_point = round(np.mean(knee_points_filtered))
                else:
                    #default value if no sample yields a valid knee point
                    print(f"Warning: No valid knee points found for metric {metric}. Defaulting to 0.")
                    knee_point = 0

                print("calculated rarefaction depth:", knee_point)
                knee_point_list.append((knee_point, metric))
    
            combined = pd.concat(combined_dfs, ignore_index=True)
            combined.insert(0, 'id', [f"row{i}" for i in range(len(combined))])
            combined = combined.set_index('id')
            combined = qiime2.Metadata(combined)
            
            kp_df = pd.DataFrame(knee_point_list, columns=['knee', 'metric'])
            kp_df.index.name = 'id'
            kp_df.index = kp_df.index.astype(str)
            knee_point_list = qiime2.Metadata(kp_df)
            metrics = list(metrics)
        else:
            combined = None
            knee_point_list = None
            metrics_alpha = None
            knee_point = 0
            
        visualization, = viz_combined_action(kmer_run=kmer_run, steps=steps, algorithm=algorithm, kp_list=knee_point_list,
                            max_reads=int(max_reads), max_range=clean_max_range, alpha_metrics=metrics_alpha, numeric_columns=numeric_columns, 
                            combined=combined, metadata_columns=metadata_columns, rps=qiime2.Metadata(reads_per_sample_merged), 
                            kp_list_beta=kp_list_beta, data_beta=data_beta, beta_metrics=metrics_beta, num_samples=num_samples, beta_zip_path=temp_zip_path) 
                            
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


# combined visualization function for alpha and beta metrics
def _combined_viz(output_dir: str, kmer_run: bool, max_range: list[float] = None, algorithm: str = "kneedle", num_samples: qiime2.Metadata = None, steps: int = None, max_reads: int = 1,
                  metadata_columns: list[str] = None, combined: qiime2.Metadata = None, rps:qiime2.Metadata = None, kp_list: qiime2.Metadata = None, beta_zip_path: Str = None,
                  kp_list_beta: qiime2.Metadata = None, data_beta: qiime2.Metadata = None, alpha_metrics: list[str] = None, beta_metrics: list[str] = None, numeric_columns: list[str] = None)->None:  
    
    #default values for the tabbed_context
    beta = False
    alpha = False
    zero_beta_metrics = []

    # beta metric specific code
    if beta_metrics is not None and len(beta_metrics) > 0:
        beta = True

        line_chart_df = data_beta.to_dataframe().reset_index()
        line_chart_df = line_chart_df.drop('id', axis=1)
        kp_list_beta = kp_list_beta.to_dataframe().reset_index()
        kp_list_beta = kp_list_beta.drop('id', axis=1)
        kp_list_beta = kp_list_beta.replace({np.nan: None})

        kp_list_beta = kp_list_beta.to_dict(orient='records')
        df_bars = num_samples.to_dataframe().reset_index()
        df_bars = df_bars.drop('id', axis=1)
        df_bars = df_bars.replace({np.nan: None})


    #alpha metric specific code
    if alpha_metrics is not None and len(alpha_metrics) > 0:
        alpha = True
        kp_list = kp_list.to_dataframe().reset_index()
        kp_list = kp_list.drop('id', axis=1)
        kp_list = kp_list.replace({np.nan: None})

        kp_list = kp_list.to_dict(orient='records')
        rps = rps.to_dataframe().reset_index()
   
    if kmer_run:
        x_title = 'Sequencing Depth [Kmers]'
    else:
        x_title = 'Sequencing Depth [Reads]'

    TEMPLATES = os.path.join(
        os.path.dirname(__file__),
        "assets"
    ) 

    #dynamically populating the vega plot
    if beta:
        with open(os.path.join(TEMPLATES, "beta_vega.json")) as f:
            spec_beta = json.load(f)
        for d in spec_beta["data"]:
            if d["name"] == "samples":
                d["values"] = df_bars.to_dict(orient='records') 
            if d["name"] == "knee_points":
                d["values"] = kp_list_beta
        
        for signal in spec_beta["signals"]:
            if signal["name"] == "metricField":
                signal["bind"]["options"] = beta_metrics

        csv_string_beta = line_chart_df.to_csv(index=True)
        for metric in kp_list_beta:
            if metric['knee'] == 0:
                zero_beta_metrics.append(metric['metric'])
     
    else:
        spec_beta = {"warning": "Warning! No beta metric was specified!"}
        csv_string_beta = ""
        

    if alpha:
        with open(os.path.join(TEMPLATES, "alpha_vega.json")) as f:
            spec = json.load(f)

        for signal in spec["signals"]:
            if signal["name"] == "groupField":
                signal["bind"]["options"] = metadata_columns
            if signal["name"] == "metricField":
                signal["bind"]["options"] = alpha_metrics
            if signal["name"] == "x_axis_title":
                signal["value"] = x_title

        for d in spec["data"]:
            if d["name"] == "samples":
                max_range = np.linspace(1, max_reads, num=steps, dtype=int)
                depths_list = [int(d) for d in max_range]
                rps.rename(columns={"sample-id": "sample"}, inplace=True)
                rps = rps.set_index("sample").reset_index()
                rps = rps.replace({np.nan: None})

                samples_records = rps.to_dict(orient='records')
                for s in samples_records:
                    s["all_depths"] = depths_list
                d["values"] = samples_records
            if d["name"] == "knee_points":
                d["values"] = kp_list

        combined = combined.to_dataframe().reset_index()
        combined = combined.drop('id', axis=1)
        csv_string_alpha = combined.to_csv(index=True)
       
    else:
        spec = {"warning": "Warning! No alpha metric was specified!"}
        csv_string_alpha = ""

    vega_json = json.dumps(spec)
    vega_json2 = json.dumps(spec_beta)

    #processing the distance matrix zip file
    beta_zip_base64 = ""
    if beta_zip_path and os.path.exists(beta_zip_path):
        destination = os.path.join(output_dir, 'beta_matrices.zip')
        copy(beta_zip_path, destination)
        with open(beta_zip_path, "rb") as zip_file:
            encoded_bytes = base64.b64encode(zip_file.read())
            beta_zip_base64 = encoded_bytes.decode('utf-8')

    tabbed_context = {
        "tabs": [
            {"title": "Alpha Rarefaction", "url": "index.html"},
            {"title": "Beta Rarefaction", "url": "index_beta.html"},
        ],
        "vega_json": vega_json,
        "vega_json2": vega_json2,
        "zero_beta_metrics": zero_beta_metrics,
        "algorithm": str(algorithm),
        "alpha": alpha,
        "beta": beta,
        "csv_data_alpha": csv_string_alpha,
        "csv_data_beta": csv_string_beta,
        "removed_numeric_columns": numeric_columns,
        'beta_zip_base64': beta_zip_base64
    }

    templates = [
        os.path.join(TEMPLATES, 'index.html'),
        os.path.join(TEMPLATES, 'index_beta.html')
    ]

    copytree(
        src=TEMPLATES,
        dst=output_dir,
        dirs_exist_ok=True 
    )

    q2templates.render(templates, output_dir, context=tabbed_context)

    