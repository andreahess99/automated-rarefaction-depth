# ----------------------------------------------------------------------------
# Copyright (c) 2024, Andrea Hess.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import automated_rarefaction_depth
from qiime2.plugin import (Plugin, Str, Choices, Int, Bool, Range, Float, List,
                            Set, Visualization, Metadata, Citations)
from q2_types.feature_table import FeatureTable, Frequency
from q2_types.feature_data import FeatureData, Sequence
#from q2_types.distance_matrix import DistanceMatrix
from automated_rarefaction_depth import __version__
from automated_rarefaction_depth._boots_pipeline import pipeline_boots,  _combined_viz


citations = Citations.load("citations.bib", package="automated_rarefaction_depth")

plugin = Plugin(
    name="rarefaction-depth",
    version=automated_rarefaction_depth.__version__,
    website="https://github.com/andreahess99/automated-rarefaction-depth",
    package="automated_rarefaction_depth",
    description="This qiime2 plugin gives you a tool to automatically calculate the ideal rarefaction depth based on the given data and some user parameters.",
    short_description="Plugin for automatically calculating the ideal rarefaction depth.",
    citations=citations 
)
 

# descriptions of parameters used in the kmerizer action are copied from the q2-kmerizer plugin
plugin.pipelines.register_function(
    function=pipeline_boots,
    inputs={'table': FeatureTable[Frequency],
            'sequence': FeatureData[Sequence]},
    outputs={'visualization': Visualization},
    parameters={'metadata': Metadata,
                'iterations': Int % Range(1, 100),
                'max_samples': Int % Range(1, None),
                'steps': Int % Range(10, 100), 
                'algorithm': Str % Choices("kneedle", "gradient"),
                'seed': Int % Range(1, None),
                'metrics': Set[Str % Choices(['observed_features', 'shannon', 'braycurtis', 'jaccard', 'simpson', 'brillouin_d', 
                                         'chao1', 'enspie', 'goods_coverage', 'michaelis_menten_fit','dominance', 'simpson_e', 'mcintosh_e',
                                         'robbins', 'berger_parker_d','hamming', 'dice', 'correlation', 'sokalmichener', 'yule',
                                         'jensenshannon', 'matching', 'rogerstanimoto', 'russellrao', 'sokalsneath',
                                         'canberra_adkins', 'cosine', 'aitchison',  'canberra'])], 
                'kmer_size': Int,
                'max_depth': Int,
                'tfidf': Bool,
                'max_df': Float % Range(0, 1, inclusive_start=True,
                            inclusive_end=True) | Int,
                'min_df': Float % Range(0, 1, inclusive_start=True,
                            inclusive_end=False) | Int,
                'max_features': Int,
                'norm': Str % Choices(['None', 'l1', 'l2']) },
    parameter_descriptions={'iterations': 'The number of times each sample gets rarefied at each depth, a positive number below 100.',
        'max_samples': 'The maximum number of samples to keep in the feature table, a positive number.',
        'steps': 'The number of depths that get evaluated between the minimum and maximum sample depth, choose a number between 5 and 100.',
        'algorithm': 'The algorithm to use for the rarefaction depth calculation, either kneedle or gradient.',
        'seed': 'The seed used for the random sampling of samples in case the table is larger than the max_samples parameter. A positive integer.',
        'max_depth': 'The maximum depth to evaluate for the rarefaction curves. If None, the maximum depth will be determined automatically.',
        'metrics': 'The different alpha and beta diversity metrics to use for the rarefaction curves. The available metrics are: observed_features, shannon, '
                'braycurtis, jaccard, simpson, brillouin_d, chao1, enspie, goods_coverage, michaelis_menten_fit, dominance, simpson_e, mcintosh_e, robbins, '
                'berger_parker_d, hamming, dice, correlation, sokalmichener, yule, jensenshannon, matching, rogerstanimoto, russellrao, sokalsneath, canberra_adkins, cosine, aitchison and canberra.',
        'kmer_size': 'Only needed for kmerizer! Length of kmers to generate.',
        'tfidf': 'Only needed for kmerizer! If True, kmers will be scored using TF-IDF and output '
             'frequencies will be weighted by scores. If False, kmers are counted without TF-IDF scores.',
        'max_df': 'Only needed for kmerizer! Ignore kmers that have a frequency strictly higher than '
              'the given threshold. If float, the parameter represents a '
              'proportion of sequences, if an integer it represents an absolute count.',
        'min_df': 'Only needed for kmerizer! Ignore kmers that have a frequency strictly lower than '
              'the given threshold. If float, the parameter represents a '
              'proportion of sequences, if an integer it represents an absolute count.',
        'max_features': 'Only needed for kmerizer! If not None, build a vocabulary that only considers '
                    'the top max_features ordered by frequency (or TF-IDF score).',
        'norm': 'Only needed for kmerizer! Normalization procedure applied to TF-IDF scores. Ignored '
            'if tfidf=False. l2: Sum of squares of vector elements is 1. '
            'l1: Sum of absolute values of vector elements is 1.'},
    input_descriptions={
        'table': 'Feature table to compute rarefaction curves from.',
        'sequence': 'FeatureData containing sequences to be used in the kmerizer. If not provided, the feature table will be used directly.'
    },
    output_descriptions={
        'visualization': 'Visualization of the optimal rarefaction depth.'
    },
    name='Automated Rarefaction Depth Pipeline',
    description=("Automatically computes an optimal rarefaction depth using q2-boots. If sequences are provided, the "
                 "kmerizer is used to generate a feature table from the sequences, which will be used for the rarefaction depth calculation. "),
    citations=citations,
)

plugin.visualizers.register_function(
    function=_combined_viz,
    inputs={},
    parameters={'combined': Metadata,
                'steps': Int,
                'max_reads': Int % Range(1, None),
                'kmer_run': Bool,
                'max_range': List[Float],
                'algorithm': Str % Choices(['kneedle', 'gradient']),
                'num_samples': Metadata,
                'metadata_columns': List[Str],
                'rps': Metadata,
                'alpha_metrics': List[Str],
                'kp_list': Metadata,
                'data_beta': Metadata,
                'kp_list_beta': Metadata,
                'beta_metrics': List[Str],
                'numeric_columns': List[Str],
                'beta_zip_path': Str
                },
    input_descriptions={},
    parameter_descriptions={
        'steps': 'The number of depths that get evaluated between the minimum and maximum sample depth.',
        'max_reads': 'The maximum amount of reads a single sample has.',
        'kmer_run': 'True if the pipeline was run with the kmerizer, False otherwise.',
        'algorithm': "The algorithm which was chosen for the knee point calculation, kneedle or gradient",
        'max_range': 'The different read depths at which the distance matrix was calculated.',
        'num_samples': "The metadata containing the number of samples at each read depth.",
        'metadata_columns': "The metadata columns that are present in the combined feature table.",
    },
    name='Automated Rarefaction Depth',
    description=("Makes the graphs and produces the visualization for both alpha and beta diversity metrics."),
    citations=citations,
)


