# ----------------------------------------------------------------------------
# Copyright (c) 2024, Andrea Hess.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import automated_rarefaction_depth
from qiime2.plugin import (Plugin, Str, Properties, Choices, Int, Bool, Range,
                           Float, Set, Visualization, Metadata, MetadataColumn,
                           Categorical, Numeric, Citations, Threads, List)
from q2_types.feature_table import FeatureTable, Frequency
from q2_types.feature_data import (FeatureData, Sequence, RNASequence, ProteinSequence)
from automated_rarefaction_depth import __version__
from automated_rarefaction_depth._boots_pipeline import pipeline_boots, _rf_visualizer_boots


citations = Citations.load("citations.bib", package="automated_rarefaction_depth")

plugin = Plugin(
    name="rarefaction-depth",
    version=automated_rarefaction_depth.__version__,
    website="https://github.com/andreahess99/automated-rarefaction-depth",
    package="automated_rarefaction_depth",
    description="This qiime2 plugin gives you a tool to automatically calculate the ideal rarefaction depth based on the given data and some user parameters.",
    short_description="A tool that automatically calculates the ideal rarefaction depth.",
    citations=citations 
)
 


plugin.visualizers.register_function(
    function=automated_rarefaction_depth.rarefaction_depth,
    inputs={'table': FeatureTable[Frequency]},
    parameters={'seed': Int % Range(1, None),
                'percent_samples': Float % Range(0, 1),
                'iterations': Int % Range(1, 100),
                'table_size': Int % Range(1, None),
                'steps': Int % Range(5, 100),
                'algorithm': Str % Choices("kneedle", "gradient")},
    input_descriptions={
        'table': ('Feature table to compute rarefaction curves from.')
    },
    parameter_descriptions={
        'seed': 'The seed used for the random subsampling.',
        'percent_samples': 'The minimal percentage of samples you want to keep, choose a decimal between 0 and 1.',
        'iterations': 'The number of times each sample gets rarefied at each depth, a positive number below 100.',
        'table_size': 'The number of samples to keep in the feature table, a positive number.',
        'steps': 'The number of points each sample is evaluated at when making the rarefaction curves, choose a number between 5 and 100.',
        'algorithm': 'The algorithm to use for the rarefaction depth calculation, either kneedle or gradient.'
    },
    name='Automated Rarefaction Depth',
    description=("Automatically computes an optimal rarefaction depth. Outputs a visualization with the rarefaction curves, the optimal rarefaction depth and a histogram of the reads per sample."),
    citations=citations,
)

# descriptions of parameters used in the kmerizer action are copied from the q2-kmerizer plugin

plugin.pipelines.register_function(
    function=pipeline_boots,
    inputs={'table': FeatureTable[Frequency],
            'sequence': FeatureData[Sequence]},
    outputs={'visualization': Visualization},
    parameters={'percent_samples': Float % Range(0, 1),
                'iterations': Int % Range(1, 100),
                'table_size': Int % Range(1, None),
                'steps': Int % Range(5, 100),
                'algorithm': Str % Choices("kneedle", "gradient"),
                'kmer_size': Int,
                'tfidf': Bool,
                'max_df': Float % Range(0, 1, inclusive_start=True,
                            inclusive_end=True) | Int,
                'min_df': Float % Range(0, 1, inclusive_start=True,
                            inclusive_end=False) | Int,
                'max_features': Int,
                'norm': Str % Choices(['None', 'l1', 'l2']) },
    parameter_descriptions={'percent_samples': 'The minimal percentage of samples you want to keep, choose a decimal between 0 and 1.',
        'iterations': 'The number of times each sample gets rarefied at each depth, a positive number below 100.',
        'table_size': 'The number of samples to keep in the feature table, a positive number.',
        'steps': 'The number of depths that get evaluated between the minimum and maximum sample depth, choose a number between 5 and 100.',
        'algorithm': 'The algorithm to use for the rarefaction depth calculation, either kneedle or gradient.',
        'kmer_size': 'Only needed for kmerizer! Length of kmers to generate.',
        'tfidf': 'Only needed for kmerizer! If True, kmers will be scored using TF-IDF and output '
             'frequencies will be weighted by scores. If False, kmers are '
             'counted without TF-IDF scores.',
        'max_df': 'Only needed for kmerizer! Ignore kmers that have a frequency strictly higher than '
              'the given threshold. If float, the parameter represents a '
              'proportion of sequences, if an integer it represents an '
              'absolute count.',
        'min_df': 'Only needed for kmerizer! Ignore kmers that have a frequency strictly lower than '
              'the given threshold. If float, the parameter represents a '
              'proportion of sequences, if an integer it represents an '
              'absolute count.',
        'max_features': 'Only needed for kmerizer! If not None, build a vocabulary that only considers '
                    'the top max_features ordered by frequency (or TF-IDF '
                    'score).',
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
    description=("Automatically computes an optimal rarefaction depth using q2-boots. If sequences are provided, "
                 "kmerizer is used to generate a feature table from the sequences, which will be used for the rarefaction depth calculation. "),
    citations=citations,
)

plugin.visualizers.register_function(
    function=_rf_visualizer_boots,
    inputs={'artifacts_list': FeatureTable[Frequency]},
    parameters={'sorted_depths': List[Int],
                'percent_samples': Float % Range(0, 1),
                'reads_per_sample': List[Int],
                'max_reads': Int % Range(1, None),
                'depth_threshold': Int % Range(1, None),
                'sample_list': List[Str],
                'steps': Int % Range(5, 100),
                'algorithm': Str % Choices("kneedle", "gradient"),
                'kmer_run': Bool
                },
    input_descriptions={
        'artifacts_list': 'A set of .qza files containing AlphaDiversity artifacts.'
    },
    parameter_descriptions={
        'sorted_depths': 'A list of sorted depths as integers.',
        'percent_samples': 'The minimal percentage of samples you want to keep, choose a decimal between 0 and 1.',
        'reads_per_sample': 'A list of how many reads each sample has.',
        'max_reads': 'The maximum amount of reads a single sample has.',
        'depth_threshold': 'The highest read_depth to still be within the accepted area.',
        'sample_list': 'The list of samples in the order they are processed.', 
        'steps': 'The number of depths that get evaluated for each sample.',
        'algorithm': 'The algorithm to use for the rarefaction depth calculation, either kneedle or gradient.',
        'kmer_run': 'True if the pipeline was run with kmerizer, False otherwise.'
    },
    name='Automated Rarefaction Depth',
    description=("Automatically computes an optimal rarefaction depth."),
    citations=citations,
)

