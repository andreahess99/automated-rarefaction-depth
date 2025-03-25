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
                           Categorical, Numeric, Citations, Threads)
import q2_diversity
from q2_types.feature_table import FeatureTable, Frequency
from q2_types.sample_data import AlphaDiversity, SampleData
from q2_types.tree import Phylogeny, Rooted
from automated_rarefaction_depth import __version__
from automated_rarefaction_depth._pipeline import pipeline_test_new, _rf_visualizer


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
        'seed': 'The seed used for random number generation.',
        'percent_samples': 'The minimal percentage of samples you want to keep, choose a decimal between 0 and 1.',
        'iterations': 'The number of times each sample gets rarefied at each depth, a positive number below 100.',
        'table_size': 'The number of samples to keep in the feature table, a positive number.',
        'steps': 'The number of depths that get evaluated between the minimum and maximum sample depth, choose a number between 5 and 100.',
        'algorithm': 'The algorithm to use for the rarefaction depth calculation, either kneedle or gradient.'
    },
    name='Automated Rarefaction Depth',
    description=("Automatically computes an optimal rarefaction depth."),
    citations=citations,
)



plugin.pipelines.register_function(
    function=pipeline_test_new,
    inputs={'table': FeatureTable[Frequency]},
    outputs={'visualization': Visualization},
    parameters={'seed': Int % Range(1, None),
                'percent_samples': Float % Range(0, 1),
                'iterations': Int % Range(1, 100),
                'table_size': Int % Range(1, None),
                'steps': Int % Range(5, 100),
                'algorithm': Str % Choices("kneedle", "gradient")},
    parameter_descriptions={'seed': 'The seed used for random number generation.',
        'percent_samples': 'The minimal percentage of samples you want to keep, choose a decimal between 0 and 1.',
        'iterations': 'The number of times each sample gets rarefied at each depth, a positive number below 100.',
        'table_size': 'The number of samples to keep in the feature table, a positive number.',
        'steps': 'The number of depths that get evaluated between the minimum and maximum sample depth, choose a number between 5 and 100.',
        'algorithm': 'The algorithm to use for the rarefaction depth calculation, either kneedle or gradient.'},
    input_descriptions={
        'table': ('Feature table to compute rarefaction curves from.')
    },
    output_descriptions={
        'visualization': 'Visualization of the optimal rarefaction depth.'
    },
    name='Automated Rarefaction Depth Pipeline',
    description=("Automatically computes an optimal rarefaction depth using q2-boots."),
    citations=citations,
)

plugin.visualizers.register_function(
    function=_rf_visualizer,
    inputs={'artifacts_list': Set[FeatureTable[Frequency]]},#'table': FeatureTable[Frequency]
    parameters={'sorted_depths': Set[Int],
                'percent_samples': Float % Range(0, 1),
                'reads_per_sample': Set[Int],
                #'artifacts_list': Set[FeatureTable[Frequency]], #potentially wrong types for artifacts
                'max_reads': Int % Range(1, None),
                'depth_threshold': Int % Range(1, None),
                'sample_list': Set[Str],
                'depths_list': Set[Int],
                'steps': Int % Range(5, 100),
                'algorithm': Str % Choices("kneedle", "gradient")
                },
    input_descriptions={
        'artifacts_list': 'A set of .qza files containing AlphaDiversity artifacts.'#'table': ('Feature table to compute rarefaction curves from.')
    },
    parameter_descriptions={
        'sorted_depths': 'A list of sorted depths as integers.',
        'percent_samples': 'The minimal percentage of samples you want to keep, choose a decimal between 0 and 1.',
        'reads_per_sample': 'A list of how many reads each sample has.',
        #'artifacts_list': 'A set of .qza files containing AlphaDiversity artifacts.',
        'max_reads': 'The maximum amount of reads a single sample has.',
        'depth_threshold': 'The highest read_depth to still be within the accepted area.',
        'sample_list': 'The list of samples in the order they are processed.', 
        'depths_list': 'The depths used for the different samples in the order they are processed.',
        'steps': 'The number of depths that get evaluated for each sample.',
        'algorithm': 'The algorithm to use for the rarefaction depth calculation, either kneedle or gradient.'
    },
    name='Automated Rarefaction Depth',
    description=("Automatically computes an optimal rarefaction depth."),
    citations=citations,
)


