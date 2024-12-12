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
#from automated_rarefaction_depth._methods import duplicate_table


citations = Citations.load("citations.bib", package="automated_rarefaction_depth")

plugin = Plugin(
    name="rarefaction-depth",
    version=automated_rarefaction_depth.__version__,
    #website="https://example.com",
    package="automated_rarefaction_depth",
    description="This qiime2 plugin gives you a tool to automatically calculate the ideal rarefaction depth based on the given data and some user parameters.",
    short_description="A tool that automatically calculates the ideal rarefaction depth.",
    # The plugin-level citation of 'Caporaso-Bolyen-2024' is provided as
    # an example. You can replace this with citations to other references
    # in citations.bib.
    citations=[citations['Caporaso-Bolyen-2024']]
)
 
plugin.visualizers.register_function(
    function=automated_rarefaction_depth.rarefy,
    inputs={'counts': FeatureTable[Frequency]}, #might be something different: technically 1 row of the feature table df
    parameters={'depth': Int % Range(1, None),
                'iteration': Int % Range(1, None),
                'seed': Int % Range(1, None)},
    input_descriptions={
        'counts': ('The feature table containing the sample to rarefy.')
    },
    parameter_descriptions={
        'depth': 'The depth to rarefy to.',
        'iteration': 'The number of iteration, how many times this specific sample has been rarified before.',
        'seed': 'The seed used for random number generation.'
    },
    name='Rarefy',
    description=("Rarefy a feature table to a specified depth."),
)


plugin.visualizers.register_function(
    function=automated_rarefaction_depth.rarefaction_depth,
    inputs={'table': FeatureTable[Frequency]},
    parameters={'seed': Int % Range(1, None),
                'percent_samples': Float % Range(0, 1),
                'iterations': Int % Range(1, 100),
                'table_size': Int % Range(1, None),
                'steps': Int % Range(5, 100),
                'algorithm': ["kneedle", "gradient"]},
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
    citations=[citations['Caporaso-Bolyen-2024']],
    #examples={'alpha_group_significance_faith_pd':
    #          ex.alpha_group_significance_faith_pd}
)

