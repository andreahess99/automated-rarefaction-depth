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
    function=automated_rarefaction_depth.rarefaction_depth,
    inputs={'table': FeatureTable[Frequency], 'phylogeny': Phylogeny[Rooted]},
    parameters={'metadata': Metadata,
                'p_samples': Float % Range(0, 1),
                'iterations': Int % Range(1, None)},
    input_descriptions={
        'table': ('Feature table to compute rarefaction curves from.'),
        'phylogeny': ('Optional phylogeny for phylogenetic metrics.')
    },
    parameter_descriptions={
        'metadata': 'The sample metadata.',
        'p_samples': 'The lower bound of the percentage of your samples want to keep.',
        'iterations': 'The number of rarefied feature tables to '
                       'compute at each step.',
                       'metadata': 'The sample metadata.'
    },
    name='Automated Rarefaction Depth',
    description=("Automatically computes an optimal rarefaction depth."),
    citations=[citations['Caporaso-Bolyen-2024']],
    #examples={'alpha_group_significance_faith_pd':
    #          ex.alpha_group_significance_faith_pd}
)

