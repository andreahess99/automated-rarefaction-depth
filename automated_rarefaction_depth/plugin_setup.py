# ----------------------------------------------------------------------------
# Copyright (c) 2024, Andrea Hess.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from qiime2.plugin import (Plugin, Str, Properties, Choices, Int, Bool, Range,
                           Float, Set, Visualization, Metadata, MetadataColumn,
                           Categorical, Numeric, Citations, Threads)
from q2_types.feature_table import FeatureTable, Frequency
import q2_diversity
from automated_rarefaction_depth import __version__
from automated_rarefaction_depth._methods import duplicate_table

citations = Citations.load("citations.bib", package="automated_rarefaction_depth")

plugin = Plugin(
    name="rarefaction-depth",
    version=__version__,
    website="https://example.com",
    package="automated_rarefaction_depth",
    description="This qiime2 plugin gives you a tool to automatically calculate the ideal rarefaction depth based on the given data and some user parameters.",
    short_description="A tool that automatically calculates the ideal rarefaction depth.",
    # The plugin-level citation of 'Caporaso-Bolyen-2024' is provided as
    # an example. You can replace this with citations to other references
    # in citations.bib.
    citations=[citations['Caporaso-Bolyen-2024']]
)

#example method
plugin.methods.register_function(
    function=duplicate_table,
    inputs={'table': FeatureTable[Frequency]},
    parameters={},
    outputs=[('new_table', FeatureTable[Frequency])],
    input_descriptions={'table': 'The feature table to be duplicated.'},
    parameter_descriptions={},
    output_descriptions={'new_table': 'The duplicated feature table.'},
    name='Duplicate table',
    description=("Create a copy of a feature table with a new uuid. "
                 "This is for demonstration purposes only. 🧐"),
    citations=[]
)

#example pipeline
plugin.pipelines.register_function(
    function=q2_diversity.alpha,
    inputs={'table':
            FeatureTable[Frequency | RelativeFrequency | PresenceAbsence]},
    parameters={'metric': Str % Choices(alpha.METRICS['NONPHYLO']['IMPL'] |
                                        alpha.METRICS['NONPHYLO']['UNIMPL'])},
    outputs=[('alpha_diversity', SampleData[AlphaDiversity])],
    input_descriptions={
        'table': ('The feature table containing the samples for which alpha '
                  'diversity should be computed.')
    },
    parameter_descriptions={
        'metric': 'The alpha diversity metric to be computed. Information '
        'about specific metrics is available at '
        'https://data.qiime2.org/a_diversity_metrics'
    },
    output_descriptions={
        'alpha_diversity': 'Vector containing per-sample alpha diversities.'
    },
    name='Alpha diversity',
    description=('Computes a user-specified alpha diversity metric for all '
                 'samples in a feature table.')
)

#example visualizer
plugin.visualizers.register_function(
    function=q2_diversity.alpha_group_significance,
    inputs={'alpha_diversity': SampleData[AlphaDiversity]},
    parameters={'metadata': Metadata},
    input_descriptions={
        'alpha_diversity': 'Vector of alpha diversity values by sample.'
    },
    parameter_descriptions={
        'metadata': 'The sample metadata.'
    },
    name='Alpha diversity comparisons',
    description=("Visually and statistically compare groups of alpha diversity"
                 " values."),
    citations=[citations['kruskal1952use']],
    examples={'alpha_group_significance_faith_pd':
              ex.alpha_group_significance_faith_pd}
)

