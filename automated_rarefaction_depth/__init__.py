# flake8: noqa
# ----------------------------------------------------------------------------
# Copyright (c) 2024, Andrea Hess.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from ._version import get_versions
from ._visualizers import rarefaction_depth
from ._pipeline import pipeline_test_new, _rf_visualizer

__all__ = ['rarefaction_depth', 'pipeline_test_new', '_rf_visualizer']	

__version__ = get_versions()["version"]
del get_versions


