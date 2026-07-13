# flake8: noqa
# ----------------------------------------------------------------------------
# Copyright (c) 2024, Andrea Hess.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from ._version import get_versions
from ._boots_pipeline import pipeline_boots, _combined_viz


__all__ = [ 'pipeline_boots', '_combined_viz']		

__version__ = get_versions()["version"]
del get_versions


