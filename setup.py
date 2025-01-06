# ----------------------------------------------------------------------------
# Copyright (c) 2024, Andrea Hess.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from setuptools import find_packages, setup

import versioneer

description = ("A template QIIME 2 plugin.")

setup(
    name="automated-rarefaction-depth",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    license="BSD-3-Clause",
    packages=find_packages(),
    author="Andrea Hess",
    author_email="andhess@ethz.ch",
    description=description,
    url="https://example.com",
    entry_points={
        "qiime2.plugins": [
            "automated_rarefaction_depth="
            "automated_rarefaction_depth"
            ".plugin_setup:plugin"]
    },
    package_data={
        "automated_rarefaction_depth": ["citations.bib", "assets/*", "assets/js/*", "assets/css/*"],
        "automated_rarefaction_depth.tests": ["data/*"],
    },
    zip_safe=False,
)
