# automated-rarefaction-depth

A [QIIME 2](https://qiime2.org) plugin [developed](https://develop.qiime2.org) by Andrea Hess (andhess@ethz.ch). 🔌

## Installation instructions

### Install Prerequisites

[Miniconda](https://conda.io/miniconda.html) provides the `conda` environment and package manager, and is currently the only supported way to install QIIME 2.
Follow the instructions for downloading and installing Miniconda.

After installing Miniconda and opening a new terminal, make sure you're running the latest version of `conda`:

```bash
conda update conda
```

###  Install development version of `automated-rarefaction-depth`

Next, you need to get into the top-level `automated-rarefaction-depth` directory.
If you already have this (e.g., because you just created the plugin), this may be as simple as running `cd automated-rarefaction-depth`.
If not, you'll need the `automated-rarefaction-depth` directory on your computer.
Download this repository.
Once you have the directory on your computer, change (`cd`) into it.

If you're in a conda environment, deactivate it by running `conda deactivate`.


Then, run:

```shell
conda env create -n automated-rarefaction-depth-dev --file ./environments/automated-rarefaction-depth-qiime2-amplicon-2024.5.yml
```

After this completes, activate the new environment you created by running:

```shell
conda activate automated-rarefaction-depth-dev
```

Finally, run:

```shell
make install
```

```Install the Kneed library
pip install kneed
```
Start by making QIIME 2's command line interface aware of `automated-rarefaction-depth` by running:

```shell
qiime dev refresh-cache
```

You should then see the plugin in the list of available plugins if you run:

```shell
qiime info
```

You should be able to review the help text by running:

```shell
qiime rarefaction-depth --help
```

Have fun! 😎

## About

The `automated-rarefaction-depth` Python package was [created from template](https://develop.qiime2.org/en/latest/plugins/tutorials/create-from-template.html).
To learn more about `automated-rarefaction-depth`, refer to the [project website](https://example.com).
To learn how to use QIIME 2, refer to the [QIIME 2 User Documentation](https://docs.qiime2.org).
To learn QIIME 2 plugin development, refer to [*Developing with QIIME 2*](https://develop.qiime2.org).

`automated-rarefaction-depth` is a QIIME 2 community plugin, meaning that it is not necessarily developed and maintained by the developers of QIIME 2.
Please be aware that because community plugins are developed by the QIIME 2 developer community, and not necessarily the QIIME 2 developers themselves, some may not be actively maintained or compatible with current release versions of the QIIME 2 distributions.
More information on development and support for community plugins can be found [here](https://library.qiime2.org).
If you need help with a community plugin, first refer to the [project website](https://example.com).
If that page doesn't provide information on how to get help, or you need additional help, head to the [Community Plugins category](https://forum.qiime2.org/c/community-contributions/community-plugins/14) on the QIIME 2 Forum where the QIIME 2 developers will do their best to help you.
