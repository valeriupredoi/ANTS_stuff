.. _running-applications:

####################
Running applications
####################

Applications are split between :ref:`core-applications` and :ref:`contrib-applications`.
The former is owned by the ANTS team and function as general purpose scripts applicable to many people/workflows.
The later concerns science owned applications which typically utilise ants to generate ancillaries.

*****************
Commandline usage
*****************

Applications are run utilising a commandline interface which provides a common UI for generating ancillaries.
The reader is referred to :class:`ants.command_parse.AntsArgParser` for documentation of this common interface.
For example usage, see :ref:`current-applications`.

Here is an example where we call an application from the commandline::

    $ ancil_general_regrid.py <source_filepath> --target-grid <target_grid_filepath> \
    -o <output_filepath> --ants-config <config_filepath>

The purpose of this call will be to utilise the `ancil_general_regrid.py`
application to regrid our specified source to the specified target grid.
`config_filepath` is a filepath to a configuration file (more details can be found in
`runtime_configuration`_).
This configuration file can be used to customise run-time characteristics of ANTS,
for example the decomposition configuration or even the regridding methods utilised
(see :mod:`ants.config` for further details).

This is an example of a configuration file which could be passed to our application above::

    [ants_decomposition]
        x_split = 10
        y_split = 10

    [ants_regridding_horizontal]
        scheme = Linear

Let's take a closer look:

Here we split our source-target pair into 10x10 pieces (100 in total)::

    x_split = 10
    y_split = 10

Lastly, we define our horizontal regridding scheme utilised by ANTS::

    scheme = Linear

**********************
Rose application usage
**********************

Here we show how we include our configuration within a `rose-app.conf` file::

    [command]
        default=ancil_generalised_regrid.py ${source} --target-grid ${target} --ants-config ${ANTS_CONFIG}

    [env]
        ANTS_CONFIG=${PWD}/rose-app-run.conf
        output=<output filepath>
        source=<source filepath>

    [decomposition]
        x_split = 10
        y_split = 10

    [ants_regridding_horizontal]
        scheme = Linear

    [ants_regridding_vertical]
        scheme = Linear
        extrapolation_mode = Nearest

Note, that here ``--ants-config`` provides a path to the rose application's own `rose-app.conf` file.
ANTS respects schedulers (SLURM, PBS and LSF) which means that it is the sheduler which informs ants of how many
processes it should use (and so whether a serial or parallel job).  However, users still have the
ability to specify the number of processes to utilise in the case where no scheduler is used.
Additionally, users have many other options to how they confiure ants run-time (inc. the regridder to use).
See :mod:`ants.config` for more details.

A good real working example of how to get up and running within rose and ants is via the :contrib:`rose-stem`.
Points of particular note include ensuring that dask plays nicely with numpy (single threaded)
(as mentioned in the dask documentation ref:`here <https://docs.dask.org/en/stable/array-best-practices.html#avoid-oversubscribing-threads>`_)::

    export OMP_NUM_THREADS=1
    export OPENBLAS_NUM_THREADS=1
    export MKL_NUM_THREADS=1
    export VECLIB_MAXIMUM_THREADS=1
    export NUMEXPR_NUM_THREADS=1

This can be seen defined in the environment setup wrapper for the :source:`ants-launch script <bin/ants-launch>`. For more details on ants-launch see :ref:`ants_launch`.

.. _runtime_configuration:

**************************
ANTS runtime configuration
**************************

For ANTS runtime configuration, the reader should refer to :mod:`ants.config`.
Illustrative usage can be found in the above two sections.
