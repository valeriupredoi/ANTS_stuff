.. _writing-applications:



############################
Developing a new application
############################
User applications for generating ancillaries are held under the :contrib:`contrib<>` project.
These applications typically utilise the ants library to aid in the generation of ancillaries.

As with :ref:`developing ANTS itself <developer_guidelines>`, it is advised
that users :anciltrac:`create a ticket <newticket?component=contrib&milestone=Ants%20Support&cc=miao@metoffice.gov.uk>`
to log their progress and increase visibility of their development/intentions.
We also advise that users get in contact with the relevant science owner for
the ancillary (where applicable) as well as getting in contact with
:ref:`us<about>`.

Applications for generating ancillaries are typically split into three
categories: pre-processing, ancillary generating and post-processors.
The reader is referred to :ref:`ants_principles` for further details of
these categories.  However, in the discussion of developing a new ancillary
application, the reader should carefully consider good workflow as highlighted
in this guide to ensure that the generation of an ancillary is both flexible
and maintainable and that the processing made is defined under the correct
application type (pre-processor vs ancillary generator).
To summarise, pre-processors are model independent processing while ancillary
generating applications are model dependent processing.  Generally speaking,
post-processors are not desirable and should be avoided where possible.

The ANTS :ref:`toolkit` (aka library) provides functionality common to many
ancillary applications.  When writing a new application try to use such
functionality where possible.  If you need guidance then please
:ref:`contact us<about>`.

.. _sample-applications:

********************
A sample application
********************

Application development should target the :contrib:`contrib<>`
project rather than ANTS itself.

An application will typically do the following:

  1. Read input, such as input and output file names, target
     grid definition.
  2. Load the data from the input files.
  3. Process (e.g regrid) the source data to produce the ancillary field(s).
     This may need to decompose the analysis into sub-domains to help with
     performance and/or reduce the memory footprint.
  4. Save the data to the output file.

ANTS has code to help with all these steps.  In particular,
:mod:`ants.command_parse` provides a common commandline interface, which is
important to use to ensure a consistent experience (UI) between all
applications.  It also provides a number of benefits which are highlighted in
the module documentation.  ANTS also utilises a run-time configuration which
is detailed in :mod:`ants.config`.  As ANTS is based on iris, you can also
use any of the iris functionality to process your data.

Because contrib is not within the ANTS unit testing framwork, an application
should contain minimal code. In general, if the code is complex enough to
require unit tests, it should probably be in the ANTS library.

The following represents a simple example application, where the reader should
take particular attention to application documentation.  The application
documentation should follow `numpydoc style
<https://numpydoc.readthedocs.io/en/latest/example.html>`_::

    """
    Application name
    ****************

    - Description of both input datasets as well as output datasets:
        - Wat does the quantity represent?
        - What is the identifying information (stash, standard_name, ...)
    - If this represents an ancillary generating application, then which
      schemes use this ancillary?
    - Bullet point outline of processing done.  This should detail all steps
      taken by the application without going into technical details.  Formula
      should also be provided where applicable, such as:

    .. math::

        a**2 = b**2 + c**2

    """
    import ants
    import ants.io.save as save
    import ants.decomposition as decomp


    def process(source, target):
        # The core processing needed to go from the source data
        # to the target grid.
        mean_cube = ants.analysis.mean(source, target)
        return mean_cube


    def load_data(source_filepath, target_filepath):
        # Load the data from the source file and the grid from the target file.
        # This example assumes only one input variable in the source file.
        source_cube = ants.load_cube(source_filepath)
        target_cube = ants.load_cube(target_filepath)
        return source_cube, target_cube


    def main(source_filepath, target_filepath, output_filepath):
        # This is the top-level function for an application
        # The filepath arguments should be read from the command line
        # before calling this routine.
        source_cube, target_cube = load_data(source_filepath, target_filepath)

        # This next line calls the decomposition framework.
        # The decomposition framework divides the target and source into
        # sub domains and then calls ``process`` on each sub-domain.
        cube_result = decomp.decompose(process, source_cube, target_cube)
        filler = ants.analysis.FillMissingPoints(cube_result, target_cube)
        filler(cube_result)

        # Save the data to a netcdf and UM ancillary file.
        save.netcdf(cube_result, output_filepath)
        save.ancil(cube_result, output_filepath)

        return cube_result


    if __name__ == '__main__':
        # Keep this __main__ small: just parse the command line
        # then call the main application function.
        parser = ants.AntsArgParser(target_lsm=True)
        args = parser.parse_args()
        main(args.sources, args.target_lsm, args.output)


It is convention that applications should contain both a `main` function
and a `load_data` function. All data required in the application should
be loaded in the `load_data` function.

**************
Pre-processors
**************

In the general case, pre-processing can be split into two categories: fixing of
source datasets and model independent processing.
Ancillary generation applications utilise the output of pre-processing
applications.  Given that it is model independent, they persist on disk and
need only be run once.  NetCDF is the format of choice for the output of
pre-processed sources as it allows data to include additional metadata as well
as well as provide fast and efficient read access.

Common fixes to source datasets include: adding/fixing coordinate systems and
associated metadata (see helper function :func:`ants.utils.cube.set_crs`);
fixing time coordinates for better defining climatologies (see helper function
:func:`ants.utils.cube.set_month_mean_for_year`);  adding a 'um_runid'
attribute (where applicable); adding 'source' attribute to point to the origin
of the dataset and also removing any metadata not applicable to the dataset.

A good example split in this pre-processing chain, is the soils generation.
Here we see two pre-processors :contrib:`ancil_soils_preproc_cosby.py<DataPreparation/soils_cosby/ancil_soils_preproc_cosby.py>`
and :contrib:`ancil_soils_preproc.py<DataPreparation/soils/ancil_soils_preproc.py>`.
The former application fixes the source dataset while the later represents
the model independent part of the processing which is involved in deriving the
soil parameters.

While the majority of applications are located in the
:contrib:`Apps directory<Apps/>`, some preprocessors are stored in the
:contrib:`DataPreparation directory<DataPreparation/>` dirctory. This directory
is for preprocessors which are kept as a record of what was done, rather than
scripts which will need to be run again - for example, scripts used to
preprocess source data into UMDIR master files. For this case, please see the
information regarding
:ancilwiki:`updating files in UMDIR <ANTS/ProjectManagement/updating_UMDIR>`.

********************
Filename conventions
********************

To aid transparency and readability, we adopt a filename convention for
applications:

`ancil_<class>[_preproc[_<source>]].py`

Where 'class' is the class of ancillary we are generating like
'topographic_index' or 'sst_seaice'.  'preproc' denotes that it's a
pre-processor and 'source' denotes the source dataset being pre-processed if
it's a pre-processor for fixing a source dataset.

Taking our example from the previous section, our pre-processor of the cosby
source dataset is named :contrib:`ancil_soils_preproc_cosby.py<DataPreparation/soils_cosby/ancil_soils_preproc_cosby.py>`.  Our more general
model independent processing is named :contrib:`ancil_soils_preproc.py<DataPreparation/soils/ancil_soils_preproc.py>`.  Our ancillary
generating application is then named :contrib:`ancil_soils.py<Apps/SoilParameters/ancil_soils.py>`.

.. _testing-applications:

***********************
Testing the application
***********************

It is suggested that each application has an end-to-end test associated with
it.  This ensures that applications actually run as expected.
These end-to-end tests take the form of a cylc workflow. Under cylc 7 this can
be lauched as::

    $ cd rose-stem
    $ rose stem --new --group=all

This workflow is located :contrib:`here<rose-stem>`.
Each rose application test uses
:contrib:`sample datasets<rose-stem/sources>`.  These should be
very low resolution versions of the original source datasets, see
:ref:`here<gen_data>`. Tests can be added to existing groups as required but
should always be added such that they will be run when using the all group.

Ancillaries are generated by the test workflow and compared against :contrib:`known good
output ancillaries<KGO>`.

Tests include:

    1. Metadata comparisons for both NetCDF and UM ancillary fields files.
    2. Unittests where applicable.

Where the application is non-trivial, it may be advised that additional testing
is warranted.  In which case, the user is referred to contrib application
:contrib:`Lai<Apps/Lai>` for an example of how this can be done.  Note that such tests
are again run within the :contrib:`rose-stem<rose-stem>` workflow.

If the application uses :ref:`decomposition<decomposition>`, the end-to-end tests
should be run with both decomposition disabled and enabled and separate known
good outputs for each case provided. Please refer to the
:contrib:`rose-stem<rose-stem>` workflow for existing examples.

Feel free to :ref:`contact us<about>` if you need some assistance/guidance.

.. _gen_data:

Hints and tips for test data production
---------------------------------------

When you are producing sample test input then the following guidelines can
help you avoid some common pitfalls.

* Keep the test data small in size (less than 600 Kilobytes where possible,
  roughly 1x1 deg resolution for a single level field).  Keeping the test data
  small means it can be kept in the repository, and the tests can run
  relatively quickly.
* Ensure the test data is representative of the source.  This can be done by
  sub-sampling the original source data and randomising the data in some way
  to avoid any issues relating to License terms (if applicable).  We should
  ensure our sample data has the same characteristics:

  * Contiguous/discontinuous grid cells.
  * Any uneven (stretched) coordinates.
  * Global or regional.
  * Coordinate reference system or lack of.
  * Has any defects contained in the original source.

* Check the license conditions on the source. If the license means the data
  can not be put in a public repository then you will need to derive an
  alternative source.


.. _contrib_release:

**************************
Contrib release management
**************************

Contrib release management runs parallel to ants.  This means that those
packages included in a tag release of contrib vx.x are designed to work
against the corresponding tag release of ants (vx.x).

Not all content within contrib will be included within a tag release.
There are a few preconditions.  Applications will:

- Have a corresponding end-to-end rose application acceptance test (see
  :ref:`testing-applications`).
- Be tested and validated, against the version of ants that contrib targets.

  - This means that contrib applications need to be updated and validated for
    each contrib release or otherwise not be included in that release.
  - Applications can be retrospectively added to a contrib tag release, so
    long as the change is additive.
  - Those applications which are one-off usage can be removed as they remain
    available at the revision they were used.

- Include module documentation (see :ref:`sample-applications`).
