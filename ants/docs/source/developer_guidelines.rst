.. _developer_guidelines:

#################################
ANTS library developer guidelines
#################################

If developing an application using ANTS and not developing the ANTS toolkit (library) itself, then the reader is referred to :ref:`writing-applications`.

Here on this page, the 'developer' refers to the person who is intending to make and submit a change to ANTS.

ANTS has a :ref:`Contributor Licence Agreement (CLA) <contribute>`.  Please ensure to sign this as a pre-request to making contributions to ANTS.

***************
Getting started
***************
Here are some important considerations to take before getting to coding:

- Will the change require some science review?
    - If your development introduces new science or changes science then you will need science review.  This will take the form of a study, a verification workflow or visual inspection.

- Does the proposed capability best fit within the ANTS core library?
    - The ANTS library represents core capability that is not particular to any one application.  Instead, ANTS represents the capability common to generating ancillaries (saving F03 ancillaries, regridding, merging etc.).  For application specific logic, this should live with the application itself (the contrib project).  The user is referred to :ref:`writing-applications` in this case.
    - Is the capability best captured in ANTS, or perhaps it is better suited to a wider audience through other libraries like Iris/Cartopy? (that is, where the capability is not limited in applicability to the domain of ancillaries).

- What existing capability is present in the library?
    - Is your usecase covered by this capability?
    - Can existing library capability be extended to meet your requirements?

- Generally small compartmentalised changes are easier to develop, manage and review.
    - Reviewing 5 small development changes stand a much better chance of getting prioritised by the reviewer than one giant behemoth of a ticket.
    - Larger developments are more prone to scope creep, hiding changes and being difficult to understand.
    - If something doesn't have to be lumped into the same pot, consider splitting this off to another ticket (branch).

A member of the :ref:`ants team <Contact>` can assist in answering these questions.  We strongly encourage people to get in contact with us prior to starting development.

Creating a ticket
*****************
Assuming that ANTS library development is determined to be the suitable pathway to follow, development begins with the creation of an ANTS ticket (:anciltrac:`click here to create ticket <newticket?component=ants&milestone=Ants%20Support&cc=miao@metoffice.gov.uk>`).

Creating a branch
*****************

ANTS developments use fcm (see :fcm:`fcm
code management <user_guide/code_management.html>` for more details).

The `ANCIL repository <https://code.metoffice.gov.uk/trac/ancil/browser>`_
contains four directories at its top level.

ants
  is the core functionality - the library and the bare minimum applications
  (e.g. regrid). This is mostly independent of particular sources or projects
  and is generally sufficiently isolated from scientists that only a technical
  review of changes is required.  
contrib
  contains end-user applications and processing code, often particular to a
  specific project or data source. Some of these are developed by the ANTS
  team and some by scientists. In general, changes here require verification
  by scientists as well as a technical review.
data
  is used to populate :mod:`$UMDIR/ancil/data`. Contains some data used for
  running models or generating ancillaries.
main
  is the CAP. This is out of scope for the ANTS team.

For a ticket involving ANTS code changes, there may be an `ants` branch, a
`contrib` branch, or both.

.. _running-tests:

Running the tests
*****************

There are end-to-end application tests in `ants` and `contrib`
which are run as part of a rose stem test workflow. Under cylc 7 this can be
launched as::

    $ cd <branch>/rose-stem
    $ rose stem --new --group=all

By default the contrib tests are run using a fixed revision of ANTS trunk.
To use a different revision or branch of ANTS, change the ``source=`` item
under ``[file:$ROSE_SUITE_DIR/share/fcm_make_ants/build]`` entry in the
``install_cold`` application configuration file in the test workflow.

When selecting the "group" to run tests for it may speed up turnaround if you
use one of the smaller, more targetted, groups for your development testing -
inspect the ``suite.rc`` file in the rose-stem suite for details of these.
However, you must ensure that full testing with the ``--group=all`` continues
to pass with your changes and that this has been run and is passing prior
to submission for review.

The application tests do include running the unittests.  For convenience, the
unittests can also be run independently of the application tests::

    $ cd <branch>
    $ pytest

The structure of the unittests reflects the hierarchy and layout of the code.
This aids transparency, discoverability and maintenance.  For example:
`ants.<module>.<Class>.<method>` has an associated unit test at
`ants/tests/<module>/<class>/test_<method>.py`.  Potentially we might also
have an integration test to ensure this class serves its purpose as a
component of a larger chain at `ants/tests/<module>/test_integration.py`.  For
a real example, take the :func:`~ants.utils.ndarray.merge_array` function
under the :mod:`ants.utils.ndarray` module.  This has an associated unit test:
`ants/tests/utils/ndarray/test_merge_array.py`.

The application tests also include running code quality assurance tests using
`black`, `flake8`, and `isort`.  Each of these tests can be run
independently from the rose suite via::

    $ cd <branch>
    $ black .
    $ flake8 .
    $ isort .

Running `black` and `isort` in your working copy fixes issues in place, so will
save time and effort over manually repairing the issues identified by `black`
and `isort` in the rose tests.

.. _buildingdocs:

Building the documentation
**************************

To build the documentation, run the :ref:`rose test workflow <running-tests>`.
This builds the documentation for applications in the core ANTS library.  The
path for the resulting documentation build can be found in the output of the
`build_docs` task.

************************
How are changes reviewed
************************

Review will require technical review and also science review where applicable.
The former is typically conducted by a member of the ANTS team.
Here we specifically address how to submit a change to ANTS.

- Before assigning a ticket for review, ensure that the full rose stem workflow has
  :ref:`been run <running-tests>`.  This will likely reduce the number of review
  iterations required.
- If your change needs science review then get this before passing for
  technical review.  Log the evidence from the science review on the ticket.
- Initiate a review by assigning the ticket for review with a member of the
  ANTS team (contact miao@metoffice.gov.uk).
- If the branch has been developed for some time, it is possible that it has
  become stale.  By this, it is meant that there is a conflict when one would
  try merging this branch to trunk.  You may be asked to resolve this using
  standard fcm workflow before going for review.  The ANTS team can be
  contacted if guidance is required.

- Documentation should be up-to-date to reflect the change
    - At minimum, docstrings should be present for all public functions,
      classes and methods.  Potentially, model documentation should also be
      provided where useful.  Also, a corresponding restructured text file may
      be required depending on the extent of the change.  See
      :ref:`here <buildingdocs>` for building the documentation.

.. _ants_release:

***********************
ANTS release management
***********************

ANTS versioning
***************

The ants version identified is available via *ants.__version__*.  This has
the following specification:

```
major.minor[.micro]
````

The "major" identifier change can result from a number of causes. It can be
used to signify a state of development or maturity of the project. It can
also be used due to a breaking change made to a public API or due to updating a
major dependency.

The "minor" identifier is used for the majority of ANTS releases. Changes
included can include backwards-incompatible changes. These will have been
deprecated according to the ANTS deprecation policy below.

The "micro" identifier is only occasionally used and represents a release result
only from a necessary bugfix applied to the previous release. There will be no
changes that require the end user to do anything differently at their end.

Trunk version identifier is simply denoted by major.minor**dev**

Deprecation Policy
******************

Typical depecation policy, where practical, is that if a breaking or
backwards-incompatible change is made to a public API, then the old behaviour
must first be deprecated. This includes adding a warning to the user if the old
behaviour is used. This should be a call to `warnings.warn` using the
`FutureWarning` category class. The warning message should include the release
number that the behaviour is deprecated from and guidance to the user regarding
how to avoid using the deprecated behaviour. Any docstring for the deprecated
behaviour should also be updated to make the deprecation clear. Deprecated
behaviour will be supported for the next two minor releases, or until the next
major release, whichever comes sooner.
