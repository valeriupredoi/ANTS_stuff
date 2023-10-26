========================================
Appendix A: Ancillary file time metadata
========================================

This
appendix outlines how ANTS deals with the time metadata when writing `F03`_ ancillary files. You will need to understand this when preparing data for processing by ANTS.  It will help you understand why ANTS will return an error when processing some files.  It will also help you understand the content of F03 ancillary files.

ANTS aims to write ancillary files that:

  1. Are usable by the Unified Model.
  2. Help all users of the ancillary files understand how appropriate an ancillary might be for their model run.
  3. Keep ancillary files, as far as possible, consistent with other files described by F03.

Not all applications that read F03 ancillaries use the metadata in the same way.   This means the ancillary files produced by ANTS will not necessarily be read cleanly by other applications.

.. _`F03`: https://code.metoffice.gov.uk/doc/um/latest/papers/umdp_F03.pdf

Time based metadata in an ancillary file
----------------------------------------

F03 ancillaries record time related metadata in both the fixed length header and the lookup header. In ANTS we use the lookup header to keep additional time metadata that is useful to all users of the data, not just the UM.  Depending on the configuration, the UM may only use the time information in the fixed length header to determine when to update the ancillary field.

When you are preparing input data for ANTS you will concentrate on getting the metadata that goes into the lookup header correct.  ANTS infers the elements of the fixed length header used by the UM from the lookup header information.  If you are preparing NetCDF files for processing by ANTS you will need to get the CF metadata correct.  If you are preparing pp [1]_ data you will need to get the pp header correct.

.. [1] We recommend converting to pp source files to NetCDF where possible first.  We recognise, however, that there is still a large use of pp. Reading of pp source data is supported by ANTS.

The UM recognises three types of time updating in ancillary fields.  These are:

   1. fixed fields, these are never updated through a run.  Fixed fields can be used for runs starting any year and with any calendar.
   2. time series fields, these are updated through a run.  Time series fields are for a set of dates and they should only be used for runs that cover those dates.  Some time series fields can be used for runs with any calendar, while others can only be used for runs with a particular calendar (more on this later). 
   3. periodic fields, these are updated during a run and the fields are thought of as an infinite repeating time series. In practice periodic fields can be used for runs starting in any year.  In the current ANTS implementation periodic fields can be used with any calendar.


Fixed, periodic and time series fields
--------------------------------------

ANTS infers whether data is fixed, a time series based or periodic based on the time data and metadata associated with a field.  ANTS will set the value of the time indicator (item 10) in the fixed length header to 0 for fixed fields, 1 for time series fields and 2 for periodic fields.

Any fields with a single time point are written as fixed fields.  All data that is neither a single time or inferred as periodic is written as a time series field.  

The algorithm used to infer that data is periodic is work in progress.  Currently, an ancillary field is inferred to be periodic if it is an annual time series of monthly mean data.  This simple algorithm can be changed as we see more real examples of periodic fields and as we implement support for multi year climatologies.  In a `CF`_ NetCDF file the field should have ``cell_methods= "time: mean"`` and the time coordinate must have bounds.  In pp and fields files the ``IB`` component of ``LBTIM`` must be set to 2 and ``LBPROC`` set to 128 to indicate a time mean and the time entries (LBYR, LBMON, LBDAT) should be set to the start and end of each month.

.. _`CF`: http://cfconventions.org/latest.html

Many periodic ancillary fields are derived from multi year climatologies.  These are currently not supported in ANTS. Until we have added support for multi year climatologies set a field's time coordinate and  metadata as if it was a single year mean.  Choose a year that is representative of the climatology.  Where possible produce a file to record the climatology period to live alongside the source data.  Record this time metadata in an 'ini' style file.  For instance if you have monthly mean data for the climatology period January 1989 to December 2000 your metadata file would be::

  [climatology]
  cell_methods: time: mean within years time: mean over years
  bounds_diff: month
  start: 1989-01-01
  end: 2001-01-01

Include the `CF`_ ``cell_methods`` appropriate for your data. The start and end are the first and last bounds of the `CF`_ climatology bounds.  This is a temporary solution until we support climatologies.

Note that, while fixed fields and periodic fields can be used for runs starting in any year, you will need to understand whether the time period of the source data is appropriate for your run.

Calendar information
--------------------

When preparing input data in pp or NetCDF files for ANTS you should always use the calendar that describes the data.  So for real world data this is most likely to be the Gregorian or Proleptic Gregorian calendars, but for data derived from other model runs this could be Proleptic Gregorian, 360 day or 365 day calendars.  ANTS will write this calendar to the lookup header of the ancillary file.  This gives other users of the ancillary a way of understanding the source of the data and help them make informed decisions about whether a particular ancillary can be used for their model run.

The calendar written to the fixed length header (fixed length header item 8) will not necessarily be the same as that in the lookup header.  For fixed fields and periodic or time series data with intervals greater than a month the fixed length header calendar will be written as IMDI.  This is because in many cases these ancillary files can be used in UM runs that use any calendar.  This avoids having to generate ancillary files that differ only in the calendar indicator of the fixed length header.  

For time series data with an interval less than a month ANTS will write the actual calendar to the fixed length header.

Date information
----------------

Just as for the calendar, when preparing data in NetCDF or pp files for ANTS you should record the dates that best describe the data.  If you are preparing periodic data that is based on an observed or model climatology then use the dates of that climatology.  ANTS will write this information to the lookup header of the ancillary file.  This then helps users of the data understand whether an ancillary is appropriate for their model run.  Note this differs from some legacy ancillaries where year 0 was used in the lookup header.  Writing the year of the data  helps other data users understand the applicability of the data to their model run.

If your data is a time mean between two dates then the bounds should be written in a way that is consistent with both UM fields files that contain diagnostic output and with the CF convention.  This means the end bounds of one interval should be identical to the start bounds of the next interval.  For example if you have monthly mean data in a pp file then the time header elements should be written as::

   Field 1 LBYR 2000
           LBMON 1
           LBDAT 1
           LBHR 0
           LBMIN 0
           
           LBYRD 2000
           LBMOND 2
           LBDATD 1
           LBHRD 0
           LBMIND 0

   Field 2 LBYR 2000
           LBMON 2
           LBDAT 1
           LBHR 0
           LBMIN 0
           
           LBYRD 2000
           LBMOND 3
           LBDATD 1
           LBHRD 0
           LBMIND 0

Note this differs from some legacy ancillaries which would use the last minute of the last hour as the end time of an interval.  Making the bounds contiguous as advised here is both consistent with other uses of time bounds and consistent with the documented specifications.

ANTS will use the time information in your file to infer the values of the start date, end date and interval to put into the fixed length header.  All fixed fields will have a start date, end date and interval of 0 (all items 21-41 will be set as 0); in practice these are ignored by the UM.

Final Note
----------

We recognise that the use of time information by ANTS is different, in places, to legacy ancillary files.  This will impact you as you may need to change some practices and deal with sets of ancillary file sets that are mixes of the ANTS use of metadata and some legacy metadata.  Where we have changed the use of metadata we have tried to do it to:

  1. make it easier to automatically infer other metadata elements, so you don't have to put that metadata in your self
  2. make is easier for other users of the data to understand the range of applicability of an ancillary file.

We hope these benefits are worth any inconvenience in the transition period.

We have also made choices for interpreting time metadata to make F03 ancillaries based on the sample data we have seen.  These will not cover all cases, so please :ref:`contact us<about>` if the interpretation we currently use does not match your expectations.


