# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.

"""
River routing application
*************************

Sort NEMO rivers into an order, giving each river outflow a unique ID number.
Put these river numbers into UM and NEMO river_number ancillary files.

* Load the NEMO land-sea mask from NEMO domain config file.
* Load the NEMO runoff file and assign a uniqie river_number ID to each
  defined river.
* At each TRIP outflow point perform a diamond_search to find the nearest
  NEMO outflow point. Copy the river_number from the NEMO grid to the TRIP
  grid for these matching points.
* Remove any unused NEMO river numbers, cascading up subsequent river numbers
  so that there are no gaps.
* Output NEMO and TRIP/UM river numbers as ancillary files.

"""
import logging
import math

import ants.io.save as save
import iris
import numpy as np
from netCDF4 import Dataset

from . import location_class

_LOGGER = logging.getLogger(__name__)


small_amount = 1.0e-5
# This is the threshold at which we use NEMO's rivers. Anything under this is not
# a river outflow point unless we can't find any nearby river outflow points
# when we will make one.


def get_neighbouring_points(runoff_nc_array, location_list, river_index_array):
    """
    Gets all the neighbouring grid boxes and fills in river_index_array with
    river numbers

    Parameters
    ----------
    runoff_nc_array : :class:`numpy.ndarray`
       Array of runoff values from NEMO's climatological runoff file.
       This array will gradually be reset to zero) by subsequent calls
       to this routine.
    location_list : list
       List of location classes to find neighbours of
    river_index_array : :class:`numpy.ndarray`
       Array of river numbers on the NEMO grid. This will be filled in
       by subsequent calls to this routine.

    Returns:
    -------
    new_location_list : list
       List of neighbours (where each list item is a location class)

    """
    new_location_list = []
    for location in location_list:

        for direction in range(1, 9):
            location.shift(direction)
            newj, newi = location.j, location.i
            location.revert()

            if not math.isnan(newj):
                if (
                    runoff_nc_array[newj, newi] <= location.amount * 1.1
                    and river_index_array[newj, newi] == 0
                    and runoff_nc_array[newj, newi] > small_amount
                ):
                    # A valid neighbouring river outflow has been found

                    # Assign river number to output array
                    river_index_array[newj, newi] = location.river_number

                    # Build up a new list of locations for this river (for the next
                    # call to this routine)
                    new_location_list.append(
                        location_class.LocationClass(
                            newj,
                            newi,
                            runoff_nc_array.shape,
                            amount=runoff_nc_array[newj, newi],
                            river_number=location.river_number,
                        )
                    )

                    # Reset river outflow to zero so this point becomes invalid and
                    # can not be found again
                    runoff_nc_array[newj, newi] = 0.0

    return new_location_list


def resort_rivers(river_index_array, accumulative_amount_list):
    """
    river_index_array currently has river numbers in order of single point amounts
    and not accumulated anounts for the whole runoff region. This routine reorders
    these river numbers to be for the whole accumulated amounts.

    Parameters
    ----------
    river_index_array : :class:`numpy.ndarray`
       Array of river numbers on the NEMO 2D spatial grid
    accumulative_amount_list : list
       List over all rivers with each item being a list of
       [river_number,accumulated amount]

    Returns
    -------
    rindex_new_array : :class:`numpy.ndarray`
       Modified array of river numbers, ordered by cumulative runoff amounts.

    """
    # Make a blank river_index_array to put the data in
    rindex_new_array = river_index_array.copy()
    rindex_new_array[:, :] = 0

    # Sort the accumulative amounts so the largest accumulated outflows appear first
    accumulative_amount_sorted = sorted(
        accumulative_amount_list,
        key=lambda accumulative_amount_list: accumulative_amount_list[1],
        reverse=True,
    )

    new_river_number = 1
    for accumulative_amount in accumulative_amount_sorted:

        # Find the old river indexes and replace them with the new ones
        where_river = river_index_array == accumulative_amount[0]
        rindex_new_array[where_river] = new_river_number
        new_river_number = new_river_number + 1

    return rindex_new_array


def remove_outflows_invalid(orca_dom_file, runoff_nc_array):
    """
    Remove NEMO outflow points that are not valid such as at land points
    or those in the wraparound or north fold points.

    Parameters
    ----------
    orca_dom_file : str
       Input file (netcdf) contining the NEMO domain on the ORCA grid
    runoff_nc_array : :class:`numpy.ndarray`
       Array of runoff values from NEMO's climatological runoff file.
       This array will gradually be reset to zero) by subsequent calls
       to this routine.
    """

    # If a mask is suppled then mask out all river outflows not in ocean points
    if orca_dom_file:
        _LOGGER.info("Removing rivers not at ocean points")
        rootgrp = Dataset(orca_dom_file, "r")
        tmask = rootgrp.variables["top_level"][0, :, :]
        rootgrp.close()

        where_land = tmask == 0
        where_rivers_land = np.logical_and(runoff_nc_array > 0.0, where_land)
        _LOGGER.info(("Number of river outflow on land = ", np.sum(where_rivers_land)))

        runoff_nc_array[where_land] = 0.0

    # Remove all river outflows from the 1 grid box halo around the whole field
    runoff_nc_array[0, :] = 0.0
    runoff_nc_array[-1, :] = 0.0
    runoff_nc_array[:, 0] = 0.0
    runoff_nc_array[:, -1] = 0.0


def assign_river_number_to_points(
    river_index_array,
    runoff_nc_array,
    twod_index,
    river_number,
    amount,
    accumulative_amount_list,
):
    """
    Assign river numbers to NEMO grid points. This is done at one central point
    and then all neighbouring points that have the same outflow amount or less
    (but not zero) are also assigned the same number.

    Parameters
    ----------
    river_index_array : :class:`numpy.ndarray`
       Array of river numbers on the NEMO 2D spatial grid
    runoff_nc_array : :class:`numpy.ndarray`
       Array of runoff values from NEMO's climatological runoff file.
    twod_index : tuple
       Indices (y, x) of river_index_array and runoff_nc_array where the
       maximum runoff is
    river_number : int
       River number to assign
    amount : float
       River outflow amount for this point
    accumulative_amount_list : list
       List over all rivers with each item being a list of
       [river_number,accumulated amount]

    """

    river_index_array[twod_index] = river_number
    runoff_nc_array[twod_index] = 0.0

    # Start a location list that will grow to include all the points for
    # this river outflow
    location_list = [
        location_class.LocationClass(
            twod_index[0],
            twod_index[1],
            runoff_nc_array.shape,
            amount=amount,
            river_number=river_number,
        )
    ]
    accumulative_amount = amount

    # Iteratively grow the number of points for this river outflow
    # adding river numbers to river_index_array as we go.
    for loop in range(1000):
        location_list = get_neighbouring_points(
            runoff_nc_array, location_list, river_index_array
        )
        if len(location_list) == 0:
            break
        for location in location_list:
            accumulative_amount = accumulative_amount + location.amount

    accumulative_amount_list.append([river_number, accumulative_amount])


def generate_nemo_river_numbers(ocean_runoff_file, orca_dom_file=None):
    """
    Generates an index for each ocean runoff region.

    Parameters
    ----------
    ocean_runoff_file : str
       Filename of netcdf file that contains the climatological
       runoff used by NEMO in standalone mode.
    orca_dom_file : str
       Filename of netcdf file containing the domain of the NEMO model.
       This is used to get the land-sea mask on the NEMO grid.

    Returns
    -------
    river_index_array : :class:`numpy.ndarray`
       Array of river numbers on the NEMO 2D spatial grid
    nav_lat : :class:`numpy.ndarray`
       Array of latitudes for the NEMO 2D spatial grid
    nav_lon : :class:`numpy.ndarray`
       Array of longitudes for the NEMO 2D spatial grid

    """
    # Set a maximum number of rivers in the NEMO file
    rmax = 4000

    # Load the NEMO river outflows
    rootgrp = Dataset(ocean_runoff_file, "r")

    # Find the variable and dimension names you need
    var_name = "sornficb"
    if "sorunoff" in list(rootgrp.variables.keys()):
        var_name = "sorunoff"

    # Get the data from the river outflow netcdf file
    runoff_nc_array = rootgrp.variables[var_name][:, :, :].mean(axis=0)
    nav_lat_ncvar = rootgrp.variables["nav_lat"]
    nav_lon_ncvar = rootgrp.variables["nav_lon"]
    nav_lat = nav_lat_ncvar[:, :]
    nav_lon = nav_lon_ncvar[:, :]
    rootgrp.close()

    # Remove NEMO outflow points that are not valid
    remove_outflows_invalid(orca_dom_file, runoff_nc_array)

    # Make a new array to hold the river indices
    river_index_array = runoff_nc_array.copy()
    river_index_array[:, :] = 0
    accumulative_amount_list = []

    # Loop over river numbers
    for river_number in range(1, rmax):

        # Find the maximum
        oned_index = runoff_nc_array.argmax()
        twod_index = np.unravel_index(oned_index, runoff_nc_array.shape)

        # If the runoff is tiny then stop here and exit out of this loop
        amount = runoff_nc_array[twod_index]
        if amount < small_amount:
            break

        # Assign river numbers to these NEMO points
        assign_river_number_to_points(
            river_index_array,
            runoff_nc_array,
            twod_index,
            river_number,
            amount,
            accumulative_amount_list,
        )

        if river_number == rmax - 1:
            error_str = (
                "ERROR: Exceeded maximum number of rivers. "
                "Increase rmax in order_nemo_rivers.py"
            )
            _LOGGER.error(error_str)
            raise Exception(error_str)

        if river_number == 44 or river_number == 55:
            _LOGGER.info(
                "River number ",
                river_number,
                " has amount ",
                amount,
                " at location ",
                twod_index,
            )

    # Resort the river numbers by their accumulative amount
    river_index_array = resort_rivers(river_index_array, accumulative_amount_list)

    return river_index_array, nav_lat, nav_lon


def blank_nemo_river_numbers(orca_dom_file):
    """
    Make a NEMO array of river indices where each index is zero

    Parameters
    ----------
    orca_dom_file : str
       Filename of netcdf file that contains the NEMO ORCA
       domain information. It must have nav_lat and nav_lon
       netcdf fields within it.

    Returns
    -------
    river_index_array : :class:`numpy.ndarray`
       Array of river numbers on the NEMO 2D spatial grid
    nav_lat : :class:`numpy.ndarray`
       Array of latitudes for the NEMO 2D spatial grid
    nav_lon : :class:`numpy.ndarray`
       Array of longitudes for the NEMO 2D spatial grid

    """
    # Load the NEMO resolution data
    rootgrp = Dataset(orca_dom_file, "r")

    # Get the data from the netcdf file
    nav_lat_ncvar = rootgrp.variables["nav_lat"]
    nav_lon_ncvar = rootgrp.variables["nav_lon"]
    nav_lat = nav_lat_ncvar[:, :]
    nav_lon = nav_lon_ncvar[:, :]
    rootgrp.close()

    # Make a new array to hold the river indices
    river_index_array = nav_lat.copy()
    river_index_array[:, :] = 0

    return river_index_array, nav_lat, nav_lon


def scale_by_resolution(nav_lat, nav_lon):
    """
    Define a scaling number that is linear with horizontal
    ocean resolution and of value 1 for ORCA025

    Parameters
    ----------
    nav_lat : :class:`numpy.ndarray`
       Array of latitudes on the NEMO 2D spatial grid
    nav_lon : :class:`numpy.ndarray`
       Array of longitudes on the NEMO 2D spatial grid

    Returns
    -------
    scaling : float
       Scaling number to apply to distances.
       1 for ORCA025 and scales linearly with horizontal resolution

    """
    nlats = nav_lat.shape[0]
    nlons = nav_lat.shape[1]

    # Calculate the grid spacing in both directions (but we will only use lons)
    lons_diff = (
        nav_lon[int(nlats / 2), 1 : nlons - 1] - nav_lon[int(nlats / 2), 0 : nlons - 2]
    )

    # Remove any stupidly large values due to wrap around
    large_neg = np.where(lons_diff < -300)
    if len(large_neg) > 0:
        lons_diff = np.delete(lons_diff, large_neg)
    large_pos = np.where(lons_diff > 300)
    if len(large_pos) > 0:
        lons_diff = np.delete(lons_diff, large_pos)

    # Generate the mean
    resolution = lons_diff.mean()

    # Generate a scaling
    scaling = 0.25 / resolution
    return scaling


def get_min_distance(j, i, sequence_cube, scaling):
    """
    Determine the distance to search for your NEMO outflow point.

    Parameters
    ----------
    j : int
       j index of TRIP data point to examine
    i : int
       i index of TRIP data point to examine
    sequence_cube : :class:`iris.cube.Cube`
       River sequence numbers as output from ancil_river_routing.river_routing
    scaling : float
       Scaling factor to apply to distances.
       1 for ORCA025 and scales linearly with horizontal resolution

    Returns
    -------
    min_distance : float
       minimum distance in NEMO grid boxes to look for a corresponding river outflow
       point

    """
    # Determine values used in this function. The initial values below are for
    # ORCA025 but they are then scaled up or down depending on the ocean resolution.
    max_value = int(round(32 * scaling))
    min_value = int(round(3 * scaling))
    st_laurance_sequence_number = int(round(34 * scaling))
    gradient = float(max_value - min_value) / float(st_laurance_sequence_number)

    # Determine the distance to search for your NEMO outflow point:
    # The St Laurance river TRIP outflow point is particularly far away from the
    # St Laurance estuary and requires looking 32 grid boxes away. Therefore we
    # will use linear interpolation to scale the distance to look from 32
    # (maximum needed for St Laurance) to 3 (minimum for all small rivers). This
    # distance is scaled by the river sequence number (which is 34 for the St
    # Laurance) using linear interpolation. This means the Ganges (which has a
    # sequence number of 28) will look 27 grid boxes away which is more than
    # enough (we used to use 12 grid boxes).
    min_distance = math.ceil(gradient * sequence_cube.data[j, i] + min_value)
    if min_distance > max_value:
        min_distance = max_value
    if min_distance < min_value:
        min_distance = min_value

    return min_distance


def find_nemo_point(
    direction_cube,
    nav_lat,
    nav_lon,
    scaling,
    river_index_array,
    orca_mask_array,
    sequence_cube,
    river_number_cube,
):
    """
    Finds the best NEMO point to link to your TRIP point. This isn't necessarily
    the nearest point as sometimes you have to pass over minor NEMO river
    outflow points to find the larger NEMO river outflow points to link to.

    Parameters
    ----------
    direction_cube : :class:`iris.cube.Cube`
       River direction numbers where 9 is an ocean outflow point and 10 is an inland
       basin outflow point
    nav_lat : :class:`numpy.ndarray`
       Array of latitudes on the NEMO 2D spatial grid
    nav_lon : :class:`numpy.ndarray`
       Array of longitudes on the NEMO 2D spatial grid
    scaling : float
       Scaling number to apply to distances.
       1 for ORCA025 and scales linearly with horizontal resolution
    river_index_array : :class:`numpy.ndarray`
       Array of river numbers on the NEMO 2D spatial grid
    orca_mask_array : :class:`numpy.ndarray`
       The land-sea mask: sea = 1, land = 0
    sequence_cube : :class:`iris.cube.Cube`
       River sequence numbers as output from ancil_river_routing.river_routing
    river_number_cube : :class:`iris.cube.Cube`
       The new river number assigned to this river
    """

    # Loop over TRIP points
    ni = direction_cube.data.shape[1]
    nj = direction_cube.data.shape[0]

    for i in range(ni):
        _LOGGER.info("Processing row ", i)
        for j in range(nj):

            if direction_cube.data[j, i] == 9 or direction_cube.data[j, i] == 10:

                # Find the latitude and longitde of your point
                latitude = direction_cube.coords("latitude")[0].points[j]
                longitude = direction_cube.coords("longitude")[0].points[i]
                if longitude > 180.0:
                    longitude = longitude - 360.0

                # Find the distance to your point from every other point
                distance_to_lat = (nav_lat - latitude) ** 2
                distance_to_lon = (nav_lon - longitude) ** 2
                distance = (distance_to_lat + distance_to_lon) ** 0.5

                # Find the closest point
                oned_index = distance.argmin()
                twod_index = np.unravel_index(oned_index, distance.shape)

                # Get the minimum distance to search for a corresponding NEMO outflow
                # point
                min_distance = get_min_distance(j, i, sequence_cube, scaling)

                # Perform a diamond search to get the river number
                # Use a min distance of 12 to cover a large enough area to find the
                # largest rivers.
                river_location = location_class.LocationClass(
                    twod_index[0],
                    twod_index[1],
                    distance.shape,
                )
                river_number_cube.data[j, i] = river_location.diamond_search(
                    river_index_array=river_index_array,
                    orca_mask_array=orca_mask_array,
                    min_distance=min_distance,
                )

                _LOGGER.info(j, i, min_distance, sequence_cube, river_number_cube)


def remove_nemo_rivers(river_index_array, river_number_cube):

    """
    Remove NEMO rivers not used in the coupling and decrease higher river numbers to
    fill the gaps

    Parameters
    ----------
    river_index_array : :class:`numpy.ndarray`
       Array of river numbers on the NEMO 2D spatial grid
    river_number_cube : :class:`iris.cube.Cube`
       Cube of river numbers on the TRIP 2D spatial grid
    """

    n_rivers = river_index_array.max()
    n_rivers_removed = 0
    for loop in range(1, int(n_rivers) + 1):
        river = loop - n_rivers_removed
        if not np.any(river_number_cube.data == river):
            where_river = river_index_array == river
            if np.any(where_river):
                river_index_array[where_river] = 0
                where_larger = river_index_array > river
                river_index_array[where_larger] = river_index_array[where_larger] - 1
                where_larger = river_number_cube.data > river
                river_number_cube.data[where_larger] = (
                    river_number_cube.data[where_larger] - 1
                )
                n_rivers_removed = n_rivers_removed + 1


def output_to_netcdf(
    ocean_river_number_file, river_index_array_out, nav_lat_orig, nav_lon_orig
):

    """
    Output river_index_array to netcdf file

    Parameters
    ----------
    ocean_river_number_file : str
       Output file for NEMO file (netcdf) containing river numbers on ORCA grid.
    river_index_array_out : :class:`numpy.ndarray`
       Array of river numbers on the NEMO grid. This has a halo of zeros surrounding it
       whilst the interior has been copied from river_index_array.
    nav_lat_orig : :class:`numpy.ndarray`
       Array of latitudes for the NEMO 2D spatial grid including original halo
    nav_lon_orig : :class:`numpy.ndarray`
       Array of longitudes for the NEMO 2D spatial grid including original halo
    """

    if ocean_river_number_file is not None:
        _LOGGER.info("Outputting to netcdf file ", ocean_river_number_file)
        output_nc = Dataset(ocean_river_number_file, "w")
        output_nc.createDimension("x", river_index_array_out.shape[1])
        output_nc.createDimension("y", river_index_array_out.shape[0])
        nav_lat_out_ncvar = output_nc.createVariable("nav_lat", "f8", ("y", "x"))
        nav_lon_out_ncvar = output_nc.createVariable("nav_lon", "f8", ("y", "x"))
        river_number_ncvar = output_nc.createVariable("river_number", "i4", ("y", "x"))
        nav_lat_out_ncvar[:, :] = nav_lat_orig
        nav_lon_out_ncvar[:, :] = nav_lon_orig
        river_number_ncvar[:, :] = river_index_array_out
        output_nc.close()


def main(
    direction_cube,
    sequence_cube,
    ocean_runoff_file=None,
    orca_dom_file=None,
    um_river_number_file=None,
    ocean_river_number_file=None,
):

    """
    Make two ancillary files (NEMO and UM) that contain the river numbers on both grids.

    Parameters
    ----------
    ocean_runoff_file : str
       Input file (netcdf) containing the climatological river outflow used in ocean
       only simulations
    direction_cube : :class:`iris.cube.Cube`
       River direction as output from ancil_river_routing.river_routing
    sequence_cube : :class:`iris.cube.Cube`
       River sequence numbers as output from ancil_river_routing.river_routing
    orca_dom_file : str
       Input file (netcdf) contining the NEMO domain on the ORCA grid
    um_river_number_file : str
       Output file for UM file (ancillary format) containing river numbers on TRIP grid.
    ocean_river_number_file : str
       Output file for NEMO file (netcdf) containing river numbers on ORCA grid.
    """

    # Load the mask. (land = 0, ocean = 1)
    rootgrp = Dataset(orca_dom_file, "r")
    tmask = rootgrp.variables["top_level"][0, :, :]
    rootgrp.close()

    # Generate an array of river indices for NEMO river outflows
    if ocean_runoff_file is not None:
        (
            river_index_array_orig,
            nav_lat_orig,
            nav_lon_orig,
        ) = generate_nemo_river_numbers(ocean_runoff_file, orca_dom_file=orca_dom_file)
    else:
        river_index_array_orig, nav_lat_orig, nav_lon_orig = blank_nemo_river_numbers(
            orca_dom_file
        )

    # Make a new cube to contain the river numbers on the TRIP grid
    river_number_cube = direction_cube.copy()
    river_number_cube.data[:, :] = 0

    # Shrink the NEMO data down to remove the 1 grid box halo
    river_index_array = river_index_array_orig[1:-1, 1:-1]
    nav_lat = nav_lat_orig[1:-1, 1:-1]
    nav_lon = nav_lon_orig[1:-1, 1:-1]
    orca_mask_array = tmask[1:-1, 1:-1]

    # Determine the scaling used to scale up/down the number of grid boxes to search
    # across
    scaling = scale_by_resolution(nav_lat, nav_lon)

    # Find the best NEMO outflow point to link to every TRIP outflow point
    find_nemo_point(
        direction_cube,
        nav_lat,
        nav_lon,
        scaling,
        river_index_array,
        orca_mask_array,
        sequence_cube,
        river_number_cube,
    )

    # Remove NEMO rivers not used in the coupling and decrease higher river numbers
    # to fill the gaps
    remove_nemo_rivers(river_index_array, river_number_cube)

    # Expand the arrays back to their original size
    river_index_array_out = np.zeros_like(river_index_array_orig)
    river_index_array_out[1:-1, 1:-1] = river_index_array

    # Output river_index_array to netcdf
    output_to_netcdf(
        ocean_river_number_file, river_index_array_out, nav_lat_orig, nav_lon_orig
    )

    # Change the cube meta data to make it suitable for river number (not river
    # direction)
    river_number_cube.var_name = "river_number"
    river_number_cube.long_name = "river_number"
    river_number_cube.attributes["STASH"] = iris.fileformats.pp.STASH(0o1, 00, 154)

    # Output river_number_cube to ancil format using ANTS
    save.ancil(river_number_cube, um_river_number_file)
    save.netcdf(river_number_cube, um_river_number_file)

    # If using ANTS it saves the wrong dimension size in the Integer Header

    _LOGGER.info("Finished main river number ancillary files")
