# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.

"""
Location class
********************

Based on a location (j, i) you can see what grid boxes are nearby, what type of
grid boxes they are (land, sea or coastal) and search outwards from your location.

Functions are:

* __init__ = Initialises this class at a grid point that is provided (j,i)
* east,west,north,south = Moves one grid box in this direction
* northeast, northwest, southeast, southwest = Moves diagonally in this direction
* shift = Moves in a direction specified by input argument
* __str__ = Controls how the location class appears when printed using the print
  function centred on your location.
* is_coastal = Returns true if this is a coastal land point
* find_neighbouring_sea = Find all the sea points a certain distance away as long as
  they follow the coastline.
* diamond_search = Gradually searches outwards in a diamond shape to find the nearest
  river outflow point or coastal land point.

"""
import logging

_LOGGER = logging.getLogger(__name__)


class LocationClass:

    """
    Class to hold information on a grid box and routines to get the neighbouring
    grid boxes

    i indices are longitude indices of a 2D spatial point
    j indices are latitude indices of a 2D spatial point

    Some routines held within this class allow you to move your indices in different
    directions. In relation to point x, these directions are:

        8  1  2
        7  x  3
        6  5  4

     e.g. 3 is the box to the east


    """

    def __init__(self, j, i, shape, amount=1.0, river_number=0):
        """
        Initialises an instance of this class

        Parameters
        ----------
        j : int
           latitude index of data point to initialise at
        i : int
           longitude index of data point to initialise at
        shape : list
           shape [ny, nx] of data array
        amount : float
           Any quantity that is associated with this point.
           For river routing applications this is the river outflow amount.
        river_number : int
           For river routing applications this is the river index number

        """

        self.i = i  # i location (longitude index)
        self.j = j  # j location (latitude index)
        self.j_overlap = j  # j_overlap is allowed to go outside the domain
        self.i_original = i  # original i location (will not change)
        self.j_original = j  # original j location (will not change)
        self.ni = shape[1]  # size of i axis
        self.nj = shape[0]  # size of j axis
        self.shape = shape  # shape of domain
        self.amount = amount  # quantity of field at this point
        self.river_number = river_number  # river index number

    def east(self):
        """
        Moves one grid box to the east
        """
        self.i = self.i + 1
        if self.i > self.ni - 1:
            self.i = 0

    def west(self):
        """
        Moves one grid box to the west

        Returns
        -------
        newj : int
           New j index
        newi : int
           New i index

        """
        self.i = self.i - 1
        if self.i < 0:
            self.i = self.ni - 1

    def north(self, amount=1):
        """
        Moves one or more grid box to the north

        Parameters
        ----------
        amount : int
           The number of grid boxes to move by

        """
        self.j_overlap = self.j_overlap + amount
        if self.j_overlap > self.nj - 1:
            self.j = self.nj - 1
        elif self.j_overlap < 0:
            self.j = 0
        else:
            self.j = self.j_overlap

    def south(self, amount=1):
        """
        Moves one or more grid box to the south

        Parameters
        ----------
        amount : int
           The number of grid boxes to move by

        """
        self.j_overlap = self.j_overlap - amount
        if self.j_overlap > self.nj - 1:
            self.j = self.nj - 1
        elif self.j_overlap < 0:
            self.j = 0
        else:
            self.j = self.j_overlap

    def northeast(self):
        """
        Moves one grid box to the northeast

        """
        self.north()
        self.east()

    def northwest(self):
        """
        Moves one grid box to the northwest

        """
        self.north()
        self.west()

    def southeast(self):
        """
        Moves one grid box to the southeast

        """
        self.south()
        self.east()

    def southwest(self):
        """
        Moves one grid box to the southwest

        """
        self.south()
        self.west()

    def shift(self, direction):
        """
        Return the coordinates of the gridbox shifted by one in one of the following
        directions:

        8  1  2
        7  x  3
        6  5  4

        e.g. 3 is the box to the east

        Parameters
        ----------
        direction : int
           The direction to move in (see grid above)

        """

        if direction == 1:
            self.north()
        elif direction == 2:
            self.northeast()
        elif direction == 3:
            self.east()
        elif direction == 4:
            self.southeast()
        elif direction == 5:
            self.south()
        elif direction == 6:
            self.southwest()
        elif direction == 7:
            self.west()
        elif direction == 8:
            self.northwest()
        else:
            raise Exception("Only direction 1 to 8 allowed")

    def revert(self):
        """
        Returns the location back to its original i and j coordinates
        """

        self.i = self.i_original
        self.j = self.j_original
        self.j_overlap = self.j_original

    def __str__(self):
        """
        Controls how the location class appears when printed using the print function

        """

        return (
            "Catchment at "
            + str(self.j)
            + ","
            + str(self.i)
            + " has amount of "
            + str(self.amount)
            + " and river number "
            + str(self.river_number)
        )

    def is_coastal(self, orca_mask_array):
        """
        Returns true if this is a coastal land point

        Parameters
        ----------
        orca_mask_array : :class:`numpy.ndarray`
           The land-sea mask: sea = 1, land = 0

        Returns
        -------
        is_coastal : bool
           Returns one of the following:
              A land point with sea next to it = True
              A land point with no sea next to it = False
              A sea point = False
        """

        # If it is a sea point return false
        if orca_mask_array[self.j, self.i] > 0.5:
            return False

        is_coastal = False
        for direction in range(1, 9):

            # Find neighbouring point
            oldj, oldi, oldj_overlap = self.j, self.i, self.j_overlap
            self.shift(direction)
            newj, newi = self.j, self.i

            # If this is a sea point return True
            if orca_mask_array[newj, newi] > 0.5:
                is_coastal = True

            # Revert back to what we had before the direction shift as we don't
            # want to affect the diamond_search that called this routine
            self.j, self.i, self.j_overlap = oldj, oldi, oldj_overlap

        return is_coastal

    def find_neighbouring_sea(self, orca_mask_array, max_distance):
        """
        Find all the sea points a certain distance away as long as they follow the
        coastline.

        Parameters
        ----------
        orca_mask_array : :class:`numpy.ndarray`
           Array of NEMO land sea mask (land = 0, ocean = 1)
        max_distance : int
           Maximum distance out to do the search

        Returns
        -------
        sea_coord_list : list
           A list of coordinates where each entry is of the form (j,i).
           Each coordinate is a coastal sea point.

        """

        # Initialise lists
        full_land_coord_list = [(self.j, self.i)]
        next_land_coord_list = [(self.j, self.i)]
        sea_coord_list = []

        for step in range(max_distance):

            # Clear out the next list of land coordinates ready for filling
            new_land_coord_list = next_land_coord_list
            next_land_coord_list = []

            for land_coord in new_land_coord_list:
                land_coord_location = LocationClass(
                    land_coord[0], land_coord[1], self.shape
                )
                for direction in range(1, 9):

                    # Find neighbouring point
                    land_coord_location.shift(direction)
                    newj, newi = land_coord_location.j, land_coord_location.i
                    land_coord_location.revert()

                    # If it is a sea point add it to the list that will be returned
                    if orca_mask_array[newj, newi] > 0.5:
                        if (newj, newi) not in sea_coord_list:
                            sea_coord_list.append((newj, newi))

                    # If it is a land point add it to the list which we will loop
                    # over next time
                    else:
                        if (newj, newi) not in full_land_coord_list:
                            full_land_coord_list.append((newj, newi))
                            next_land_coord_list.append((newj, newi))

            # If we have run out of land points then break out
            if len(next_land_coord_list) == 0:
                break

        if len(sea_coord_list) == 0:
            _LOGGER.warning(
                "Could not find any sea coordinates near land point at i=",
                self.i,
                " j=",
                self.j,
            )

        return sea_coord_list

    def diamond_search(
        self, river_index_array=None, orca_mask_array=None, min_distance=None
    ):
        """
        Gradually searches outwards in a diamond shape to find the best NEMO river
        outflow point or coastal land point to link to your river. As the code
        travels outwards in a diamond shape it adds all NEMO river outflow points
        and NEMO coastal points to a pair of lists (outflow_point_list and
        land_point_list). If an outflow point is found this should already have a
        river number so return that. If a coastal point is found then generate a
        new set of river outflow points along the coast with a new river number
        and return that. The diamond search alternates between starting its
        diamond shape to the north and to the south of the original location. This
        is to stop the program having either a northward or southward bias.

        A diamond shape was chosen over a square shaped spiral search as the
        majority of river misalignments are due to the rivers being shifted in
        either north, south, east or west directions and not the diagonal
        directions (i.e. northeast). Therefore a spiral search that looks further
        in the north, south, east or west directions instead of the diagonal
        directions is better (i.e. a diamond search).

        Parameters
        ----------
        river_index_array : :class:`numpy.ndarray`
           Array of river numbers on the NEMO grid. All river outflow points
           have a river number greater than 0.
        orca_mask_array : :class:`numpy.ndarray`
           Array of NEMO land sea mask (land = 0, ocean = 1)
        min_distance : int
           Minimum distance out to do the search

        Returns
        -------
        river_number : int
           The river number that has been found from this search
        """

        outflow_point_list = []
        land_point_list = []
        go_north = True  # Start off by going north to start your diamond

        # To summarise the following algorithm, consider two iterations of
        # "step", searching a 6x6 array starting from "X":

        # [[ 0  0  0  0  8  0  0]
        #  [ 0  0  0  7  4  9  0]
        #  [ 0  0  6  3  X  1 10]
        #  [ 0  0  0  5  2 11  0]
        #  [ 0  0  0  0 12  0  0]
        #  [ 0  0  0  0  0  0  0]
        #  [ 0  0  0  0  0  0  0]]

        # For the first iteration of "step", the search is shifted one grid
        # box north.  From here, we enter the iterations of "ministep", and
        # the first grid box searched is the grid box to the SE of the shifted
        # start position (i.e. the grid box labelled 1).  The second iteration
        # of ministep goes SW from here to grid box 2, and so on through to
        # grid box 4.
        #
        # The next iteration of "step" goes 2 grid boxes south from the initial
        # position, and the "ministep" iterations go clockwise, searching grid
        # boxes 5 through 12.
        #
        # This pattern then repeats - i.e. shifts south 3 grid box searching
        # counterclockwise, then south 4 grid boxes searching clockwise etc.
        #
        # The search wraps around the east/west boundary, and is truncated at
        # the north/south boundary:

        # [[ 1  10 0  0  6  3  X]
        #  [ 11 0  0  0  0  5  2]
        #  [ 0  0  0  0  0  0 12]
        #  [ 0  0  0  0  0  0  0]
        #  [ 0  0  0  0  0  0  0]
        #  [ 0  0  0  0  0  0  0]
        #  [ 0  0  0  0  0  0  0]]
        #
        # Because of the truncation at the north boundary, some grid boxes are
        # searched multiple times.  The 4th grid box searched is X; and the 7th,
        # 8th, and 9th grid boxes searched are 3, X and 1.

        for step in range(1, 1000):

            # Go back to the original point for the next iteration of the diamond
            self.revert()

            # Go north or south
            if go_north:
                self.north(amount=step)
                go_north = False  # Make it go south on the next step
                direction_order = [
                    4,
                    6,
                    8,
                    2,
                ]  # The order of directions to go round the diamond
            else:
                self.south(amount=step)
                go_north = True  # Make it go north on the next step
                direction_order = [8, 2, 4, 6]

            # Then go diagonal the step size in each of the 4 diagonal directions
            for direction in direction_order:
                for mini_step in range(1, step + 1):
                    self.shift(direction)

                    # Have we found a river outflow point. If so add it to the list.
                    if river_index_array[self.j, self.i] > 0:
                        outflow_point_list.append((self.j, self.i))

                    # Have we found a coastal land point. If so add it to the list.
                    if self.is_coastal(orca_mask_array):
                        land_point_list.append((self.j, self.i))

            # Break out if we have found something (either an outflow point or a
            # coastal land point) and we have done the required number of loops
            if step >= min_distance:
                if len(outflow_point_list) > 0:

                    # Return the river number you have found
                    river_number_list = [
                        river_index_array[outflow_point]
                        for outflow_point in outflow_point_list
                    ]
                    river_number = min(river_number_list)
                    break

                if len(land_point_list) > 0:

                    # Find any neighbouring sea points
                    coord_list = []
                    for land_point in land_point_list:
                        river_location = LocationClass(
                            land_point[0],
                            land_point[1],
                            self.shape,
                        )
                        this_coord_list = river_location.find_neighbouring_sea(
                            orca_mask_array, 3
                        )
                        coord_list.extend(this_coord_list)

                    # If we can't find any sea points then assume this is an inland
                    # basin flow point and assign this as a river number of -2
                    if len(coord_list) == 0:
                        error_str = (
                            "WARNING: Could not find any sea coordinates near any "
                            + "nearby land point i="
                            + repr(self.i)
                            + " j="
                            + repr(self.j)
                            + ". Assuming this is an inland basin flow point and "
                            + "assigning as river_number=-2. If any water tries to "
                            + "flow out to the ocean at this point it will make "
                            + "the model fail."
                        )
                        _LOGGER.error(error_str)
                        river_number = -2
                        break

                    # Remove duplicates
                    coord_list = list(set(coord_list))

                    # Do the neighbouring sea points already have a river number
                    # assigned to them. If so we will use that.
                    river_number_list = []
                    for coord in coord_list:
                        if river_index_array[coord] > 0:
                            river_number_list.append(river_index_array[coord])
                    if len(river_number_list) > 0:
                        river_number = min(river_number_list)
                        break

                    # If we have not found a river then we need to make one
                    new_outflow_count = 0
                    river_number = river_index_array.max() + 1

                    for coord in coord_list:
                        if river_index_array[coord] == 0:
                            river_index_array[coord] = river_number
                            new_outflow_count = new_outflow_count + 1

                    if new_outflow_count == 0:
                        error_str = "Could not find any sea points for outflow"
                        _LOGGER.error(error_str)
                        raise Exception(error_str)

                    break

        return int(river_number)
