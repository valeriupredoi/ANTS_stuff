# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.
import ants.tests
import numpy as np
from proc_ants.location_class import LocationClass
from proc_ants.order_nemo_rivers import get_neighbouring_points


class TestAll(ants.tests.TestCase):
    def assert_location(self, final_location_list):

        # Initialise lists for end coordinates
        ej_list = []
        ei_list = []

        for location in final_location_list:
            ej_list.append(location.j)
            ei_list.append(location.i)

        # Convert to numpy arrays
        ej_array = np.array(ej_list)
        ei_array = np.array(ei_list)

        # Do the tests
        self.assertArrayEqual(ej_array, np.array([3, 3, 5, 6]))
        self.assertArrayEqual(ei_array, np.array([6, 5, 3, 3]))

    def test_all(self):

        # Array of river runoff
        runoff_nc_array = np.array(
            [
                [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                [0.0, 0.1, 0.1, 0.0, 0.0, 0.0, 0.0],
                [0.0, 0.1, 0.1, 0.0, 0.0, 0.0, 0.3],
                [0.0, 0.0, 0.0, 0.0, 0.0, 0.3, 0.3],
                [0.0, 0.0, 0.0, 0.0, 0.0, 0.3, 0.3],
                [0.0, 0.0, 0.2, 0.2, 0.2, 0.3, 0.3],
                [0.0, 0.2, 0.2, 0.2, 0.3, 0.3, 0.3],
            ]
        )

        # What will be the starting location. We start one away from in the
        # bottom right corner.
        sj = 5
        si = 5
        river_number = 1

        # Array of river numbers (starts off all zeros)
        river_index_array = np.array(
            [
                [0, 0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 0],
            ]
        )

        # Populate the starting position with the river number
        river_index_array[sj, si] = river_number

        # Set the starting location as an item for the location_list
        start_location_list = [
            LocationClass(
                sj,
                si,
                runoff_nc_array.shape,
                amount=runoff_nc_array[sj, si],
                river_number=river_number,
            )
        ]

        # Call get_neighbouring_points twice to fully test the routine
        interim_location_list = get_neighbouring_points(
            runoff_nc_array, start_location_list, river_index_array
        )
        final_location_list = get_neighbouring_points(
            runoff_nc_array, interim_location_list, river_index_array
        )

        # Target array of river numbers (once modified) All locations with 0.3
        # runoff have been identified as river number 1 apart from the
        # northernmost location
        target_river_index_array = np.array(
            [
                [0, 0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 1, 1],
                [0, 0, 0, 0, 0, 1, 1],
                [0, 0, 0, 1, 1, 1, 1],
                [0, 0, 0, 1, 1, 1, 1],
            ]
        )

        self.assertArrayEqual(river_index_array, target_river_index_array)
        self.assert_location(final_location_list)


if __name__ == "__main__":
    ants.tests.main()
