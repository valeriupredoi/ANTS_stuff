# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.
import ants.tests
import iris
import numpy as np
from proc_ants.order_nemo_rivers import find_nemo_point


class TestAll(ants.tests.TestCase):
    def test_all(self):

        # Generate some rivers flowing east
        river_direction = np.array(
            [
                [0, 0, 0, 0, 0],
                [3, 3, 9, 0, 0],
                [3, 3, 9, 0, 0],
                [3, 4, 0, 0, 0],
                [0, 0, 9, 0, 0],
            ]
        )

        lat_coord = iris.coords.DimCoord(
            [-0.2, -0.1, 0.0, 0.1, 0.2], standard_name="latitude"
        )
        lon_coord = iris.coords.DimCoord(
            [89.8, 89.9, 90.0, 90.1, 90.2], standard_name="longitude"
        )
        direction_cube = iris.cube.Cube(
            river_direction,
            long_name="river direction",
            dim_coords_and_dims=[(lat_coord, 0), (lon_coord, 1)],
        )

        nav_lat = np.array(
            [
                [0.2, 0.2, 0.2, 0.2, 0.2],
                [0.1, 0.1, 0.1, 0.1, 0.1],
                [0.0, 0.0, 0.0, 0.0, 0.0],
                [-0.1, -0.1, -0.1, -0.1, -0.1],
                [-0.2, -0.2, -0.2, -0.2, -0.2],
            ]
        )

        nav_lon = np.array(
            [
                [89.8, 89.9, 90.0, 90.1, 90.2],
                [89.8, 89.9, 90.0, 90.1, 90.2],
                [89.8, 89.9, 90.0, 90.1, 90.2],
                [89.8, 89.9, 90.0, 90.1, 90.2],
                [89.8, 89.9, 90.0, 90.1, 90.2],
            ]
        )

        scaling = 1.0

        river_index_array = np.array(
            [
                [4, 0, 0, 0, 0],
                [4, 2, 2, 0, 0],
                [4, 2, 2, 0, 0],
                [4, 5, 6, 0, 0],
                [4, 5, 6, 0, 0],
            ]
        )

        orca_mask_array = np.array(
            [
                [1, 0, 0, 0, 0],
                [1, 1, 1, 0, 0],
                [1, 1, 1, 0, 0],
                [1, 1, 1, 0, 0],
                [1, 1, 1, 0, 0],
            ]
        )

        river_sequence = np.array(
            [
                [0, 0, 0, 0, 0],
                [7, 6, 5, 0, 0],
                [9, 8, 7, 0, 0],
                [4, 3, 0, 0, 0],
                [0, 0, 2, 0, 0],
            ]
        )

        sequence_cube = iris.cube.Cube(river_sequence, long_name="river sequence")

        river_number_cube = iris.cube.Cube(river_sequence, long_name="river number")

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

        target_river_number_array = np.array(
            [
                [0, 0, 0, 0, 0],
                [7, 6, 2, 0, 0],
                [9, 8, 2, 0, 0],
                [4, 3, 0, 0, 0],
                [0, 0, 2, 0, 0],
            ]
        )

        self.assertArrayEqual(river_number_cube.data, target_river_number_array)


if __name__ == "__main__":
    ants.tests.main()


if __name__ == "__main__":
    ants.tests.main()
