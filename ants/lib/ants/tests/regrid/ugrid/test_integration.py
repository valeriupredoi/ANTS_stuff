# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.
import ants.tests
import iris
import numpy as np
from ants.regrid import GeneralRegridScheme


class _Common(object):
    def setUp(self):
        self.mesh = ants.tests.stock.mesh_C4()
        scheme = GeneralRegridScheme(horizontal_scheme="_UGrid")
        self.result = self.grid.regrid(self.mesh, scheme)
        self.expected_target = self.mesh

    def test_returns_cube(self):
        expected = iris.cube.Cube

        actual = self.result

        self.assertIsInstance(actual, expected)

    def test_result_is_mesh(self):
        expected = "UGRID-1.0"

        actual = self.result.attributes["Conventions"]

        self.assertIn(expected, actual)

    def test_returned_latitudes(self):
        expected = self.mesh.coord("latitude")

        actual = self.result.coord("latitude")

        self.assertEqual(expected, actual)

    def test_returned_longitudes(self):
        expected = self.mesh.coord("longitude")

        actual = self.result.coord("longitude")

        self.assertEqual(expected, actual)

    def test_target_is_unchanged(self):
        expected = self.expected_target

        actual = self.mesh

        self.assertEqual(expected.metadata, actual.metadata)
        self.assertMaskedArrayEqual(expected.data, actual.data)

    def test_returned_data(self):
        # See subclasses for how expected data is derived.
        expected = self.expected_data

        data = self.result.data
        actual = np.sum(data[~np.isnan(data)])

        self.assertAlmostEqual(expected, actual)

    def test_attributes(self):
        expected = self.mesh.attributes

        actual = self.result.attributes

        self.assertEqual(expected, actual)

    def test_source_is_unchanged(self):
        expected = self.expected_source

        actual = self.grid

        self.assertEqual(expected, actual)

    def test_all_coordinates_from_source_on_target(self):
        # TODO: when reference surfaces are also regridded (#1070), the
        # reference surface and dependent factories need to be added to the
        # result too.
        not_handled_yet = ("altitude", "surface_altitude")
        expected = sorted(
            [
                coord.name()
                for coord in self.grid.coords()
                if coord.name() not in not_handled_yet
            ]
        )

        actual = sorted([coord.name() for coord in self.result.coords()])
        self.assertEqual(expected, actual)


@ants.tests.skip_esmf
class Test4DGridToMesh(_Common, ants.tests.TestCase):
    def setUp(self):
        self.grid = ants.tests.stock.simple_4d_with_hybrid_height()
        self.grid.data = np.float32(self.grid.data)
        super(Test4DGridToMesh, self).setUp()
        # Common factor of three levels, so expected result is sum of
        # simple_4d_with_hybrid_height (64,620), divided by numer of cells per
        # level in that cube (30) multiplied by number of cells per level in
        # C4 mesh (96), i.e. (206,784)
        self.expected_data = 206784
        self.expected_source = ants.tests.stock.simple_4d_with_hybrid_height()


@ants.tests.skip_esmf
class Test4DGridTo3DMesh(_Common, ants.tests.TestCase):

    # Deliberately does not call superclass setup
    def setUp(self):
        template = ants.tests.stock.mesh_C4()
        self.grid = ants.tests.stock.simple_4d_with_hybrid_height()
        self.grid.data = np.float32(self.grid.data)
        non_horizontal_dims = self.grid.shape[:-2]
        mesh_shape = non_horizontal_dims + template.shape
        self.mesh = iris.cube.Cube(np.ones(mesh_shape))

        # Add horizontal coords from template
        for coord in template.coords(dim_coords=False):
            coord_dim = len(self.mesh.shape) - 1
            self.mesh.add_aux_coord(coord, coord_dim)

        # Add non-horizontal coords from grid
        for coord in self.grid.coords(dim_coords=True, dimensions=non_horizontal_dims):
            try:
                coord_dim = self.grid.coord_dims(coord)[0]
            except IndexError:
                # scalar coordinate - ignore for now
                continue
            self.mesh.add_dim_coord(coord, coord_dim)
        for coord in self.grid.coords(dim_coords=False, dimensions=non_horizontal_dims):
            try:
                coord_dim = self.grid.coord_dims(coord)[0]
            except IndexError:
                # scalar coordinate - ignore for now
                continue
            self.mesh.add_aux_coord(coord, coord_dim)

        # And finally add UGrid topology attributes from template
        preserved_attributes = [
            "mesh_topology",
            "face_node_connectivity",
            "Conventions",
            "nodes",
        ]
        for key in preserved_attributes:
            self.mesh.attributes[key] = template.attributes[key]

        scheme = GeneralRegridScheme(horizontal_scheme="_UGrid")
        self.result = self.grid.regrid(self.mesh, scheme)
        # Common factor of three levels, so expected result is sum of
        # simple_4d_with_hybrid_height (64,620), divided by numer of cells per
        # level in that cube (30) multiplied by number of cells per level in
        # C4 mesh (96), i.e. (206,784)
        self.expected_data = 206784
        self.expected_source = ants.tests.stock.simple_4d_with_hybrid_height()
        self.expected_target = self.mesh


@ants.tests.skip_esmf
class TestSingleValueDimension(_Common, ants.tests.TestCase):
    def setUp(self):
        # UKCA data has a single valued model level number dimension for
        # surface fields
        grid = iris.cube.Cube(np.ones((1, 5, 6)))
        ref = ants.tests.stock.simple_3d_time_varying()
        grid.add_aux_coord(ref.coord(dimensions=2), 2)
        grid.add_aux_coord(ref.coord(dimensions=1), 1)
        grid.add_aux_coord(
            iris.coords.DimCoord(
                0,
                standard_name="model_level_number",
            ),
            0,
        )
        self.grid = grid
        super(TestSingleValueDimension, self).setUp()
        self.expected_data = 96
        self.expected_source = grid


@ants.tests.skip_esmf
class Test4DAnonymousDimensionGridToMesh(_Common, ants.tests.TestCase):

    # Some netCDF source files don't have a dimension coordinate variable, but
    # just an auxcoord along the dimension (e.g. pseudo levels).  This is
    # similar to how we're handling lat/lon for UGrid files, but needs to be
    # handled in our regrid:
    def setUp(self):
        self.grid = ants.tests.stock.simple_4d_with_hybrid_height()
        self.grid.data = np.float32(self.grid.data)
        tc = iris.coords.AuxCoord.from_coord(self.grid.coord("time"))
        self.grid.remove_coord("time")
        self.grid.add_aux_coord(tc, data_dims=0)
        super(Test4DAnonymousDimensionGridToMesh, self).setUp()
        # Common factor of three levels, so expected result is sum of
        # simple_4d_with_hybrid_height (64,620), divided by numer of cells per
        # level in that cube (30) multiplied by number of cells per level in
        # C4 mesh (96), i.e. (206,784)
        self.expected_data = 206784
        self.expected_source = ants.tests.stock.simple_4d_with_hybrid_height()


@ants.tests.skip_esmf
class TestMasked2DGridToMesh(_Common, ants.tests.TestCase):
    # Extra tests in this class to account for masking

    def setUp(self):
        # Need 5 test regimes:
        # 1. source and target masked
        # 2. source masked, target unmasked
        # 3. source unmasked, target masked
        # 4. source and target unmasked
        # 5. partially masked source, target unmasked
        # Each test regime will have constant data to simplify testing(but
        # leave data outside of test regimes as standard geodetic)

        self.grid = self._prepare_source(ants.tests.stock.geodetic((12, 12)))
        self.mesh = self._mask_target(ants.tests.stock.mesh_C4())
        scheme = GeneralRegridScheme(horizontal_scheme="_UGrid")

        # Expected data is empirical rather than having been derived - it's a
        # KGO, effectively.
        self.expected_data = 5804.8479298136963
        self.result = self.grid.regrid(self.mesh, scheme)
        self.expected_target = self.mesh
        self.expected_source = self._prepare_source(ants.tests.stock.geodetic((12, 12)))

    @staticmethod
    def _prepare_source(cube):
        cube.data = np.ma.masked_array(cube.data)
        cube.data[7:9, 4:6] = 1  # Region 1
        cube.data[7:9, 4:6] = np.ma.masked  # Region 1
        cube.data[7:9, 6:8] = 2  # Region 2
        cube.data[7:9, 6:8] = np.ma.masked  # Region 2
        cube.data[3:5, 4:6] = 3  # Region 3
        cube.data[3:5, 6:8] = 4  # Region 4
        cube.data[7:9, 9:11] = 5  # Region 5
        cube.data[8, 9] = np.ma.masked  # Region 5
        return cube

    @staticmethod
    def _mask_target(cube):
        cube.data = np.ma.masked_array(cube.data)
        cube.data[0] = np.ma.masked  # Region 1
        cube.data[12] = np.ma.masked  # Region 3
        return cube

    def test_region_one_masked(self):
        self.assertTrue(np.ma.is_masked(self.result.data[0]))

    def test_region_two_unmasked(self):
        self.assertFalse(np.ma.is_masked(self.result.data[3]))

    def test_region_three_masked(self):
        self.assertTrue(np.ma.is_masked(self.result.data[12]))

    def test_region_four_unmasked(self):
        self.assertFalse(np.ma.is_masked(self.result.data[15]))

    def test_region_five_unmasked(self):
        self.assertFalse(np.ma.is_masked(self.result.data[19]))

    # Only test data for unmasked regions:
    def test_region_two_data(self):
        # All source cells are masked, want a NaN (or masked) output:
        actual = self.result.data[3]

        self.assertTrue(np.isnan(actual))

    def test_region_four_data(self):
        expected = 4

        actual = self.result.data[15]

        self.assertEqual(expected, actual)

    def test_region_five_data(self):
        expected = 5

        actual = self.result.data[19]

        self.assertEqual(expected, actual)


if __name__ == "__main__":
    ants.tests.main()
