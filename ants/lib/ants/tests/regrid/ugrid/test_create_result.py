# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.
from unittest.mock import patch

import ants.tests
import iris
from ants.regrid._ugrid import _create_result as create_result


class TestAll(ants.tests.TestCase):
    def setUp(self):
        self.target = ants.tests.stock.mesh_C4()
        self.source = ants.tests.stock.geodetic((2, 2))
        self.mesh_attributes = [
            "Conventions",
            "mesh_topology",
            "face_node_connectivity",
            "nodes",
        ]

    @patch("ants.regrid._ugrid._extrude_cube")
    def test_extrude_cube_call(self, mock_extrude):
        create_result(self.source, self.target)
        mock_extrude.assert_called_once_with(self.target, [96], self.source.dtype)

    def test_lat_from_target(self):
        expected = self.target.coord("latitude")

        actual = create_result(self.source, self.target).coord("latitude")

        self.assertEqual(expected, actual)

    def test_lon_from_target(self):
        expected = self.target.coord("longitude")

        actual = create_result(self.source, self.target).coord("longitude")

        self.assertEqual(expected, actual)

    def test_scalar_coordinates(self):
        self.source.add_aux_coord(iris.coords.AuxCoord(0, long_name="foo"))
        expected = self.source.coord("foo")

        result = create_result(self.source, self.target)
        actual = result.coord("foo")

        self.assertEqual(actual, expected)

    def test_additional_coordinates(self):
        source = ants.tests.stock.simple_3d_time_varying()
        expected = source.coord(axis="T")

        result = create_result(source, self.target)
        actual = result.coord(axis="T")

        self.assertEqual(actual, expected)

    def test_shape_with_multidimensional_source(self):
        expected = (3, 96)

        source = ants.tests.stock.simple_3d_time_varying()
        actual = create_result(source, self.target).shape

        self.assertEqual(expected, actual)

    def test_attributes_describing_mesh(self):
        # Results dict small enough that hiding it on failure is
        # counter-productive:
        self.maxDiff = None
        # Attributes that are needed to describe the mesh come from the target:

        expected = {
            attribute: self.target.attributes[attribute]
            for attribute in self.mesh_attributes
        }

        result = create_result(self.source, self.target)
        actual = {
            attribute: result.attributes[attribute]
            for attribute in self.mesh_attributes
        }

        self.assertEqual(expected, actual)

    def test_nonmesh_metadata_from_source(self):
        # metadata includes name, units and attributes.  Mesh attributes need
        # to come from target, while other attributes come from the source.
        self.source.attributes["foo"] = "bar"
        expected = self.source.metadata

        actual = create_result(self.source, self.target).metadata
        [
            actual.attributes.pop(attribute)
            for attribute in self.mesh_attributes
            if attribute in actual.attributes
        ]

        self.assertEqual(expected, actual)

    def test_contradictory_attributes(self):
        # In future, this may become combining attributes.  For now, though,
        # same attribute on source and target will throw exception:
        self.source.attributes["mesh_topology"] = "bar"
        with self.assertRaisesRegex(ValueError, "One or more attributes are"):
            create_result(self.source, self.target)


if __name__ == "__main__":
    ants.tests.main()
