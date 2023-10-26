# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.
import collections
import difflib
import os
import subprocess
import time
import unittest
import warnings
from contextlib import contextmanager
from io import StringIO

import ants
import iris
import mule
import numpy as np
import numpy.testing
from ants.config import CONFIG

from . import stock

try:
    from um_utils import pumf

    _PUMF_IMPORT_ERROR = False
except Exception as err:
    _PUMF_IMPORT_ERROR = err


__all__ = ["stock"]


# FIXME: Iris does not expose a public interface for controlling the output of
# its test methods.
_RESULT_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")


_RESOURCE_PATH = os.path.join(os.path.split(__file__)[0], "resources")


SKIP_OPTIONAL_TESTS = True


# We are dependent on some iris private methods for now, so we set the iris
# environmental variable based on the ANTS one, during this transition away
# from the iris testing code.
os.environ["IRIS_TEST_CREATE_MISSING"] = os.environ.get("ANTS_TEST_CREATE_MISSING", "")


#: Default perceptual hash size.
_HASH_SIZE = 16
#: Default maximum perceptual hash hamming distance.
_HAMMING_DISTANCE = 0


def get_data_path(relative_path):
    """
    Given the test data resource, returns the full path to the file.

    This should not be needed often in tests - but there are cases where it's
    required.

    """
    if not isinstance(relative_path, str):
        relative_path = os.path.join(*relative_path)
    return os.path.abspath(os.path.join(_RESOURCE_PATH, relative_path))


def _skip_importable(module, name):
    skip = unittest.skipIf(
        condition=not module and SKIP_OPTIONAL_TESTS,
        reason="Test requires '{}'.".format(name),
    )
    return skip


def skip_gdal(fn):
    """
    Decorator to choose whether to run tests, based on the availability of the
    libgdal library.

    Example usage:
        @skip_gdal
        class MygdalTest(test.IrisTest):
            ...

    """
    return _skip_importable(ants.fileformats.raster.gdal, "gdal")(fn)


def skip_f90nml(fn):
    """
    Decorator to choose whether to run tests, based on the availability of the
    f90nml library.

    Example usage:
        @skip_f90nml
        class Myf90nmlTest(test.IrisTest):
            ...

    """
    return _skip_importable(ants.fileformats.namelist.f90nml, "f90nml")(fn)


def skip_stratify(fn):
    """
    Decorator to choose whether to run tests, based on the availability of the
    stratify library.

    Example usage:
        @skip_stratify
        class MyStratifyTest(test.IrisTest):
            ...

    """
    return _skip_importable(ants.regrid.interpolation.stratify, "stratify")(fn)


def skip_esmf(fn):
    """
    Decorator to choose whether to run tests, based on the availability of the
    ESMPy library.

    Example usage:
        @skip_esmf
        class MyESMFTests(test.IrisTest):
            ...

    """
    return _skip_importable(ants.regrid.esmf.ESMF, "ESMF")(fn)


def skip_spiral(fn):
    """
    Decorator to choose whether to run tests, based on the availability of the
    compiled 'spiral' search.

    Example usage:
        @skip_spiral
        class MySprialTests(test.IrisTest):
            ...

    """
    return _skip_importable(ants.analysis._merge.spiral, "spiral")(fn)


def _method_path(meth, cls):
    return ".".join([cls.__module__, cls.__name__, meth.__name__])


class TestCase(unittest.TestCase):
    _assertion_counts = collections.defaultdict(int)

    @staticmethod
    def assertArrayEqual(actual, expected, error_message="", verbose=True):
        """
        Test that two numpy arrays are equal.

        Consult the `numpy testing docs
        <https://numpy.org/doc/stable/reference/routines.testing.html#test-support-numpy-testing>`_
        for more details.

        Parameters
        ----------
        actual : np.ma.masked_array
            The first array to compare.
        expected : np.ma.masked_array
            The second array to compare.
        error_message : str, optional
            Message to print on test failure
        verbose : bool, optional
            If True, prints differences between the arrays.

        """
        numpy.testing.assert_array_equal(
            actual, expected, err_msg=error_message, verbose=verbose
        )

    @staticmethod
    def assertArrayAlmostEqual(
        actual, expected, decimal=6, error_message="", verbose=True
    ):
        """
        Test that two numpy arrays are equal.

        Consult the `numpy testing docs
        <https://numpy.org/doc/stable/reference/routines.testing.html#test-support-numpy-testing>`_
        for more details.

        The behaviour of this test may change in a future version of Ants,
        when the numpy version is upgraded.

        Parameters
        ----------
        actual : np.ma.masked_array
            The first array to compare.
        expected : np.ma.masked_array
            The second array to compare.
        decimal : int, optional
            Precision for comparison, defaults to 6.
        error_message : str, optional
            Message to print on test failure
        verbose : bool, optional
            If True, prints differences between the arrays.

        """
        # TODO: should be updated to assert_array_almost_equal_nulp.  Using
        # assert_array_almost_equal for now as a stepping stone.
        # https://code.metoffice.gov.uk/trac/ancil/ticket/1480
        numpy.testing.assert_array_almost_equal(
            actual, expected, decimal=decimal, err_msg=error_message, verbose=verbose
        )

    @classmethod
    def assertMaskedArrayEqual(cls, actual, expected):
        """
        Test that two masked arrays are equal.

        Two checks are performed.  First, the data is checked for equality.
        Secondly, the masks are compared.

        Parameters
        ----------
        actual : np.ma.masked_array
            The first array to compare.
        expected : np.ma.masked_array
            The second array to compare.

        """
        actual = _expand_mask(actual)
        expected = _expand_mask(expected)

        # Check unmasked data:
        cls.assertArrayEqual(actual.data[~actual.mask], expected.data[~expected.mask])
        # Check masks are equal
        cls.assertArrayEqual(actual.mask, expected.mask)

    @classmethod
    def assertMaskedArrayAlmostEqual(cls, actual, expected, decimal=6):
        """
        Test that two masked arrays are equal.

        Two checks are performed.  First, the data is checked for equality.
        Secondly, the masks are compared.

        Parameters
        ----------
        actual : np.ma.masked_array
            The first array to compare.
        expected : np.ma.masked_array
            The second array to compare.
        decimal : int
            Precision for comparison, defaults to 6.

        """
        # Default precision of 6 to match numpy/iris

        actual = _expand_mask(actual)
        expected = _expand_mask(expected)

        # Check unmasked data is equal
        cls.assertArrayAlmostEqual(
            actual.data[~actual.mask],
            expected.data[~expected.mask],
            decimal=decimal,
        )
        # Check masks are equal
        cls.assertArrayEqual(actual.mask, expected.mask)

    def assertCDL(self, netcdf_filename, reference_filename=None):
        """
        Test that the CDL for the given netCDF file matches the contents
        of the reference file.

        If the environment variable ANTS_TEST_CREATE_MISSING is non-empty, the
        reference file is created if it doesn't exist.

        Parameters
        ----------
        netcdf_filename : basestring
            The path to the netCDF file.
        reference_filename : basestring or iterable of basestrings, optional
            The relative path (relative to the test results directory).
            If omitted, the result is generated from the calling method's name,
            class, and module using :meth:`ants.tests.TestCase.result_path`.

        """
        if reference_filename is None:
            reference_path = self.result_path(None, "cdl")
        else:
            reference_path = self.get_result_path(reference_filename)

        flags = ["-h"]
        with subprocess.Popen(
            ["ncdump"] + flags + [netcdf_filename],
            stderr=subprocess.PIPE,
            stdout=subprocess.PIPE,
            universal_newlines=True,
        ) as proc:
            # Ingest the CDL for comparison, excluding first line.
            lines = proc.stdout.readlines()[1:]

        # Sort the dimensions (except for the first, which can be unlimited).
        # This gives consistent CDL across different platforms.
        def sort_key(line):
            return ("UNLIMITED" not in line, line)

        dimension_lines = slice(
            lines.index("dimensions:\n") + 1, lines.index("variables:\n")
        )
        lines[dimension_lines] = sorted(lines[dimension_lines], key=sort_key)
        # Remove NetCDF version attributes from the comparison.
        lines = [line for line in lines if "_NCProperties" not in line]
        cdl = "".join(lines)

        self._check_same(cdl, reference_path, type_comparison_name="CDL")

    def assertCML(self, cubes, reference_filename=None, checksum=True):
        """
        Test that the CML for the given cubes matches the contents of
        the reference file.
        If the environment variable IRIS_TEST_CREATE_MISSING is
        non-empty, the reference file is created if it doesn't exist.
        Args:
        * cubes:
            Either a Cube or a sequence of Cubes.
        Kwargs:
        * reference_filename:
            The relative path (relative to the test results directory).
            If omitted, the result is generated from the calling
            method's name, class, and module using
            :meth:`iris.tests.IrisTest.result_path`.
        * checksum:
            When True, causes the CML to include a checksum for each
            Cube's data. Defaults to True.
        """
        # TODO: Can we remove CML tests?
        # https://code.metoffice.gov.uk/trac/ancil/ticket/1481
        if isinstance(cubes, iris.cube.Cube):
            cubes = [cubes]
        if reference_filename is None:
            reference_filename = self.result_path(None, "cml")

        if isinstance(cubes, (list, tuple)):
            xml = iris.cube.CubeList(cubes).xml(
                checksum=checksum, order=False, byteorder=False
            )
        else:
            xml = cubes.xml(checksum=checksum, order=False, byteorder=False)
        reference_path = self.get_result_path(reference_filename)
        self._check_same(xml, reference_path)

    def _check_reference_file(self, reference_path):
        reference_exists = os.path.isfile(reference_path)
        if not (reference_exists or os.environ.get("ANTS_TEST_CREATE_MISSING")):
            tip = (
                'Set environmental variable "ANTS_TEST_CREATE_MISSING" to '
                "generate the reference file."
            )
            msg = "Missing test result: {}\n{}".format(reference_path, tip)
            raise AssertionError(msg)
        return reference_exists

    def _check_same(self, item, reference_path, type_comparison_name="CML"):
        if self._check_reference_file(reference_path):
            with open(reference_path, "rb") as reference_fh:
                reference = "".join(
                    part.decode("utf-8") for part in reference_fh.readlines()
                )
            self._assert_str_same(reference, item, reference_path, type_comparison_name)
        else:
            self._ensure_folder(reference_path)
            warnings.warn("Creating result file: {}".format(reference_path))
            with open(reference_path, "wb") as reference_fh:
                reference_fh.writelines(part.encode("utf-8") for part in item)

    def _assert_str_same(
        self,
        reference_str,
        test_str,
        reference_filename,
        type_comparison_name="Strings",
    ):
        if reference_str != test_str:
            diff = "".join(
                difflib.unified_diff(
                    reference_str.splitlines(1),
                    test_str.splitlines(1),
                    "Reference",
                    "Test result",
                    "",
                    "",
                    0,
                )
            )
            self.fail(
                "%s do not match: %s\n%s"
                % (type_comparison_name, reference_filename, diff)
            )

    def _ensure_folder(self, path):
        dir_path = os.path.dirname(path)
        if not os.path.exists(dir_path):
            os.makedirs(dir_path)

    @staticmethod
    @contextmanager
    def assertDurationLess(seconds):
        t1 = time.time()
        yield
        t2 = time.time()
        duration = t2 - t1
        if duration > seconds:
            msg = "Test expected to take less than {}s but took {}s"
            raise RuntimeError(msg.format(seconds, duration))

    @staticmethod
    def assertAntsVersionLessEqual(msg, major, minor):
        version = ants.__version__.split("dev")[0].split(".")
        vmaj, vmin = version[:2]
        if int(vmaj) > major or int(vmin) > minor:
            raise RuntimeError(msg)

    @classmethod
    def setUpClass(cls):
        # Ensure that tests are not sensitive to user configuration except
        # for the testing section.
        CONFIG.__init__()
        os.environ["ANTS_NPROCESSES"] = str(1)
        CONFIG.config["ants_decomposition"]["x_split"] = 3
        CONFIG.config["ants_decomposition"]["y_split"] = 3

    @staticmethod
    def get_result_path(relative_path):
        """
        Returns the absolute path to a result file when given the relative path
        as a string, or sequence of strings.
        """
        if not isinstance(relative_path, str):
            relative_path = os.path.join(*relative_path)
        return os.path.abspath(os.path.join(_RESULT_PATH, relative_path))

    def assertAncil(self, actual_filename, relative_reference):
        """
        Ancillary comparison to reference file.

        Compares the provided ancillary with a reference ancillary metadata
        file and where the two do not match, a human readable diff is
        provided.  Raises an assertionError if they are not comparable.
        If the environment variable ANTS_TEST_CREATE_MISSING is non-empty, the
        reference file is created if it doesn't exist.

        Parameters
        ----------
        actual : str
            Filepath to an ancillary fileformat file
        relative_reference : str or sequence of str
            Path of the reference file, relative to the reference results
            directory as a string, or sequence of strings.

        """

        def _get_error_message(actual, expected, reference_filename):
            """
            Formats the error message when an ancil comparison test fails.

            Performs a line by line diff of the mule pumf output for the file
            under test and the reference ancil.

            Parameters
            ----------

            actual : list
                List of mule_pumf output strings for the file under test.
            expected : list
                List of mule_pumf output strings for the reference file.
            reference_filename: str
                Filename for the reference ancil.

            """
            diff_messages = [
                "ANCILs do not match for: {}".format(reference_filename),
            ]

            def _parse_result(lines):
                result = collections.OrderedDict()
                header = None
                for line in lines:
                    if line.startswith("*"):
                        header = line
                        result[header] = list()
                    elif header:
                        result[header].append(line)
                return result

            actual_parts = _parse_result(actual)
            expected_parts = _parse_result(expected)

            # Diffs of content of each section
            common_sections = set(actual_parts) & set(expected_parts)
            for section in common_sections:
                actual_section_length = len(actual_parts[section])
                expected_section_length = len(expected_parts[section])
                actual_entries = actual_parts[section]
                expected_entries = expected_parts[section]

                if actual_section_length != expected_section_length:
                    diff_messages.append(
                        "Section: {} has {} actual entries, but {} entries "
                        "were expected.".format(
                            section,
                            actual_section_length,
                            expected_section_length,
                        )
                    )
                elif actual_entries != expected_entries:
                    diff_messages.append("Section: {}".format(section))
                    content = zip(actual_parts[section], expected_parts[section])
                    for (actual_value, expected_value) in content:
                        # Ignore white space for whether line is different:
                        if str(actual_value).replace(" ", "") != str(
                            expected_value
                        ).replace(" ", ""):
                            # But if line is different, include whitespace in
                            # output for readability:
                            diff_messages.append(
                                "\tActual:   {}\n\tExpected: {}".format(
                                    actual_value, expected_value
                                )
                            )

            # Diff_Messages for sections missing from one or other file
            unexpected_sections = set(actual_parts) - set(expected_parts)
            for section in unexpected_sections:
                diff_messages.append(
                    'Section "{}" was not in reference ancil.'.format(
                        section.lstrip("* ").rstrip(" *")
                    )
                )
            missing_sections = set(expected_parts) - set(actual_parts)
            for section in missing_sections:
                diff_messages.append(
                    'Section "{}" was missing compared to reference '
                    "ancil.".format(section.lstrip("* ").rstrip(" *"))
                )

            error_message = "\n".join(diff_messages)
            return error_message

        if pumf is None:
            raise _PUMF_IMPORT_ERROR

        # Generate pumf from ancillary file provided.
        # pumf can send this to stdout so we capture it using a buffer.
        actual_ancil = mule.AncilFile.from_file(actual_filename)
        actual_buffer = StringIO()
        pumf.pprint(actual_ancil, stdout=actual_buffer)
        actual = actual_buffer.getvalue().split("\n")
        actual_buffer.close()
        # Remove headers and empty lines and undesirable entries.
        actual = actual[actual.index("* fixed_length_header *") :]
        actual = [
            line
            for line in actual
            if not line.startswith("%%")
            and line not in ""
            and "model_version" not in line
        ]

        reference_path = self.get_result_path(relative_reference)
        if not self._check_reference_file(reference_path):
            # Reference file doesn't exist - use actual to create reference
            self._ensure_folder(reference_path)
            warnings.warn("Creating result file: {}".format(reference_path))
            with open(reference_path, "w") as fh:
                fh.write("\n".join(actual))
        else:
            # Load reference file and compare
            with open(reference_path, "r") as fh:
                expected = fh.read().splitlines()
            try:
                # Check for non-whitespace differences only:
                self.assertEqual(
                    str(actual).replace(" ", ""), str(expected).replace(" ", "")
                )
            except AssertionError:
                # Actual != expected, so format the diff in a more helpful way:
                error_message = _get_error_message(actual, expected, reference_path)
                raise AssertionError(error_message)


def _expand_mask(array):
    """
    Return array as masked array with a mask value for each data point.

    Converts a masked array with a mask of a single boolean False to an array
    of False values, or converts a numpy array to a masked array.

    Parameters
    ----------
    array : :class:`np.ma.masked_array`
    Array that may have a single boolean for a mask.

    Returns
    -------
    : :class:`np.ma.masked_array`
    Masked array with a masked value for every data point.

    """
    return np.ma.masked_array(array, mask=(np.zeros(array.shape, dtype=bool)))


def main():
    iris.tests.main()
