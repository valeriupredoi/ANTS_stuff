# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.
import abc
import os
import re
import unittest

import ants

LICENSE_TEMPLATE = """# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details."""


class Common(metaclass=abc.ABCMeta):
    def setUp(self):
        # Contruct the full path to all the files in the current directory.
        # Can be run from any subdirectory in the ANTS tree, but only searches
        # the current directory and it's children:
        current_directory = os.getenv("PWD")
        self.all_filepaths = [
            os.path.join(dirpath, filename)
            for dirpath, _, filenames in os.walk(current_directory)
            for filename in filenames
        ]
        self._exclude = [
            os.path.join("KGO"),
            os.path.join("tests", "sources"),
            os.path.join("bin", "ants-launch$"),
            ".*.yml$",
            ".*.cfg$",
            ".*.cfg.*",
            ".+.ini$",
            ".+.org$",
            ".*.pyc$",
        ]
        # Determine the files in the 'ants' package that do not contain
        # any items in the 'exclude' list:
        self._include_files = [
            filepath
            for filepath in self.all_filepaths
            if not self._string_matches_pattern(filepath, self._exclude)
        ]

    @abc.abstractmethod
    def get_files(self):
        pass

    @staticmethod
    def _string_matches_pattern(string, patterns):
        # Return whether the string matches any pattern in the patterns
        # list.
        result = False
        for pattern in patterns:
            if re.search(pattern, string):
                result = True
                break
        return result

    def _get_files(self, include_list):
        """
        Return all files within the ``ants`` package that match a
        pattern in the include list.

        Parameters
        ----------
        include_list : list
            The regexp patterns to use to filter the files.

        Returns
        -------
        : list
            All files within the ``ants`` package that match a pattern
            in the include list.
        """
        output = [
            filepath
            for filepath in self._include_files
            if self._string_matches_pattern(filepath, include_list)
        ]
        return output


class TestLicenseHeaders(Common, unittest.TestCase):
    def get_files(self):
        """
        Return list of names of python files.

        Python files are defined by a '.py' suffix.

        Returns
        -------
        : list of str
        List of filenames for python files.
        """
        return self._get_files(include_list=[r".*\.py$"])

    @staticmethod
    def check_license_header(fnme):
        """Check license header and add where missing."""
        with open(fnme, "r") as fh:
            cont = fh.read()
        license_re = re.compile(r"((\#\!.*|\/\*)\n)?" + re.escape(LICENSE_TEMPLATE))
        match = re.match(license_re, cont)
        messages = []
        if not match:
            if re.search("copyright", cont, flags=re.IGNORECASE):
                msg = "{}: Corrupted copyright/license notice"
                messages.append(msg.format(fnme))
            else:
                msg = "{}: Missing copyright/license notice.  Attempting to add it."
                messages.append(msg.format(fnme))

                if cont.startswith("#!"):
                    insert_pos = cont.index("\n") + 1
                    cont = "{}{}\n{}".format(
                        cont[:insert_pos], LICENSE_TEMPLATE, cont[insert_pos:]
                    )
                else:
                    cont = "{}\n{}".format(LICENSE_TEMPLATE, cont)
                with open(fnme, "w") as fh:
                    fh.write(cont)
        return messages

    def test_license_headers(self):
        files = self.get_files()
        messages = []
        for fnme in files:
            messages.extend(self.check_license_header(fnme))
        msg = "There were license header failures."
        msg = "{}\n{}".format(msg, "\n".join(messages))
        self.assertFalse(messages, msg)


if __name__ == "__main__":
    ants.tests.main()
