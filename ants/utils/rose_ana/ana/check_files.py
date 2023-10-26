#!/usr/bin/env python
# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.

import os
import os.path
import re
import subprocess as subp

from rose.apps.rose_ana import AnalysisTask
from rose.env import env_var_process


class FilesTest(AnalysisTask):
    """
    Base class for file comparison tests

    This class cannot be used on its own. A subclass must be created with a
    specified `_py_file` which denotes the python module being used for the
    file comparison. The attribute `_py_cmd`, and the methods `create_cmd` and
    `create_args`, can be re-defined if required for the given `_py_file`.
    """

    _py_file = None
    _py_cmd = "{0} {1} {2}"

    def run_analysis(self):
        """Loop over sources and targets and generate comparisons"""
        config = self.process_config()
        (sources, targets) = self.process_opt_files(config)

        self.passed = True
        for (source, target) in zip(sources, targets):
            if os.path.isfile(source) and os.path.isfile(target):
                args = self.create_args(source, target)
                cmd = self.create_cmd(args)
                self.reporter("Command to run: {0}".format(cmd))
                ret = subp.call(cmd, shell=True)
                self.reporter("Return Code: {0:d}".format(ret))
                if ret > 0:
                    self.passed = False
            else:
                msg = "{0} {1} does not exist"
                if not os.path.isfile(source):
                    self.reporter(msg.format("source", source))
                if not os.path.isfile(target):
                    self.reporter(msg.format("target", target))
                self.passed = False

    def process_config(self):
        """The config input from rose_ana has not been processed in the
        same way as the options, so doing that here."""
        config = {}
        for key, value in self.config.items():
            # If the value contains newlines, split it into a list
            # and either way remove any quotation marks and process
            # any environment variables
            value = env_var_process(value)
            values = value.split("\n")
            for ival, value in enumerate(values):
                values[ival] = re.sub(r"^((?:'|\")*)(.*)(\1)$", r"\2", value)
            if len(values) == 1:
                values = values[0]
            config[key] = values
        return config

    def process_opt_files(self, config):
        """Process the script options"""

        # Get config items into options if they don't already exist
        for (key, value) in config.items():
            if key not in self.options:
                self.options[key] = os.path.expandvars(value)

        # Get configuration options
        source_root = self.options.pop("source_root", None)
        target_root = self.options.pop("target_root", None)
        source_files = self.options.pop("source_file", None)
        target_files = self.options.pop("target_file", None)

        # Check consistency of options
        assert source_files, "source_file must be specified"
        if not target_files:
            msg = "target_root must be specified if target_file is not"
            assert target_root, msg
        if source_root:
            msg = "There can be only one source root"
            assert isinstance(source_root, str), msg
        if target_root:
            msg = "There can be only one target root"
            assert isinstance(target_root, str), msg

        # If only one entry in file options then make it a single element list
        if isinstance(source_files, str):
            source_files = [source_files]
        if isinstance(target_files, str):
            target_files = [target_files]

        # Deal with source files, new files to compare
        if source_root:
            for i, source_file in enumerate(source_files):
                if not source_file.startswith("/"):
                    source_files[i] = os.path.join(source_root, source_file)

        # Deal with target files, old files to compare against
        if not target_files:
            target_files = [os.path.basename(x) for x in source_files]
        if target_root:
            for i, target_file in enumerate(target_files):
                if not target_file.startswith("/"):
                    target_files[i] = os.path.join(target_root, target_file)

        # Return lists of source and target files
        msg = "There must be equal numbers of source and target files"
        assert len(target_files) == len(source_files), msg
        return source_files, target_files

    def create_cmd(self, args):
        """Create command that is going to be run for comparison"""
        return self._py_cmd.format(self.options["python_setup"], self._py_file, args)


class CumfTest(FilesTest):
    """
    Subclass comparing fields files using mule-cumf
    """

    _py_file = "check_cumf.py"

    def create_args(self, source, target):
        """Create argument list for command based on files to compare"""
        return "{0} {1}".format(source, target)


class NccmpTest(FilesTest):
    """
    Subclass comparing NetCDF files using nccmp
    """

    _py_file = "check_nccmp.py"

    def create_args(self, source, target):
        return "{0} {1} --exclude-attr history,_NCProperties,Conventions".format(
            source, target
        )
