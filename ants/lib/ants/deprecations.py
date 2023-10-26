# (C) Crown Copyright, Met Office. All rights reserved.
#
# This file is part of ANTS and is released under the BSD 3-Clause license.
# See LICENSE.txt in the root of the repository for full licensing details.
import warnings


def issue_save_deprecation(interface_name):
    msg = save_deprecation_message(interface_name)
    warnings.warn(msg, FutureWarning, stacklevel=2)


def save_deprecation_message(interface_name):
    return (
        f"'{interface_name}' is deprecated as of ANTS 1.1.0 and will be "
        "removed in a future release. Please use the appropriate function "
        "from the 'ants.io.save' module."
    )


def save_deprecation_message_for_docstring(interface_name):
    msg = save_deprecation_message(interface_name).replace("'", "``")
    formatted_msg = f"{8*' '}.. deprecated:: 1.1.0\n{11*' '}{msg}\n\n"
    return formatted_msg
