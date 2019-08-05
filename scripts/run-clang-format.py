#!/usr/bin/env python3
## vi: tabstop=4 shiftwidth=4 softtabstop=4 expandtab
## ---------------------------------------------------------------------
##
## Copyright 2019 Michael F. Herbst
##
## Licensed under the Apache License, Version 2.0 (the "License");
## you may not use this file except in compliance with the License.
## You may obtain a copy of the License at
##
##     http://www.apache.org/licenses/LICENSE-2.0
##
## Unless required by applicable law or agreed to in writing, software
## distributed under the License is distributed on an "AS IS" BASIS,
## WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
## See the License for the specific language governing permissions and
## limitations under the License.
##
## ---------------------------------------------------------------------
import os
import glob
import common
import argparse
import subprocess


def read_ignored_files():
    ignore_file = os.path.join(common.get_repo_root_path(), ".clang-format.ignore")
    if not os.path.isfile(ignore_file):
        return []
    with open(ignore_file) as ignf:
        return [line.strip() for line in ignf if not line.startswith("#")]


def main():
    parser = argparse.ArgumentParser(
        description="Run clang-format on all currently committed c++ files"
    )
    parser.add_argument("-clang-format", metavar="PATH", default="clang-format",
                        help="Name or path of clang-format executable to use.")
    args = parser.parse_args()

    cpp_extensions = [".cpp", ".hpp", ".cxx", ".hxx", ".hh", ".cc", ".h", ".c"]
    root = common.get_repo_root_path()
    cpp_files = [
        os.path.join(root, file) for file in common.list_committed_files()
        if os.path.splitext(file)[1] in cpp_extensions
    ]
    ignore_globs = common.read_ignore_globs(".clang-format.ignore")
    cpp_files = [
        os.path.join(root, name) for name in common.list_committed_files()
        if os.path.splitext(name)[1] in cpp_extensions and
        not any(glob.fnmatch.fnmatch(name, ignore) for ignore in ignore_globs)
    ]

    commandline = [args.clang_format, "-style=file", "-i"]
    for cfile in cpp_files:
        subprocess.check_call(commandline + [cfile])


if __name__ == "__main__":
    main()
