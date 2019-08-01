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
import subprocess


def list_committed_files():
    """List all files which are commited in the repo located
    in directory"""
    command = "git ls-tree --full-tree --name-only -r HEAD".split(" ")
    output = subprocess.check_output(command, universal_newlines=True)
    return [line for line in output.split("\n") if line != ""]


def get_repo_root_path():
    return subprocess.check_output("git rev-parse --show-toplevel".split(" "),
                                   universal_newlines=True).strip()


def read_ignore_globs(ignore_file):
    """
    Read a file, given relative to the root of the repo
    and return the list of globs of files to be ignored.

    If the file does not exist, returns an empty list.
    """
    if not os.path.isabs(ignore_file):
        ignore_file = os.path.join(get_repo_root_path(), ignore_file)

    if not os.path.isfile(ignore_file):
        return []
    with open(ignore_file) as ignf:
        return [line.strip() for line in ignf if not line.startswith("#")]
