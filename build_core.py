#!/usr/bin/env python3
"""Advanced CMake wrapper script
for building the project and doxygen Documentation

Contributed by M. F. Herbst
"""
import os
import argparse
import subprocess

from os.path import join


class GetPyBindInclude:
    """Helper class to determine the pybind11 include path

    The purpose of this class is to postpone importing pybind11
    until it is actually installed, so that the ``get_include()``
    method can be invoked. """

    def __init__(self, user=False):
        self.user = user

    def __str__(self):
        import pybind11

        return pybind11.get_include(self.user)


try:
    import multiprocessing

    NCPU = multiprocessing.cpu_count()
except ImportError:
    NCPU = 1


def configure(build_dir, source_dir, install_dir, build_type=None, features=[]):
    args = ["cmake", "-DCMAKE_INSTALL_PREFIX=" + install_dir]

    if build_type in ["Release", "Debug", "MinSizeRel", "RelWithDebInfo"]:
        args += ["-DCMAKE_BUILD_TYPE=" + build_type]
    elif build_type == "SanitizeAddress":
        cpflags = "-O1 -g -fsanitize=address -fno-omit-frame-pointer"
        ldflags = "-fsanitize=address"
        args += ["-DCMAKE_CXX_FLAGS_DEBUG={}".format(cpflags),
                 "-DCMAKE_C_FLAGS_DEBUG={}".format(cpflags),
                 "-DCMAKE_EXE_LINKER_FLAGS_DEBUG=\"{}\"".format(ldflags),
                 "-DCMAKE_BUILD_TYPE=Debug",
                 ]
    elif build_type:
        raise SystemExit("Unknown build type: " + build_type)

    if "python_iface" in features:
        args += ["-DENABLE_PYTHON_INTERFACE=ON"]
        args += ["-DINSTALL_DEVEL_HEADERS=OFF"]
        args += ["-DPYMOD_INSTALL_LIBDIR={}".format(install_dir)]
    subprocess.check_call(args + [source_dir], cwd=build_dir)


def install(build_dir, n_jobs=NCPU, verbose=False):
    args = ["make", "-j" + str(n_jobs), "install"]
    if verbose:
        args += ["VERBOSE=1"]
    subprocess.check_call(args, cwd=build_dir)


def build_documentation(docout_dir, latex=True, html=True, xml=False):
    # Build documentation into the docout_dir directory
    docout_dir = os.path.abspath(docout_dir)
    os.makedirs(docout_dir, exist_ok=True)

    try:
        subprocess.check_output(["doxygen", "-v"])
    except (FileNotFoundError, subprocess.CalledProcessError):
        raise OSError("Doxygen not installed on the system. "
                      "Please install 'doxygen'")

    this_dir = os.path.abspath(os.path.dirname(__file__))
    doxydata = [
        "PROJECT_NAME     = CPPE",
        "PROJECT_NUMBER   = ",  # TODO input git tag
        "PROJECT_BRIEF    = \"CPPE core C++ library\"",
        #
        "INPUT            = {}".format(join(this_dir, "src")),
        "RECURSIVE        = YES",
        "OUTPUT_DIRECTORY = {}".format(docout_dir),
        # "EXCLUDE_PATTERNS = */tests/* */adcman/* */legacy/* *.cc",
        #
        "GENERATE_LATEX   = {}".format("YES" if latex else "NO"),
        "GENERATE_HTML    = {}".format("YES" if html else "NO"),
        "GENERATE_XML     = {}".format("YES" if xml else "NO"),
        #
        "LOOKUP_CACHE_SIZE = 3",
    ]
    doxyfile = join(docout_dir, "Doxyfile")
    if not os.path.isfile(join(docout_dir, "Doxyfile")):
        with open(doxyfile, "w") as fp:
            fp.write("\n".join(doxydata))

    subprocess.check_call("doxygen", cwd=docout_dir)


def build_install(build_dir, install_dir, n_jobs=NCPU, build_type=None,
                  verbose=False, features=[]):
    install_dir = os.path.abspath(install_dir)
    build_dir = os.path.abspath(build_dir)
    source_dir = os.path.abspath(os.path.dirname(__file__))

    if os.path.isdir(build_dir):
        if not os.path.isfile(join(build_dir, "Makefile")):
            raise SystemExit("Something went wrong setting up cmake\n"
                             "Please delete " + build_dir + " and try again.")
        # TODO: remove
        configure(build_dir, source_dir, install_dir, build_type=build_type,
                  features=features)
        install(build_dir, n_jobs=n_jobs, verbose=verbose)
    else:
        os.mkdir(build_dir)
        configure(build_dir, source_dir, install_dir, build_type=build_type,
                  features=features)
        install(build_dir, n_jobs=n_jobs, verbose=verbose)


def main():
    parser = argparse.ArgumentParser(
        description="Simple wrapper script around CMake to configure common "
        "build modes."
    )

    this_dir = os.path.abspath(os.path.dirname(__file__))
    # if not os.path.isfile(join(parent_dir, "adccore/build_adccore.py")) or \
    #    not os.path.isfile(join(parent_dir, "adcc/__init__.py")) or \
    #    not os.path.isfile(join(parent_dir, "extension/ExportAdcc.cc")):
    #     raise SystemExit("Could not find install location for adccore "
    #                      "binaries. The script expects the adccore files to "
    #                      "be located inside the files of adcc. See the adcc "
    #                      "website for further instructions.")
    install_dir = join(this_dir, "pycppe")

    parser.add_argument("--verbose", "-v", default=False, action="store_true",
                        help="Run make in verbose mode.")
    parser.add_argument("--directory", "-d",  default="build", metavar="DIR",
                        help="The directory in which files are built.")
    parser.add_argument("--jobs", "-j", default=NCPU, metavar="N",
                        help="Number of jobs to use during build.")
    parser.add_argument("--type", "-t", default=None, metavar="BUILD_TYPE",
                        choices=["Release", "Debug", "SanitizeAddress",
                                 "MinSizeRel", "RelWithDebInfo"],
                        help="The build type to configure.")
    parser.add_argument("--features", default=[],
                        nargs="+", help="Select optional features for build.",
                        choices=["python_iface"])
    parser.add_argument("--documentation", default=False, action="store_true",
                        help="Build documentation using doxygen.")

    args = parser.parse_args()
    build_install(args.directory, install_dir, n_jobs=args.jobs,
                  build_type=args.type, verbose=args.verbose,
                  features=args.features)

    if args.documentation:
        build_documentation(join(args.directory, "docs"))


if __name__ == "__main__":
    main()
