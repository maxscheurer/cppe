#!/usr/bin/env python3

"""Setup for cppe"""
import os
import sys
import glob
import setuptools

from setuptools import Extension, find_packages, setup
from setuptools.command.test import test as TestCommand
from setuptools.command.build_ext import build_ext as BuildCommand

try:
    from sphinx.setup_command import BuildDoc as BuildSphinxDoc
except ImportError:
    # No sphinx found -> make a dummy class
    class BuildSphinxDoc(setuptools.Command):
        user_options = []

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass


__version__ = "0.2.1"


def strip_readme():
    with open("README.md") as fp:
        return "".join([line for line in fp if not line.startswith("<img")])


#
# Pybind11 BuildExt
#
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


# As of Python 3.6, CCompiler has a `has_flag` method.
# cf http://bugs.python.org/issue26689
def has_flag(compiler, flagname, opts=[]):
    """Return a boolean indicating whether a flag name is supported on
    the specified compiler.
    """
    import tempfile

    with tempfile.NamedTemporaryFile("w", suffix=".cpp") as f:
        f.write("int main (int argc, char **argv) { return 0; }")
        try:
            extra_postargs = ["-Werror", flagname] + opts
            compiler.compile([f.name], extra_postargs=extra_postargs)
        except setuptools.distutils.errors.CompileError:
            return False
    return True


def cpp_flag(compiler, opts=[]):
    """Return the -std=c++[11/14] compiler flag.

    The c++14 is preferred over c++11 (when it is available).
    """
    if has_flag(compiler, "-std=c++14", opts):
        return "-std=c++14"
    elif has_flag(compiler, "-std=c++11", opts):
        return "-std=c++11"
    else:
        raise RuntimeError("Unsupported compiler -- at least C++11 support "
                           "is needed!")


class BuildExt(BuildCommand):
    """A custom build extension for adding compiler-specific options."""
    def build_extensions(self):
        opts = []
        potential_opts = ['-fopenmp']
        potential_linker_args = ['-lgomp']
        if sys.platform == "darwin":
            potential_opts += ["-stdlib=libc++", "-mmacosx-version-min=10.9"]
        if self.compiler.compiler_type == "unix":
            opts.append(cpp_flag(self.compiler, opts))
            potential_opts += ["-fvisibility=hidden", "-Wall", "-Wextra"]
        opts.extend([newopt for newopt in potential_opts
                     if has_flag(self.compiler, newopt, opts)])

        link_args = [newopt for newopt in potential_linker_args
                     if has_flag(self.compiler, newopt, opts)]

        for ext in self.extensions:
            ext.extra_compile_args = opts
            ext.extra_link_args = link_args
        BuildCommand.build_extensions(self)


#
# Pytest integration
#
class PyTest(TestCommand):
    user_options = [
        ("pytest-args=", "a", "Arguments to pass to pytest"),
    ]

    def initialize_options(self):
        TestCommand.initialize_options(self)
        self.pytest_args = ""

    def run_tests(self):
        import shlex

        # import here, cause outside the eggs aren't loaded
        import pytest

        args = []
        args += shlex.split(self.pytest_args)
        errno = pytest.main(args)
        sys.exit(errno)


#
# Main setup code
#
# Setup RPATH on Linux and MacOS
if sys.platform == "darwin":
    extra_link_args = ["-Wl,-rpath,@loader_path",
                       # "-Wl,-rpath,@loader_path/adcc/lib"
                       ]
    runtime_library_dirs = []
elif sys.platform == "linux":
    extra_link_args = []
    runtime_library_dirs = ["$ORIGIN",
                            # "$ORIGIN/adcc/lib"
                            ]
else:
    raise OSError("Unsupported platform: {}".format(sys.platform))

# Setup source directories
sources = glob.glob("cppe/*.cc")
sources += glob.glob("cppe/core/*.cc")
sources += glob.glob("cppe/core/tensors/*.cc")
sources += glob.glob("cppe/utils/*.cc")
sources += glob.glob("cppe/python_iface/*.cc")

# Setup build of the libadcc extension
ext_modules = [
    Extension(
        "cppe", sources=sources,
        include_dirs=[
            # Path to pybind11 headers
            GetPyBindInclude(),
            GetPyBindInclude(user=True),
            "external/eigen3",
        ],
        extra_link_args=extra_link_args,
        runtime_library_dirs=runtime_library_dirs,
        language="c++",
    ),
]

setup(
    name="cppe",
    description="C++ and Python library for Polarizable Embedding",
    long_description=strip_readme(),
    long_description_content_type="text/markdown",
    keywords=[
        "polarizable", "embedding", "excited", "states", "QM/MM",
        "electronic", "structure", "computational", "chemistry", "quantum",
        "spectroscopy",
    ],
    author="Maximilian Scheurer, Peter Reinholdt,"
           " Michael F. Herbst, Lori A. Burns",
    autor_email="maximilian.scheurer@iwr.uni-heidelberg.de",
    license="LGPL v3",
    url="https://github.com/maxscheurer/cppe",
    project_urls={
        "Source": "https://github.com/maxscheurer/cppe",
        "Issues": "https://github.com/maxscheurer/cppe/issues",
    },
    #
    version=__version__,
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "License :: OSI Approved :: GNU Lesser "
        "General Public License v3 (LGPLv3)",
        "License :: Free For Educational Use",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Education",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: POSIX :: Linux",
    ],
    package_data={"": ["LICENSE*"]},
    ext_modules=ext_modules,
    zip_safe=False,
    platforms=["Linux", "Mac OS-X"],
    python_requires=">=3.6",
    setup_requires=["pybind11 >= 2.2"],
    install_requires=["pybind11 >= 2.2"],
    tests_require=[
        "pytest", "numpy", "h5py", "numba", "scipy",
    ],
    # extras_require={
    #     "build_docs": ["sphinx>=2", "breathe", "sphinxcontrib-bibtex",
    #                    "sphinx-automodapi"],
    # },
    cmdclass={"build_ext": BuildExt,
              "test": PyTest,
              # "build_docs": BuildDocs,
              },
)
