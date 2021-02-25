#!/usr/bin/env python3

"""Setup for cppe"""
import os
import sys
import setuptools

try:
    from skbuild import setup
except ImportError:
    print("Please update pip, you need pip 10 or greater,\n"
          " or you need to install the PEP 518 requirements"
          " in pyproject.toml yourself",
        file=sys.stderr,
    )
    raise

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

try:
    import pybind11
except ImportError:
    print("Please install pybind11 before installing cppe.",
          file=sys.stderr,
    )
    raise

__version__ = "0.3.1"


def strip_readme():
    with open("README.md") as fp:
        return "".join([line for line in fp if not line.startswith("<img")])


def is_conda_build():
    return (
        os.environ.get("CONDA_BUILD", None) == "1"
        or os.environ.get("CONDA_EXE", None)
    )


setup(
    name="cppe",
    description="C++ and Python library for Polarizable Embedding",
    long_description=strip_readme(),
    long_description_content_type="text/markdown",
    keywords=[
        "polarizable", "embedding", "excited", "states", "QM/MM",
        "electronic", "structure", "computational", "chemistry", "quantum",
        "spectroscopy", "fast", "multipole", "method",
    ],
    author="Maximilian Scheurer, Peter Reinholdt,"
           " Michael F. Herbst, Lori A. Burns",
    author_email="maximilian.scheurer@iwr.uni-heidelberg.de",
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
    packages=["cppe"],
    package_data={"": ["LICENSE*"]},
    zip_safe=False,
    platforms=["Linux", "Mac OS-X"],
    python_requires=">=3.6",
    setup_requires=["pybind11 >= 2.2"],
    install_requires=["pybind11 >= 2.2"],
    cmake_args=[
        f'-Dpybind11_DIR={pybind11.get_cmake_dir()}',
        '-DENABLE_PYTHON_INTERFACE=ON',
        '-DCMAKE_INSTALL_LIBDIR:PATH=.'
    ],
)
