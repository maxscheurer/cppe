#!/usr/bin/env python3
"""Setup for cppe

Contributed by M. F. Herbst
"""
import os
import sys
import setuptools

from os.path import join

from setuptools import find_packages, setup
from setuptools.command.test import test as TestCommand

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


# Version of the python bindings and adcc python package.
__version__ = '0.0.8'


#
# Compile and install cppe library
#
def trigger_cppe_build():
    """
    Trigger a build of the adccore library, if it exists in source form.
    """
    if os.path.isfile("build_core.py"):
        abspath = os.path.abspath(".")
        if abspath not in sys.path:
            sys.path.insert(0, abspath)

        import build_core

        build_dir = "build_core"
        install_dir = "pycppe/core"
        build_core.build_install(build_dir, install_dir,
                                 features=["python_iface"])


class BuildDocs(BuildSphinxDoc):
    def run(self):
        this_dir = os.path.dirname(__file__)
        if not os.path.isfile(join(this_dir, "build_core.py")):
            raise SystemExit("Can only build documentation if build_core.py is"
                             "available.")
        else:
            abspath = os.path.abspath(".")
            if abspath not in sys.path:
                sys.path.insert(0, abspath)

        import build_core

        coredoc_dir = join(this_dir, "docs/core")
        build_core.build_documentation(coredoc_dir, latex=False,
                                       html=False, xml=True)
        try:
            import sphinx  # noqa F401

            import breathe  # noqa F401
            import recommonmark  # noqa F401
        except ImportError:
            raise SystemExit("Sphinx or or one of its required plugins not "
                             "found.\nTry 'pip install -U adcc[build_docs]")
        super().run()


#
# Pytest integration
#
class PyTest(TestCommand):
    user_options = [
        ('mode=', 'm', 'Mode for the testsuite (fast or full)'),
        ('skip-update', 's', 'Skip updating testdata'),
        ("pytest-args=", "a", "Arguments to pass to pytest"),
    ]

    def initialize_options(self):
        TestCommand.initialize_options(self)
        self.pytest_args = ""
        self.mode = "fast"
        self.skip_update = False

    def finalize_options(self):
        if self.mode not in ["fast", "full"]:
            raise Exception("Only test modes 'fast' and 'full' are supported")

    def run_tests(self):
        import shlex

        # import here, cause outside the eggs aren't loaded
        import pytest

        if not os.path.isdir("adcc/testdata"):
            raise RuntimeError("Can only test from git repository, "
                               "not from installation tarball.")

        args = ["adcc"]
        args += ["--mode", self.mode]
        if self.skip_update:
            args += ["--skip-update"]
        args += shlex.split(self.pytest_args)
        errno = pytest.main(args)
        sys.exit(errno)


trigger_cppe_build()
setup(
    name='cppe',
    description='CPPE: ',
    long_description="",
    #
    url='https://github.com/maxscheurer/cppe',
    author="Maximilian Scheurer",
    author_email='info@maxscheurer.com',
    license="LGPL v3",
    #
    version=__version__,
    classifiers=[
        'Development Status :: 4 - Beta',
        'License :: OSI Approved :: '
        'GNU Lesser General Public License v3 (LGPLv3)',
        'License :: Free For Educational Use',
        'Intended Audience :: Science/Research',
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Education",
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Operating System :: Unix',
    ],
    #
    packages=find_packages(exclude=["*.test*", "test"]),
    package_data={'cppe': ["pycppe/core/*.so"],
                  '': ["LICENSE*"]},
    zip_safe=False,
    #
    platforms=["Linux", "Mac OS-X", "Unix"],
    python_requires='>=3.5',
    install_requires=[
        'numpy (>= 1.13)',  # Maybe even higher?
    ],
    tests_require=["pytest", "h5py"],
    extras_require={
        "build_docs": ["sphinx>=2", "recommonmark>=0.5.0", "breathe"],
    },
    #
    cmdclass={"build_docs": BuildDocs},
)
