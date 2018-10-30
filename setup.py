import sys
from setuptools import setup
from setuptools.command.test import test as TestCommand

class PyTest(TestCommand):
    user_options = [('pytest-args=', 'a', "Arguments to pass to pytest")]

    def initialize_options(self):
        TestCommand.initialize_options(self)
        self.pytest_args = ''

    def run_tests(self):
        import shlex
        #import here, cause outside the eggs aren't loaded
        import pytest
        errno = pytest.main(shlex.split(self.pytest_args))
        sys.exit(errno)

def find_version(package):
    """Read package version string from __init__.py"""
    import os
    with open(os.path.join(package, '__init__.py')) as fid:
        for line in fid.readlines():
            if line.startswith('__version__'):
                line = line[:line.find('#')]  # strip comment, if any
                version = line.split('=', 2)[1].strip()
                if version[0] != version[-1] or version[0] not in "\"'":
                    break
                return version[1:-1]
    raise RuntimeError("Could not read version from %s/__init__.py"%package)

install_requires = ['numpy', 'scipy']

if sys.platform=='win32' or sys.platform=='cygwin':
    install_requires.append('tinycc')

setup(
    name='sasmodels',
    version=find_version('sasmodels'),
    description="sasmodels package",
    long_description=open('README.rst').read(),
    author="SasView Collaboration",
    author_email="management@sasview.org",
    url="http://www.sasview.org",
    keywords="small-angle x-ray and neutron scattering analysis",
    download_url="https://github.com/SasView/sasmodels",
    classifiers=[
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: Public Domain',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Programming Language :: C',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Topic :: Scientific/Engineering :: Physics',
    ],
    packages=[
        'sasmodels',
        'sasmodels.models',
        'sasmodels.custom'
    ],
    package_data={
        'sasmodels.models': ['*.c', 'lib/*.c'],
        'sasmodels': ['*.c', '*.cl'],
    },
    install_requires=install_requires,
    extras_require={
        'full': ['docutils', 'bumps', 'matplotlib'],
        'server': ['bumps'],
        'OpenCL': ["pyopencl"],
    },
    build_requires=['setuptools'],
    test_requires=['pytest'],
    cmdclass = {'test': PyTest},
)
