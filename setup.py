try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

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

setup(
    name='sasmodels',
    version=find_version('sasmodels'),
    description="sasmodels package",
    long_description=open('README.md').read(),
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
        'sasmodels': ['*.c', '*.cl', 'convert.json'],
    },
    install_requires = [
    ],
    extras_require = {
        'OpenCL': ["pyopencl"],
        'Bumps': ["bumps"],
        'TinyCC': ["tinycc"],
    },
)
