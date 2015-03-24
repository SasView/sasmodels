from setuptools import setup,find_packages

packages = find_packages(exclude=['contrib', 'docs', 'tests*'])
package_data = {
    'sasmodels.models': ['*.c'],
}
required = []

setup(
    name="sasmodels",
    version = "1.0.0a",
    description = "sasmodels package",
    long_description=open('README.rst').read(),
    author = "SasView Collaboration",
    author_email = "management@sasview.org",
    url = "http://www.sasview.org",
    keywords = "small-angle x-ray and neutron scattering analysis",
    download_url = "https://github.com/SasView/sasmodels",
    classifiers=[
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: Public Domain',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Topic :: Scientific/Engineering :: Physics',
    ],
    packages=packages,
    package_data=package_data,
    install_requires = required,
    extras_require = {
        'OpenCL': ["pyopencl"],
        'Bumps': ["bumps"],
        }

    )