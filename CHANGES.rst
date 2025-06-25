Release notes
=============

v1.0.10 2025-06-25
------------------
* Reorganisation of project structure for wheel support
* Round to nearest integer for star polymer arms
* Fix pickle errors with product models
* Fix random model generator for Hayter MSA
* Doc update: Make equations in core shell sphere docs match the code
* Package dependency: Use tccbox instead of tinycc

v1.0.9 2025-02-28
-----------------
* Fix precision and improve stability at low Q for both the star_polymer and binary_hard_sphere models
* Fix C wrapper generation using proper parameter table when using inline C code
* Remove two unused and redundant functions from the magnetic function library
* Package dependency: NumPy 2+ compatibility
* Package dependency: Fix incompatibility between SciPy 1.15 and pytest
* Package dependency: Increment github action `upload-artifact` version
* Package dependency: Add support for python 3.12
* Doc update: Fixed broken links in README
* Doc update: Clarity on the Polymer Micelle model
* Doc update: Fix typo in equation in Lamellar Stack Paracrystal documentation
* Doc update: Complex number support

v1.0.8 2024-09-26
-----------------
* New model: Bulk ferromagnets model from marketplace
* Doc update: Archive built docs on Github
* Doc update: Display math correctly
* Fix error in FCC paracrystaline models
* Fix parameter name checking in kernel call

v1.0.7 2023-03-23
------------------
* Doc upate: corefunc and optimizer documentation
* Doc update: various models (cylinder, gel_fit, paracrystal, core_shell_ellipsoid)
* Fix error in BCC and FCC paracrystaline models
* Fix rare race condition causing errors
* Fix to allow multiple scattering scripts to run
* Fix problem with models containing complex amplitudes not compiling on the fly
* Restructuring and cross-linking of sasmodels docs
* Update web links and contributor list

v1.0.6 2022-03-17
------------------
* implements generalized 3D description of magnetic SANS
* adds Boucher-type SLD profile to the selection in the model spherical_sld
* improves naming conventions of constants in magnetic SANS
* composite mixture models with multiplicity are now allowed
* improvements to the superball model
* fix to bug in magnetic SLD calculations
* updates to documentation

v1.0.5 2021-01-26
------------------
* fix probelm with P(Q) scaling when using a shell type form factors
* fix problem with wide Gaussian approximation to resolution function at
  very low Q which could sometimes extend into negative Q by truncating the
  resolution contribution at Q = zero
* added support for complex math using openCL
* a number of doc updates

v1.0.4 2020-07-24
------------------
* doc update: how to test a plugin model in 5.x version

v1.0.3 2020-07-06
------------------
* fix polydispersity for python models
* fix dll cache so multiple versions can coexist
* fix caching bug on pyopencl for Intel on OS/X
* fix tinycc compilation bug on windows indicating bad backslash escape
* warn if ER function is present in plugin model but not used
* restore python 2.7 support

v1.0.2 2020-04-24
-----------------
* doc updates: Porod
* doc updates: change "Source intensity" to "Scale factor or Volume fraction"
* doc updates: fix latex errors
* change MODELS to DISTRIBUTIONS in sasmodels.weights
* close open handles when done reading files
* support latest numpy
* improve model tests
* faster generation of model docs
* support compiler flags coming from environment (for .deb packaging)
* build docs as part of continuous integration test suite

v1.0.1 2019-10-03
-----------------
* reparameterize models to preserve constraints within polydisperse systems
* doc update: fix units on lorentzian exponent
* barbell: use consistent equations in docs, code and references
* validation: check triaxial ellipsoid, barbell against MC calculation


v1.0.0 2019-05-20
-----------------
* see git log for details
