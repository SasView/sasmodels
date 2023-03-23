Release notes
=============

v1.0.7 2023-03-24
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
