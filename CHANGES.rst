Release notes
=============

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
