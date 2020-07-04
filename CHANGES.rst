Release notes
=============

v1.0.3 2020-07-06
------------------
* fix dll cache so multiple versions can coexist
* fix polydispersity for python models
* warn if ER function is present in model but not used

v1.0.2 2020-04-14
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
see git log for details
