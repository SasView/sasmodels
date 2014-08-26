cylinder.c + cylinder.py is the cylinder model with renamed variables and
sld scaled by 1e6 so the numbers are nicer.  The model name is "cylinder"

cylinder_clone.c + cylinder_clone.py is the cylinder model using the
same interface as the sasview, including calling the model CylinderModel,
so that it can be used as a drop-in replacement for the sasview cylinder
model.

cylindermodel_onefile.py is cylinder.c+cylinder.py merged into one file.
This doesn't actually run yet since sasmodels.gen has not been updated
to support it.  It exists as a proposal.  Note that the function declaration
has been removed since there is enough information in the parameter
definitions to generate it automatically.  Note also that "source" which
used to be all the source has been renamed "includes".

One-file models could coexist with the py+c file models by checking for the
existence of c_blah and creating the appropriate function wrappers.  These
would be appended after any include files.  You shouldn't mix the two forms
within a single model since form_volume needs to be defined before
Iq/Iqxy but after the libraries.

Note: may want to rename form_volume to calc_volume and Iq/Iqxy to
calc_Iq/calc_Iqxy in both the c+py and the one file forms so that the
names are more predictable.  Similarly ER/VR go to calc_ER/calc_VR.

Note: It is possible to translate python code automatically to opencl, using
something like numba, clyther, shedskin or pypy.

Magnetism hasn't been implemented yet.  We may want a separate Imagnetic
calculator with the extra parameters and calculations.   We should
return all desired spin states together so we can share the work of
computing the form factors for the different magnetic contrasts.  This
will mean extending the data handler to support multiple cross sections
in the same data set.
