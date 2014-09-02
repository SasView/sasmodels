cylinder.c + cylinder.py is the cylinder model with renamed variables and
sld scaled by 1e6 so the numbers are nicer.  The model name is "cylinder"

cylinder_clone.c + cylinder_clone.py is the cylinder model using the
same interface as the sasview, including calling the model CylinderModel,
so that it can be used as a drop-in replacement for the sasview cylinder
model.

lamellar.py is an example of a single file model with embedded C code.

Note: may want to rename form_volume to calc_volume and Iq/Iqxy to
calc_Iq/calc_Iqxy. Similarly ER/VR go to calc_ER/calc_VR.

Note: It is possible to translate python code automatically to opencl, using
something like numba, clyther, shedskin or pypy, so maybe the kernel functions
could be implemented without any C syntax.

Magnetism hasn't been implemented yet.  We may want a separate Imagnetic
calculator with the extra parameters and calculations.   We should
return all desired spin states together so we can share the work of
computing the form factors for the different magnetic contrasts.  This
will mean extending the data handler to support multiple cross sections
in the same data set.

Need to write code to generate the polydispersity loops in python for
kernels that are only implemented in python.  Also, need to implement
an example kernel directly in python.