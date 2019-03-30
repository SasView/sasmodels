.. fitting_sq.rst

.. Much of the following text was scraped from product.py

.. ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ

.. _Product_Models:

Fitting Models with Structure Factors
-------------------------------------

.. note::

   This help document is under development

.. figure:: p_and_s_buttons.png

**Product models**, or $P@S$ models for short, multiply the structure factor
$S(Q)$ by the form factor $P(Q)$, modulated by the **effective radius** of the
form factor.

Many of the parameters in $P@S$ models take on specific meanings so that they
can be handled correctly inside SasView:

* *scale*:

  In simple $P(Q)$ models **scale** often represents the volume fraction of
  material.
  
  In $P@S$ models **scale** should be set to 1.0, as the $P@S$ model contains a
  **volfraction** parameter.

* *volfraction*:

  The volume fraction of material.

  For hollow shapes, **volfraction** still represents the volume fraction of
  material but the $S(Q)$ calculation needs the volume fraction *enclosed by*
  *the shape.* SasView scales the user-specified volume fraction by the ratio
  form:shell computed from the average form volume and average shell volume
  returned from the $P(Q)$ calculation (the original volfraction is divided
  by the shell volume to compute the number density, and then $P@S$ is scaled
  by that to get the absolute scaling on the final $I(Q)$).

* *radius_effective*:

  The radial distance determining the range of the $S(Q)$ interaction.
  
  This may, or may not, be the same as any "size" parameters describing the
  form of the shape. For example, in a system containing freely-rotating
  cylinders, the volume of space each cylinder requires to tumble will be
  much larger than the volume of the cylinder itself. Thus the effective
  radius will be larger than either the radius or half-length of the
  cylinder. It may be sensible to tie or constrain **radius_effective** to one
  or other of these "size" parameters.

  If just part of the $S(Q)$ calculation, the value of **radius_effective** may
  be polydisperse. If it is calculated by $P(Q)$, then it will be the weighted
  average of the effective radii computed for the polydisperse shape
  parameters.

* *structure_factor_mode*:

  If the $P@S$ model supports the $\beta(Q)$ *decoupling correction* [1] then
  **structure_factor_mode** will appear in the parameter table after the $S(Q)$
  parameters.
  
  If **structure_factor_mode = 0** then the *local monodisperse approximation*
  will be used, i.e.:

    $I(Q)$ = $(scale$ / $volume)$ x $P(Q)$ x $S(Q)$ + $background$

  If **structure_factor_mode = 1** then the $\beta(q)$ correction will be
  used, i.e.:

    $I(Q)$ = $(scale$ x $volfraction$ / $volume)$ x $( <F(Q)^2>$ + $<F(Q)>^2$ x $(S(Q)$ - $1) )$ + $background$

    where $P(Q)$ = $<|F(Q)|^2>$.
    
  This is equivalent to:
  
    $I(Q)$ = $(scale$ / $volume)$ x $P(Q)$ x $( 1$ + $\beta(Q)$ x $(S(Q)$ - $1) )$ + $background$

  The $\beta(Q)$ decoupling approximation has the effect of damping the
  oscillations in the normal (local monodisperse) $S(Q)$. When $\beta(Q)$ = 1
  the local monodisperse approximation is recovered.

  More mode options may appear in future as more complicated operations are
  added.

References
^^^^^^^^^^

.. [#] Kotlarchyk, M.; Chen, S.-H. *J. Chem. Phys.*, 1983, 79, 2461

.. ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ

*Document History*

| 2019-03-30 Paul Kienzle & Steve King
