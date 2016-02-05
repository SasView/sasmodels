**********
SAS Models
**********

Small angle X-ray and Neutron (SAXS and SANS) scattering examines the
scattering patterns produced by a beam travelling through the sample
and scattering at low angles.  The scattering is computed as a function
of $q_x$ and $q_y$, which for a given beam wavelength corresponds to
particular scattering angles. Each pixel on the detector corresponds to
a different scattering angle. If the sample is unoriented, the scattering
pattern will appear as rings on the detector.  In this case, a circular
average can be taken with 1-dimension data at $q = \surd (q_x^2 + q_y^2)$
compared to the orientationally averaged SAS scattering pattern.

Models have certain features in common.

Every model has a *scale* and a *background*.

Talk about orientation, with diagrams for orientation so that we don't need
a link on every model page?

.. _orientation:

.. figure: img/orientation1.jpg

    Orientation in 3D

.. figure: img/orientation2.jpg

    Orientation cross sections

Talk about polydispersity.

Talk about magnetism, converting the magnetism help file to inline text here,
with links so that models can point back to it.

Need to talk about structure factors even though we don't have any
implemented yet.