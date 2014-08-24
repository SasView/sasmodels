"""
Drop-in replacement for sasview cylinder model.

No rescaling or renaming of the parameters.
"""
from sasmodels.sasview_model import make_class
CylinderModel = make_class('cylinder_clone','CylinderModel')
