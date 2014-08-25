"""
Drop-in replacement for sasview cylinder model.

No rescaling or renaming of the parameters.
"""
from sasmodels.sasview_model import make_class
from sasmodels.models import cylinder_clone as cylinder
CylinderModel = make_class(cylinder, dtype='single')
