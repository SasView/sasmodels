import sys, os
sys.path.insert(0, os.path.abspath('..'))
from sasmodels import generate, core

# Convert ../sasmodels/models/name.py to name
model_name = os.path.basename(sys.argv[1])[:-3]

# Load the doc string from the module definition file and store it in rst
docstr = generate.make_doc(core.load_model_info(model_name))
open(sys.argv[2],'w').write(docstr)
