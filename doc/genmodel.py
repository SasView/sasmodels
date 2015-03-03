import sys
import os
sys.path.insert(0,'..')

# Convert ../sasmodels/models/name.py to sasmodels.models.name
module_name = sys.argv[1][3:-3].replace('/','.').replace('\\','.')
print module_name
module = __import__(module_name)
for part in module_name.split('.')[1:]:
    module = getattr(module, part)
print module

# Load the doc string from the module definition file and store it in rst
from sasmodels import generate
docstr = generate.doc(module)
open(sys.argv[2],'w').write(docstr)
