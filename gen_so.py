# Script to compile all C/C++ models
# Part of the build process, referenced from setup.py
import sys, os

def main():
    sasmodels = os.path.dirname(os.path.realpath(__file__))
    sys.path.insert(0, sasmodels)
    from sasmodels import generate, core

    # Convert ../sasmodels/models/name.py to name
    model_loc = os.path.join(sasmodels, "sasmodels", "models") 
    for i in os.listdir(model_loc):
        # Choose only relevant python files
        if i.endswith(".py") and not i.startswith("_"):
            model_name = os.path.basename(i)[:-3]
            model_info = core.load_model_info(model_name)
            # Run the conversion but don't delete the .so
            model = core.build_model(model_info)

if __name__ == "__main__":
    main()



