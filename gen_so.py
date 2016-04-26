# Script to compile all C/C++ models
# Part of the build process, referenced from setup.py
import sys, os

def main():
    sasmodels = os.path.dirname(os.path.realpath(__file__))
    sys.path.insert(0, sasmodels)
    from sasmodels import generate, core

    # Convert ../sasmodels/models/name.py to name
    for model_name in core.list_models():
        # Choose only relevant python files
            # model_info = core.load_model_info(model_name)
            # Run the conversion but don't delete the .so
            model = core.precompile_dll(model_name)

if __name__ == "__main__":
    main()



