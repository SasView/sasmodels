# hook for pyinstaller to find all the extra parts of sasmodels

try:
    from PyInstaller.utils.hooks import collect_data_files

    from sasmodels import data_files

    datas = collect_data_files('sasmodels')

    # module has a convenience function to provide this information
    for target, filenames in data_files():
        for filename in filenames:
            # CRUFT: sasmodels and sasview look in 2 different places for these
            datas.append((filename, target))
            datas.append((filename, target.replace("sasmodels-data", "sasmodels")))

    print(f"sasmodels added datas: {datas}")

    hiddenimports = [
        "pyopencl",
        "sasmodels.compare_many",
        "sasmodels.guyou",
        "sasmodels.jitter",
        "sasmodels.list_pars",
        "sasmodels.multiscat",
        "sasmodels.special",
    ]
    module_collection_mode = "py"

    print(f"sasmodels added hiddenimports: {hiddenimports}")

except ImportError:
    # Runtime or test-time errors for PyInstaller are not important
    pass
