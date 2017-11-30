import sys

# Fake the existence of 'sas' module
class SasDataloaderLoader:
    pass
class SasDataloader:
    loader = SasDataloaderLoader()
class Sas:
    dataloader = SasDataloader()
sas = Sas()

sys.modules['sas'] = sas
sys.modules['sas.dataloader'] = sas.dataloader
sys.modules['sas.dataloader.loader'] = sas.dataloader.loader

def register(linter):
    pass
