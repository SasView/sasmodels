'''
Script to run a batch fit in a series of files and plot the fitted parameters.

Usage syntax:

    python batch_fit.py model.py "sample1.dat, sample2.dat, ..., other_sample.dat"
    (files named sample1.dat, sample2.dat, ..., other_sample.dat)
    
    or if the file names are numbers (and the extension is .dat):

                python batch_fit.py model.py 93190 93210 
                (files named 093190.dat, 093191.dat, ..., 093210.dat)
    
    or for Grasp-like naming:
    
                python batch_fit.py model.py 93190 93210 200
                (files named 093190_200.dat, 093191_201.dat, ..., 093210_202.dat)
    
    The script reads a series of files and fits the model defined by model.py.
    E.g. python batch_fit.py  model_ellipsoid_hayter_msa.py fits a model 
    consisting in an ellipsoid form factor multiplied by a Hayter MSA structure factor.  
    
    The file model.py must load the data using data = load_data('data.txt'), as the script
    replaces 'data.txt' by the files given here.
    
    Modify the call to bumps and the options (minimizer, steps, etc.) as desired.
    
    For each file a directory named Fit_filename is created. There the file fit.par contains
    the fitted parameters.
    
    Finally the fitted parameters are shown for the full series.
    
Example:     

    python batch_fit.py model_ellipsoid_hayter_msa.py 93191 93195 201
    
Note: 

    If sasmodels is not in the path, edit the line sys.path.append to provide the 
    right path to sasmodels.    

'''

from __future__ import print_function
import sys
import os
import numpy as np
import matplotlib.pyplot as plt

''' GET INPUT AND ENSURE MODEL AND DATA FILES ARE DEFINED'''

nargs = len(sys.argv) - 1
if (nargs < 2) or (nargs > 4):
    print ("Error in the list of arguments! \n")
    sys.exit()

modelName = sys.argv[1]
f = open(modelName, 'r')
fileModel = f.read()
f.close()

if nargs == 2:
    dataFiles = sys.argv[2].split(',')
else: 
    numorFirst = int(sys.argv[2])
    numorLast = int(sys.argv[3])
    dataFiles = []
    if nargs == 3:
        for i in range(numorFirst, numorLast+1):
            name = str(i).zfill(6) + '.dat'
            dataFiles.append(name)
    else:
        numorExt = int(sys.argv[4])    
        for i in range(numorFirst, numorLast+1):
            name = str(i).zfill(6) + '_' + str(numorExt) + '.dat'
            numorExt += 1
            dataFiles.append(name)

for file in dataFiles:
    if not os.path.isfile(file.strip()):
        print ("File %s does not exist! \n" % file.strip())
        sys.exit()

        
''' CALL TO BUMPS AND DEFINITION OF FITTING OPTIONS '''

msg0 = "python -m bumps.cli fit.py"
options = " --fit=lm --steps=200 --ftol=1.5e-8 --xtol=1.5e-8 --batch"


''' LOOP OVER FILES AND CALL TO BUMPS FOR EACH OF THEM'''   

for file in dataFiles:
    currentModel = fileModel.replace('data.txt', file.strip())
    f = open('fit.py', 'w')
    f.write(currentModel)
    f.close()
    store = " --store=Fit_" + file.strip()
    msg = msg0 + options + store
    os.system(msg)


''' SHOW FITTED PARAMETERS '''    

parDict = {}
for file in dataFiles:
    parFile = os.path.join('Fit_' + file.strip(), 'fit.par')
    f = open(parFile, 'r')
    lines = f.readlines()
    for line in lines:
        parName = line.split()[0]
        if parName in parDict.keys():
            parList = parDict[parName]
            parList.append(line.split()[1])
            parDict[parName] = parList
        else:
            parDict[parName] = [line.split()[1]]
    
for parName in parDict.keys():
    values = np.array(map(float, parDict[parName]))   
    plt.plot(values)
    plt.title(parName)
    plt.xlabel('Dataset #')
    plt.ylabel('Value')
    plt.show()

