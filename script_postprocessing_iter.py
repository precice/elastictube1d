from __future__ import division
import numpy
import sys
import os

sourceDir = str(sys.argv[1])
R = 0

allFiles = [f for f in os.listdir(sourceDir) if f.startswith('iter')]
N=0
tau=0
kappa=0

substrings = allFiles[0].split("_")
tableName = "iterations"
for s in substrings[1:-4]:
    tableName += "_" + s

table_out = open(tableName+".dat","w")

table = numpy.zeros(shape=(3,3))

for l in allFiles:
    filename = os.path.join(sourceDir, l)

    substrings = l.split("_")
    tau = float(substrings[-2])
    lhs, rhs = substrings[-3].split("[")
    N = int(rhs)
    lhs, rhs = substrings[-1].split("]")
    kappa = int(lhs)
        
    #try:
    tsteps, its, delcols = numpy.loadtxt(filename,skiprows=1,usecols=(0,1,-1), unpack=True)
    #except:
    #    print('cannot load '+str(filename))
    #    continue

    if (kappa == 1000 and tau == 0.1):
        table[0,0] = its[-1]/tsteps[-1]
    elif(kappa == 100 and tau == 0.1):
        table[0,1] = its[-1]/tsteps[-1]
    elif(kappa == 10 and tau == 0.1):
        table[0,2] = its[-1]/tsteps[-1]
    elif(kappa == 1000 and tau == 0.01):
        table[1,0] = its[-1]/tsteps[-1]
    elif(kappa == 100 and tau == 0.01):
        table[1,1] = its[-1]/tsteps[-1]
    elif(kappa == 10 and tau == 0.01):
        table[1,2] = its[-1]/tsteps[-1]
    elif(kappa == 1000 and tau == 0.001):
        table[2,0] = its[-1]/tsteps[-1]
    elif(kappa == 100 and tau == 0.001):
        table[2,1] = its[-1]/tsteps[-1]
    elif(kappa == 10 and tau == 0.001):
        table[2,2] = its[-1]/tsteps[-1]            

        
table_out.write("Tau" + "   " + "1000" + "   " + "100" + "   " + "10" + "\n")
print "Tau" + "   " + "1000" + "   " + "100" + "   " + "10" + "\n"

taus = [0.1, 0.01, 0.001]
for row, t in map(None, table, taus):
    table_out.write(str(t))
    for elem in row:
        table_out.write("  " + str(elem))
    table_out.write("\n")

table_out.close()
           
