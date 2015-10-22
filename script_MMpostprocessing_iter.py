from __future__ import division
import numpy
import sys
import os

sourceDir = str(sys.argv[1])
R = 0

allFiles = [f for f in os.listdir(sourceDir) if f.startswith('iter')]
N_fine=0
N_coarse=0
tau=0
kappa=0

substrings = allFiles[0].split("_")
tableName = "iterations"
for s in substrings[1:-3]:
    tableName += "_" + s

lhs, rhs = substrings[-3].split("[")
tableName += "_N-" + rhs
    
table_out = open(tableName+".dat","w")
pp_log = open("log.pp", "w")

table = numpy.zeros(shape=(3,3))
tableCoarse = numpy.zeros(shape=(3,3))

Tsteps = numpy.zeros(shape=(3,3))

for l in allFiles:
    filename = os.path.join(sourceDir, l)

    print l
    substrings = l.split("_")
    tau = float(substrings[-2])
    lhs, rhs = substrings[-3].split("[")
    l, r = rhs.split("-")
    N_fine = int(l)
    N_coarse = int(r)
    lhs, rhs = substrings[-1].split("]")
    kappa = int(lhs)
        
    try:
        tsteps, itsFine, itsCoarse, delcols = numpy.loadtxt(filename,skiprows=1,usecols=(0,1, 3, -1), unpack=True)
    except:
        print('cannot load '+str(filename))
        continue

    if (kappa == 1000 and tau == 0.1):
        table[0,0] = itsFine[-1]/tsteps[-1]
        tableCoarse[0,0] = itsCoarse[-1]/tsteps[-1]
        Tsteps[0,0] = tsteps[-1]
    elif(kappa == 100 and tau == 0.1):
        table[0,1] = itsFine[-1]/tsteps[-1]
        tableCoarse[0,1] = itsCoarse[-1]/tsteps[-1]
        Tsteps[0,1] = tsteps[-1]
    elif(kappa == 10 and tau == 0.1):
        table[0,2] = itsFine[-1]/tsteps[-1]
        tableCoarse[0,2] = itsCoarse[-1]/tsteps[-1]
        Tsteps[0,2] = tsteps[-1]
    elif(kappa == 1000 and tau == 0.01):
        table[1,0] = itsFine[-1]/tsteps[-1]
        tableCoarse[1,0] = itsCoarse[-1]/tsteps[-1]
        Tsteps[1,0] = tsteps[-1]
    elif(kappa == 100 and tau == 0.01):
        table[1,1] = itsFine[-1]/tsteps[-1]
        tableCoarse[1,1] = itsCoarse[-1]/tsteps[-1]
        Tsteps[1,1] = tsteps[-1]
    elif(kappa == 10 and tau == 0.01):
        table[1,2] = itsFine[-1]/tsteps[-1]
        tableCoarse[1,2] = itsCoarse[-1]/tsteps[-1]
        Tsteps[1,2] = tsteps[-1]
    elif(kappa == 1000 and tau == 0.001):
        table[2,0] = itsFine[-1]/tsteps[-1]
        tableCoarse[2,0] = itsCoarse[-1]/tsteps[-1]
        Tsteps[2,0] = tsteps[-1]
    elif(kappa == 100 and tau == 0.001):
        table[2,1] = itsFine[-1]/tsteps[-1]
        tableCoarse[2,1] = itsCoarse[-1]/tsteps[-1]
        Tsteps[2,1] = tsteps[-1]
    elif(kappa == 10 and tau == 0.001):
        table[2,2] = itsFine[-1]/tsteps[-1]
        tableCoarse[2,2] = itsCoarse[-1]/tsteps[-1]
        Tsteps[2,2] = tsteps[-1]

        
# write info to terminal
print "Tau" + "   " + "1000fine" + "   " + "1000coarse" + "   " + "100fine" + "   " + "100coarse" + "   " + "10fine" + "   " + "10coarse" + "\n"
print "Table fine iterations: "
print  table
print "Table coarse iterations: "
print  tableCoarse
print "Table timesteps: "
print  Tsteps

# write info to log file
pp_log.write("\n    ### FSI-"+str(N_fine)+":"+str(N_coarse)+" ###\n\n\n")
pp_log.write("Table fine iterations: ")
pp_log.write("\n" + str(table) + "\n\n")
pp_log.write("Table coarse iterations: ")
pp_log.write("\n"+str(tableCoarse)+"\n\n")
pp_log.write("Table timesteps: ")
pp_log.write("\n"+str(Tsteps)+"\n\n")


# write pgfplots table
taus = [0.1, 0.01, 0.001]
table_out.write("# Tau" + "   " + "1000fine" + "   " + "1000coarse" + "   " + "100fine" + "   " + "100coarse" + "   " + "10fine" + "   " + "10coarse" + "\n")
for rowFine, rowCoarse, t in map(None, table, tableCoarse, taus):
    table_out.write(str(t))
    for elemFine, elemCoarse in map(None, rowFine, rowCoarse):
        table_out.write("  " + str(elemFine) + " " + str(elemCoarse))
    table_out.write("\n")

table_out.close()
pp_log.close()
           
