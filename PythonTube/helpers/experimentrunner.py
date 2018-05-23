import subprocess
import multiprocessing
import os


def solve_fluid(config):
    print "STARTING FLUID"
    subprocess.Popen(["python", os.path.join(os.getcwd(),"../FluidSolver.py"), config]).wait()
    print "FINISHING FLUID"


def solve_structure(config):
    print "STARTING STRUCTURE"
    subprocess.Popen(["python", os.path.join(os.getcwd(),"../StructureSolver.py"), config]).wait()
    print "FINISHING STRUCTURE"


configfolder = os.path.join(os.getcwd(),'../configs')

for configname in os.listdir(configfolder):
    print "processing "+configname

    configpath = os.path.join(configfolder, configname)

    fluid = multiprocessing.Process(name='fluid', target=solve_fluid, args=(configpath,))
    structure = multiprocessing.Process(name='structure', target=solve_structure, args=(configpath,))

    fluid.start()
    structure.start()

    structure.join()
    fluid.join()
