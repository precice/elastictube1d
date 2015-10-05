import os, sys

vars = Variables(None, ARGUMENTS)
vars.Add(BoolVariable("petsc", "Enable use of the Petsc linear algebra library.", True))
vars.Add(BoolVariable("boost_inst", "Enable if Boost is available compiled and installed.", False))


env = Environment(variables = vars, ENV = os.environ)
Help(vars.GenerateHelpText(env))

conf = Configure(env)


# ====== precice ======
preciceRoot = os.getenv ('PRECICE_ROOT')
if (preciceRoot == None):
   print 'ERROR: Environment variable PRECICE_ROOT not defined!'
   sys.exit(1)
else:
   print 'Using environment variable PRECICE_ROOT =', preciceRoot

env.Append(CPPPATH = [os.path.join(preciceRoot, 'src')])
env.Append(LIBPATH = [os.path.join(preciceRoot, 'build/last')])
env.Append(CPPDEFINES = ['PRECICE_USE_MPI'])

conf.CheckLib("precice")
conf.CheckLib("python2.7")

# ====== petsc ======
#PETSC_DIR = env["ENV"]["PETSC_DIR"]
#PETSC_ARCH = env["ENV"]["PETSC_ARCH"]

if env["petsc"]:
   env.Append(CPPPATH = [os.path.join( PETSC_DIR, "include"),
                         os.path.join( PETSC_DIR, PETSC_ARCH, "include")])
   env.Append(LIBPATH = [os.path.join( PETSC_DIR, PETSC_ARCH, "lib")])
   conf.CheckLib("petsc")
else:
    env.Append(CPPDEFINES = ['PRECICE_NO_PETSC'])

# ====== lapack ======
conf.CheckLib("lapack")

# cxx = 'mpicxx.mpich2' # For systems offering mpicxx compiler
env["CXX"] = 'mpic++'      # For systems offering mpic++ compiler
# env["CXX"] = 'g++-4.8'

env.Append(CCFLAGS = ["-g3", "-O3", "-Wall", "-std=c++11"])

# ===== boost ======
if env["boost_inst"]:
   conf.CheckLib("boost_system")
   conf.CheckLib("boost_filesystem")

env = conf.Finish()
   
env.Program('StructureSolver', ['structure_solver/structure_solver.cpp'])
env.Program('FluidSolver', ['fluid_solver/fluid_solver.cpp', 'fluid_solver/fluid_nl.cpp'])
