import os, sys
from os.path import join

######### FUNCTIONS ########
def uniqueCheckLib(conf, lib):
   """ Checks for a library and appends it to env if not already appended. """
   if conf.CheckLib(lib, autoadd=0, language="C++"):
      conf.env.AppendUnique(LIBS = [lib])
      return True
   else:
      print("ERROR: Library '" + lib + "' not found!")
      Exit(1)

def checkset_var(varname, default):
    """ Checks if environment variable is set, use default otherwise and print the value. """
    var = os.getenv(varname)
    if not var:
        var = default
        vprint(varname, var)
    else:
        vprint(varname, var, False)
    return var

def print_options(vars):
    """ Print all build option and if they have been modified from their default value. """
    for opt in vars.options:
        try:
            is_default = vars.args[opt.key] == opt.default
        except KeyError:
            is_default = True
        vprint(opt.key, env[opt.key], is_default, opt.help)

def vprint(name, value, default=True, description = None):
    """ Pretty prints an environment variabe with value and modified or not. """
    mod = "(default)" if default else "(modified)"
    desc = "   " + description if description else ""
    print("{0:10} {1:25} = {2!s:8}{3}".format(mod, name, value, desc))

############################

vars = Variables(None, ARGUMENTS)
vars.Add(BoolVariable("parallel", "Compile source-code for parallel version of 1D Example for preCICE", False))
vars.Add(BoolVariable("petsc", "Enable use of the Petsc linear algebra library.", True))
vars.Add(BoolVariable("python", "Enable use of python", False))
vars.Add(PathVariable("libprefix", "Path prefix for libraries", "/usr", PathVariable.PathIsDir))
vars.Add(BoolVariable("supermuc", "Compile tutorial on SuperMUC", False))

env = Environment(variables = vars, ENV = os.environ)
Help(vars.GenerateHelpText(env))

conf = Configure(env)

print
print("Build options ...")
print_options(vars)

prefix = env["libprefix"]

if prefix is not "/usr":  # explicitely add standard search paths
    env.Append(CPPPATH = join( prefix, 'include'))
    env.Append(LIBPATH = join( prefix, 'lib'))

# ====== precice ======

preciceRoot = os.getenv ('PRECICE_ROOT')

if preciceRoot:
    print("PRECICE_ROOT defined, preCICE was probably build from source")
    print('Using environment variable PRECICE_ROOT = ' + preciceRoot)
    env.Append(CPPPATH = [os.path.join(preciceRoot, 'src')])
    env.Append(LIBPATH = [os.path.join(preciceRoot, 'build/last')]) 
    env.Append(CPPDEFINES = ['PRECICE_USE_MPI'])
    uniqueCheckLib(conf, "precice") 
else:
    print("PRECICE_ROOT not defined. Using pkg_config to find libprecice.")
    try:
        uniqueCheckLib(conf, "precice")
    except Exception():
        print("Did you forget to define PRECICE_ROOT?")
        Exit(-1)

# ====== python ======
if env["python"]:
   uniqueCheckLib(conf, "python2.7")
else:
   env.Append(CPPDEFINES = ['PRECICE_NO_PYTHON'])

# ====== petsc ======
if env["petsc"]:
   PETSC_DIR = env["ENV"]["PETSC_DIR"]
   PETSC_ARCH = env["ENV"]["PETSC_ARCH"]
   if not env["mpi"]:
       print("PETSc requires MPI to be enabled.")
       Exit(1)
   env.Append(CPPPATH = [join(prefix, PETSC_DIR, "include"),
                         join(prefix, PETSC_DIR, PETSC_ARCH, "include")])
   env.Append(LIBPATH = [join(prefix, PETSC_DIR, PETSC_ARCH, "lib"),
                         join(prefix, PETSC_DIR, "lib")])
   conf.CheckLib("petsc")
else:
   env.Append(CPPDEFINES = ['PRECICE_NO_PETSC'])

# ======= compiler ======
if env["supermuc"]:
   env["CXX"] = 'mpicc'      # For SuperMUC
else:
   env["CXX"] = 'mpic++'      # For systems offering mpic++ compiler

env.Append(CCFLAGS = ["-g3", "-O0", "-Wall", "-std=c++11"])

# ====== boost ======
uniqueCheckLib(conf, "boost_system")
uniqueCheckLib(conf, "boost_filesystem")
uniqueCheckLib(conf, "boost_thread")
uniqueCheckLib(conf, "boost_log")
uniqueCheckLib(conf, "boost_log_setup")
uniqueCheckLib(conf, "boost_program_options")

# ====== boost ======
uniqueCheckLib(conf, "xml2")
 
# ====== lapack ======
if env["supermuc"]:
   if not os.genenv("MKL_LIB"):
      print('ERROR: Environment variable MKL_LIB not defined! Please load the MKL library for Lapack by executing "module load mkl".')
      sys.exit(1)
   env['GEN_LIB_BUILD_STATIC'] = env["ENV"]["MKL_LIB"]
   env.Append(SHLINKCOM = ' $GEN_LIB_BUILD_STATIC')
   env.Append(LINKCOM = ' $GEN_LIB_BUILD_STATIC')
else:
   if env["parallel"]:
      env.Append(LIBPATH = ['./lib'])
      uniqueCheckLib(conf, "lapack")
   else:
      uniqueCheckLib(conf, "lapack")


env = conf.Finish()
   
if env["parallel"]:
   env.Program('StructureSolver', ['StructureSolver_Parallel/structureDataDisplay.cpp', 'StructureSolver_Parallel/StructureSolver.cpp', 'StructureSolver_Parallel/structureComputeSolution.cpp'])
   env.Program('FluidSolver', ['FluidSolver_Parallel/fluidDataDisplay.cpp', 'FluidSolver_Parallel/FluidSolver.cpp', 'FluidSolver_Parallel/fluidComputeSolution.cpp'])
else:
   env.Program('StructureSolver', ['StructureSolver_Serial/structure_solver.cpp'])
   env.Program('FluidSolver', ['FluidSolver_Serial/fluid_solver.cpp', 'FluidSolver_Serial/fluid_nl.cpp'])   
