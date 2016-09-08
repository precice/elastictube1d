import os, sys

######### FUNCTIONS ########
def uniqueCheckLib(conf, lib):
   """ Checks for a library and appends it to env if not already appended. """
   if conf.CheckLib(lib, autoadd=0, language="C++"):
      conf.env.AppendUnique(LIBS = [lib])
      return True
   else:
      print "ERROR: Library '" + lib + "' not found!"
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
    print "{0:10} {1:25} = {2!s:8}{3}".format(mod, name, value, desc)

############################

vars = Variables(None, ARGUMENTS)
vars.Add(BoolVariable("boost_inst", "Enable if Boost is available compiled and installed.", False))
vars.Add(BoolVariable("parallel", "Compile source-code for parallel version of 1D Example for preCICE", False))
vars.Add(BoolVariable("petsc", "Enable use of the Petsc linear algebra library.", True))
vars.Add(BoolVariable("python", "Enable use of python", False))
vars.Add(BoolVariable("supermuc", "Compile tutorial on SuperMUC", False))

env = Environment(variables = vars, ENV = os.environ)
Help(vars.GenerateHelpText(env))

conf = Configure(env)

print
print "Build options ..."
print_options(vars)

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

uniqueCheckLib(conf, "precice")

# ====== python ======
if env["python"]:
   uniqueCheckLib(conf, "python2.7")
else:
   env.Append(CPPDEFINES = ['PRECICE_NO_PYTHON'])

# ====== petsc ======
if env["petsc"]:
   PETSC_DIR = env["ENV"]["PETSC_DIR"]
   PETSC_ARCH = env["ENV"]["PETSC_ARCH"]
   env.Append(CPPPATH = [os.path.join( PETSC_DIR, "include"),
                         os.path.join( PETSC_DIR, PETSC_ARCH, "include")])
   env.Append(LIBPATH = [os.path.join( PETSC_DIR, PETSC_ARCH, "lib")])
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
if env["supermuc"]:
   boostOnSupermuc = env["ENV"]["BOOST_LIBDIR"]
   env.Append(LIBPATH = [boostOnSupermuc])
uniqueCheckLib(conf, "boost_system")
uniqueCheckLib(conf, "boost_filesystem")
uniqueCheckLib(conf, "boost_thread")
uniqueCheckLib(conf, "boost_log")
uniqueCheckLib(conf, "boost_log_setup")
uniqueCheckLib(conf, "boost_program_options")
   
# ====== lapack ======
if env["supermuc"]:
   if not os.genenv("MKL_LIB"):
      print 'ERROR: Environment variable MKL_LIB not defined! Please load the MKL library for Lapack by executing "module load mkl".'
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
   
env.Program('GSNewtonCoupling', ['Newton_Coupling/GS_coupling.cpp', 'Newton_Coupling/fluid_nl.cpp', 'utils/AD.cpp', 'utils/Linsolve.cpp'])   
