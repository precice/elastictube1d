import numpy as np
from enum import Enum

import configuration_file as conf
from pythonCouplingSolver import solve_1DTube


class Experiments(Enum):
    ImplicitEulerMonolithic = 1  # solve Tube FSI problem with monolithic approach and implicit Euler timestepping
    ImplicitEulerPartitionedImplicit = 2  # solve Tube FSI problem with partitioned approach and implicit Euler timestepping. Implicit coupling is used.
    ImplicitEulerPartitionedExplicit = 7  # solve Tube FSI problem with partitioned approach and implicit Euler timestepping. Explicit coupling is used.
    ImplicitTrapezoidalRuleMonolithic = 3  # solve Tube FSI problem with monolithic approach and trapezoidal rule timestepping
    ImplicitTrapezoidalRulePartitionedImplicit = 4  # solve Tube FSI problem with partitioned approach and trapezoidal rule timestepping. Implicit coupling is used.
    ImplicitTrapezoidalRulePartitionedCustomized = 5  # solve Tube FSI problem with partitioned approach and trapezoidal rule timestepping. Waveform Relaxation coupling is used
    ImplicitTrapezoidalRulePartitionedExplicit = 6  # solve Tube FSI problem with partitioned approach and trapezoidal rule timestepping. Explicit coupling is used

# experimental setup
n_tau = conf.n_tau  # refinement of tau reaching from tau0 * .5 ** 0...n_k
tau0 = conf.tau0  # largest timestep in use

# error for experiments is saved here
e_custom_trapz = n_tau * [None]
e_mono_trapz = n_tau * [None]
e_trapz = n_tau * [None]
e_ie = n_tau * [None]
e_mono_ie = n_tau * [None]
e_ee = n_tau * [None]

# experiments to be performed; if experiments should be skipped, just remove them from this list
experiments = [Experiments.ImplicitTrapezoidalRuleMonolithic,
               Experiments.ImplicitTrapezoidalRulePartitionedCustomized,
               Experiments.ImplicitTrapezoidalRulePartitionedImplicit,
               Experiments.ImplicitEulerPartitionedImplicit,
               Experiments.ImplicitEulerMonolithic]

print "computing experiments"
for k in range(n_tau):
    tau = tau0 * (.5 ** k)
    print "tau = %f" % tau

    if Experiments.ImplicitTrapezoidalRulePartitionedImplicit in experiments:
        print "computing ImplicitTrapezoidalRulePartitionedNormal"
        velocity_trapz, pressure_trapz, crossSectionLength_trapz = solve_1DTube(tau=tau, coupling_mode=conf.CouplingAlgorithm.PartitionedPythonImplicit, time_stepping_scheme=conf.TimeStepping.TrapezoidalRule)

    if Experiments.ImplicitTrapezoidalRulePartitionedCustomized in experiments:
        print "computing ImplicitTrapezoidalRulePartitionedCustomized"
        velocity_custom_trapz, pressure_custom_trapz, crossSectionLength_custom_trapz = solve_1DTube(tau=tau, coupling_mode=conf.CouplingAlgorithm.PartitionedPythonCustomized, time_stepping_scheme=conf.TimeStepping.TrapezoidalRule)

    if Experiments.ImplicitTrapezoidalRuleMonolithic in experiments:
        print "computing ImplicitTrapezoidalRuleMonolithic"
        velocity_mono_trapz, pressure_mono_trapz, crossSectionLength_mono_trapz = solve_1DTube(tau=tau, coupling_mode=conf.CouplingAlgorithm.Monolitic, time_stepping_scheme=conf.TimeStepping.TrapezoidalRule)

    if Experiments.ImplicitEulerPartitionedImplicit in experiments:
        print "computing ImplicitEulerPartitioned"
        velocity_ie, pressure_ie, crossSectionLength_ie = solve_1DTube(tau=tau, coupling_mode=conf.CouplingAlgorithm.PartitionedPythonImplicit, time_stepping_scheme=conf.TimeStepping.ImplicitEuler)

    if Experiments.ImplicitEulerMonolithic in experiments:
        print "computing ImplicitEulerMonolithic"
        velocity_mono_ie, pressure_mono_ie, crossSectionLength_mono_ie = solve_1DTube(tau=tau, coupling_mode=conf.CouplingAlgorithm.Monolitic, time_stepping_scheme=conf.TimeStepping.ImplicitEuler)
