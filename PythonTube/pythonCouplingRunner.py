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
experiments = [Experiments.ImplicitTrapezoidalRuleMonolithic, Experiments.ImplicitTrapezoidalRulePartitionedCustomized, Experiments.ImplicitTrapezoidalRulePartitionedImplicit, Experiments.ImplicitEulerPartitionedImplicit]

print "computing reference solution"
velocity_ref, pressure_ref, crossSectionLength_ref = solve_1DTube(tau=tau0 * (.5 ** (n_tau)), coupling_mode=conf.CouplingAlgorithm.Monolitic, time_stepping_scheme=conf.TimeStepping.TrapezoidalRule)
print "done."

print "computing experiments"
for k in range(n_tau):
    tau = tau0 * (.5 ** k)
    print "tau = %f" % tau

    if Experiments.ImplicitTrapezoidalRulePartitionedImplicit in experiments:
        print "computing ImplicitTrapezoidalRulePartitionedNormal"
        velocity_trapz, pressure_trapz, crossSectionLength_trapz = solve_1DTube(tau=tau, coupling_mode=conf.CouplingAlgorithm.PartitionedPythonImplicit, time_stepping_scheme=conf.TimeStepping.TrapezoidalRule)
        e_trapz[k] = np.sum(np.abs(velocity_trapz - velocity_ref)+np.abs(pressure_trapz - pressure_ref))

    if Experiments.ImplicitTrapezoidalRulePartitionedCustomized in experiments:
        print "computing ImplicitTrapezoidalRulePartitionedCustomized"
        velocity_custom_trapz, pressure_custom_trapz, crossSectionLength_custom_trapz = solve_1DTube(tau=tau, coupling_mode=conf.CouplingAlgorithm.PartitionedPythonCustomized, time_stepping_scheme=conf.TimeStepping.TrapezoidalRule)
        e_custom_trapz[k] = np.sum(np.abs(velocity_custom_trapz - velocity_ref)+np.abs(pressure_custom_trapz - pressure_ref))

    if Experiments.ImplicitTrapezoidalRuleMonolithic in experiments:
        print "computing ImplicitTrapezoidalRuleMonolithic"
        velocity_mono_trapz, pressure_mono_trapz, crossSectionLength_mono_trapz = solve_1DTube(tau=tau, coupling_mode=conf.CouplingAlgorithm.Monolitic, time_stepping_scheme=conf.TimeStepping.TrapezoidalRule)
        e_mono_trapz[k] = np.sum(np.abs(velocity_mono_trapz - velocity_ref)+np.abs(pressure_mono_trapz - pressure_ref))

    if Experiments.ImplicitEulerPartitionedImplicit in experiments:
        print "computing ImplicitEulerPartitioned"
        velocity_ie, pressure_ie, crossSectionLength_ie = solve_1DTube(tau=tau, coupling_mode=conf.CouplingAlgorithm.PartitionedPythonImplicit, time_stepping_scheme=conf.TimeStepping.ImplicitEuler)
        e_ie[k] = np.sum(np.abs(velocity_ie - velocity_ref)+np.abs(pressure_ie - pressure_ref))

    if Experiments.ImplicitEulerMonolithic in experiments:
        print "computing ImplicitEulerMonolithic"
        velocity_mono_ie, pressure_mono_ie, crossSectionLength_mono_ie = solve_1DTube(tau=tau, coupling_mode=conf.CouplingAlgorithm.Monolitic, time_stepping_scheme=conf.TimeStepping.ImplicitEuler)
        e_mono_ie[k] = np.sum(np.abs(velocity_mono_ie - velocity_ref)+np.abs(pressure_mono_ie - pressure_ref))

# create a loglog plot of error for different schemes and different timestep width
import matplotlib.pyplot as plt

plt.figure(2)
plt.xlabel("tau")
plt.ylabel("error")

thehandles = []
thelabels = []

if Experiments.ImplicitEulerPartitionedImplicit in experiments:
    p, = plt.loglog(tau0 * .5 ** np.arange(n_tau), e_ie, '+')
    thehandles.append(p)
    thelabels.append(Experiments.ImplicitEulerPartitionedImplicit.name)
    print Experiments.ImplicitEulerPartitionedImplicit.name
    print e_ie

if Experiments.ImplicitEulerMonolithic in experiments:
    p, = plt.loglog(tau0 * .5 ** np.arange(n_tau), e_mono_ie, '^')
    thehandles.append(p)
    thelabels.append(Experiments.ImplicitEulerMonolithic.name)
    print Experiments.ImplicitEulerMonolithic.name
    print e_mono_ie

if Experiments.ImplicitTrapezoidalRulePartitionedImplicit in experiments:
    p, = plt.loglog(tau0 * .5 ** np.arange(n_tau), e_trapz, 'x')
    thehandles.append(p)
    thelabels.append(Experiments.ImplicitTrapezoidalRulePartitionedImplicit.name)
    print Experiments.ImplicitTrapezoidalRulePartitionedImplicit.name
    print e_trapz

if Experiments.ImplicitTrapezoidalRulePartitionedCustomized in experiments:
    p, = plt.loglog(tau0 * .5 ** np.arange(n_tau), e_custom_trapz, '*')
    thehandles.append(p)
    thelabels.append(Experiments.ImplicitTrapezoidalRulePartitionedCustomized.name)
    print Experiments.ImplicitTrapezoidalRulePartitionedCustomized.name
    print e_custom_trapz

if Experiments.ImplicitTrapezoidalRuleMonolithic in experiments:
    p, = plt.loglog(tau0 * .5 ** np.arange(n_tau), e_mono_trapz, '.')
    thehandles.append(p)
    thelabels.append(Experiments.ImplicitTrapezoidalRuleMonolithic.name)
    print Experiments.ImplicitTrapezoidalRuleMonolithic.name
    print e_mono_trapz

plt.grid('on')
plt.legend(handles=thehandles, labels=thelabels)
plt.show()
