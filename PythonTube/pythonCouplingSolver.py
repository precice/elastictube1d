# based on https://github.com/precice/elastictube1d and
# [1] J. Degroote, P. Bruggeman, R. Haelterman, and J. Vierendeels. Stability of a coupling technique for partitioned solvers in FSI applications. Computers & Structures, 2008.
# for time integration details see
# [2] Gresho, P. M., & Sani, R. L. (2000). Incompressible Flow and the Finite Element Method, Isothermal Laminar Flow. John Wiley & Sons. Retrieved from http://books.google.de/books?id=m_tQAAAAMAAJ

from __future__ import division

import matplotlib
import numpy as np

import configuration_file as config
from pythonCouplingHelpers.fixedPointIteration import IQNILSScheme
from pythonCouplingHelpers.solid import solve_solid
from thetaScheme import perform_monolithic_implicit_trapezoidal_rule_step, perform_partitioned_implicit_trapezoidal_rule_step, perform_partitioned_implicit_euler_step, perform_monolithic_implicit_euler_step

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.animation as manimation

import datetime

from output import create_output, create_video


def is_partitioned_approach(coupling_mode):
    """
    this is a convenience method for classification of concrete coupling algorithms into partitioned or monolithic approach
    :param coupling_mode:
    :return:
    """
    if (coupling_mode is config.CouplingAlgorithm.PartitionedPythonExplicit) or \
                (coupling_mode is config.CouplingAlgorithm.PartitionedPythonImplicit) or \
                (coupling_mode is config.CouplingAlgorithm.PartitionedPythonCustomized):
        return True
    elif (coupling_mode is config.CouplingAlgorithm.Monolitic):
        return False
    else:
        raise Exception("cannot classify coupling mode [%s]", (coupling_mode.name))


def performs_iterations(coupling_mode):
    """
    this is a convenience method for classification of concrete coupling algorithms into algorithms with and without iteration
    :param coupling_mode:
    :return:
    """
    if (coupling_mode is config.CouplingAlgorithm.PartitionedPythonImplicit) or \
        (coupling_mode is config.CouplingAlgorithm.PartitionedPythonCustomized):
        return True
    elif coupling_mode is config.CouplingAlgorithm.PartitionedPythonExplicit:
        return False
    else:
        raise Exception("cannot classify coupling mode [%s]", (coupling_mode.name))


def solve_1DTube(N=config.n_elem, tau=config.tau0, T_max=config.T_max, L=config.L, velocity_in=config.velocity_in, coupling_mode=config.CouplingAlgorithm.Monolitic, time_stepping_scheme=config.time_stepping_scheme):

    sim_start_time = datetime.datetime.now()
    plotting_mode = config.PlottingModes.OFF
    output_mode = config.OutputModes.NETCDF


    dx = L/N  # element length

    fixed_point_solver = IQNILSScheme(config.underrelaxation_factor)

    # initialize unknowns
    N_steps = int(T_max / tau)

    if plotting_mode == config.PlottingModes.DEBUG:
        fig, ax = plt.subplots(2)
    elif plotting_mode == config.PlottingModes.VIDEO:
        fig, ax = plt.subplots(1)
        FFMpegWriter = manimation.writers['ffmpeg']
        metadata = dict(title='1DElasticTube', artist='BenjaminRueth',
            comment='ForPhD')
        writer = FFMpegWriter(fps=N_steps / T_max, metadata=metadata)
        writer.setup(fig, "writer_test_" + coupling_mode.name + ".mp4", 100)

    # here we perform a few iterations towards a quasi-stationary case to obtain proper initial conditions. A stationary
    # case solution (i.e. using velocity_in(0)) does not suffice, but we have to use a linearized version.
    # See [2],p.739,eq.3.16-34, not only the continuity equation w.r.t to the velocity, but also w.r.t the first
    # derivative of the velocity has to be fulfilled.


    pressure0 = config.p0 * np.ones(N + 1)  # initialize with equilibrium pressure
    crossSection0 = config.crossSection0(N)
    velocity0 = config.velocity_in(0) * crossSection0[0] * np.ones(N + 1) / crossSection0  # initialize such that mass conservation is fulfilled

    success = True
    # todo currently we do not do any presteps!
    """
    steps_pre = pp.steps_pre # perform some implicit Euler steps to create proper initial conditions
    crossSection0 = crossSection_over_time(pressure0, (-steps_pre + 1) * pp.tau_pre)
    

    for n in range(steps_pre):
        if not success:
            break
        if not pp.fsi_active:
            crossSection1 = crossSection_over_time(pressure0, (n - steps_pre + 1) * pp.tau_pre)  # update extrapolation of cross section
            velocity1, pressure1, success = perform_partitioned_implicit_euler_step(velocity0, pressure0, crossSection0, crossSection1, dx, pp.tau_pre, pp.velocity_in_linearized((n-steps_pre+1)*pp.tau_pre))
        elif pp.fsi_active:  # always use monolithic solver for initialization!
            velocity1, pressure1, crossSection1, success = perform_monolithic_implicit_euler_step(velocity0, pressure0, crossSection0, dx, pp.tau_pre, pp.velocity_in_linearized((n-steps_pre+1)*pp.tau_pre))
        crossSection0 = crossSection1
        velocity0 = velocity1
        pressure0 = pressure1
    """
    for n in range(N_steps):
        if not success:
            break

        t = n*tau
        fixed_point_solver.clear(config.underrelaxation_factor * tau)

        if is_partitioned_approach(coupling_mode):
            if performs_iterations(coupling_mode):
                k_max = config.k_max_coupling
            else:
                k_max = 1

            k = 0
            error = np.inf

            crossSection1 = np.copy(crossSection0)  # initial guess for new cross section

            while error > config.e_coupling and k < k_max:  # implicit coupling: iteratively improve crossSection1
                k += 1
                if time_stepping_scheme is config.TimeStepping.ImplicitEuler:
                    velocity1, pressure1, success = perform_partitioned_implicit_euler_step(velocity0, pressure0, crossSection0, crossSection1, dx, tau, velocity_in(t+tau))
                elif time_stepping_scheme is config.TimeStepping.TrapezoidalRule:
                    velocity1, pressure1, success = perform_partitioned_implicit_trapezoidal_rule_step(velocity0, pressure0, crossSection0, crossSection1, dx, tau, velocity_in(t + tau), custom_coupling=coupling_mode is config.CouplingAlgorithm.PartitionedPythonCustomized)
                else:
                    raise Exception("unknown time stepping scheme [%s]!", time_stepping_scheme.name)

                crossSection1_tilde = solve_solid(pressure1)  # new cross section corresponding to computed pressure
                if performs_iterations(coupling_mode):
                    crossSection1, error = fixed_point_solver.iterate(crossSection1, crossSection1_tilde)
                else:
                    crossSection1 = crossSection1_tilde
            if k == config.k_max_coupling:
                raise Exception("Implicit coupling break! Error: %.4g" % error)
                success = False
        else:
            if time_stepping_scheme is config.TimeStepping.ImplicitEuler:
                velocity1, pressure1, crossSection1, success = perform_monolithic_implicit_euler_step(velocity0, pressure0, crossSection0, dx, tau, velocity_in(t+tau))
            elif time_stepping_scheme is config.TimeStepping.TrapezoidalRule:
                velocity1, pressure1, crossSection1, success = perform_monolithic_implicit_trapezoidal_rule_step(velocity0, pressure0, crossSection0, dx, tau, velocity_in(t+tau))
            else:
                raise Exception("unknown time stepping scheme [%s]!", (time_stepping_scheme.name))

        ## swap at end of timestep
        pressure0 = pressure1
        velocity0 = velocity1
        crossSection0 = crossSection1

        ## postprocessing
        metadata = {
            'created_on': str(sim_start_time),
            'tau': tau,
            'dx': dx,
            'timestepping': time_stepping_scheme.name,
            'coupling': coupling_mode.name,
            'elasticity_module': config.E,
            'length': config.L,
            'n_elem': config.n_elem,
            'inflow_frequency': config.frequency,
            'inflow_amplitude': config.ampl,
            'inflow_mean': config.u0
        }

        create_output(t, velocity0, pressure0, crossSection0, metadata, output_mode)

        if plotting_mode is config.PlottingModes.VIDEO:
            create_video(t, velocity0, pressure0, crossSection0, metadata, writer)

        ## timestep end

    if plotting_mode is config.PlottingModes.VIDEO:
        writer.finish()

    return velocity0, pressure0, crossSection0
