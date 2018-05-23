# based on https://github.com/precice/elastictube1d and
# [1] J. Degroote, P. Bruggeman, R. Haelterman, and J. Vierendeels. Stability of a coupling technique for partitioned solvers in FSI applications. Computers & Structures, 2008.
# for time integration details see
# [2] Gresho, P. M., & Sani, R. L. (2000). Incompressible Flow and the Finite Element Method, Isothermal Laminar Flow. John Wiley & Sons. Retrieved from http://books.google.de/books?id=m_tQAAAAMAAJ

from __future__ import division

import matplotlib
import numpy as np

import configuration_file as config
from pythonCouplingHelpers.fixedPointIteration import IQNILSScheme
from pythonCouplingHelpers.solid import solve_solid, dsolve_solid

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.animation as manimation


def perform_monolithic_theta_scheme_step(velocity0, pressure0, crossSection0, dx, tau, velocity_in, theta=1):

    k = 0

    # initial guess for Newtons method
    pressure1 = np.copy(pressure0)
    velocity1 = np.copy(velocity0)

    N = pressure0.shape[0]-1

    alpha = 0 #pp.a0 / (pp.u0 + dx/tau)
    success = True

    while success:  # perform Newton iterations to solve nonlinear system of equations

        crossSection1 = solve_solid(pressure1)
        dcrossSection1 = dsolve_solid(pressure1)

        # compute residual
        res = np.zeros(2 * N + 2)

        for i in range(1,N):
            # Momentum
            res[i] = (velocity0[i] * crossSection0[i] - velocity1[i] * crossSection1[i]) * dx / tau

            res[i] += .25 * theta * (- crossSection1[i + 1] * velocity1[i] * velocity1[i + 1] - crossSection1[i] * velocity1[i] * velocity1[i + 1])
            res[i] += .25 * (1-theta) * (- crossSection0[i + 1] * velocity0[i] * velocity0[i + 1] - crossSection0[i] * velocity0[i] * velocity0[i + 1])

            res[i] += .25 * theta * (- crossSection1[i + 1] * velocity1[i] * velocity1[i] - crossSection1[i] * velocity1[i] * velocity1[i] + crossSection1[i] * velocity1[i - 1] * velocity1[i] + crossSection1[i - 1] * velocity1[i - 1] * velocity1[i])
            res[i] += .25 * (1-theta) * (- crossSection0[i + 1] * velocity0[i] * velocity0[i] - crossSection0[i] * velocity0[i] * velocity0[i] + crossSection0[i] * velocity0[i - 1] * velocity0[i] + crossSection0[i - 1] * velocity0[i - 1] * velocity0[i])

            res[i] += .25 * theta * (+ crossSection1[i - 1] * velocity1[i - 1] * velocity1[i - 1] + crossSection1[i] * velocity1[i - 1] * velocity1[i - 1])
            res[i] += .25 * (1-theta) * (+ crossSection0[i - 1] * velocity0[i - 1] * velocity0[i - 1] + crossSection0[i] * velocity0[i - 1] * velocity0[i - 1])

            res[i] += .25 * theta * (+ crossSection1[i - 1] * pressure1[i - 1] + crossSection1[i] * pressure1[i - 1] + crossSection1[i - 1] * pressure1[i] - crossSection1[i + 1] * pressure1[i] - crossSection1[i] * pressure1[i + 1] - crossSection1[i + 1] * pressure1[i + 1])
            res[i] += .25 * (1-theta) * (+ crossSection0[i - 1] * pressure0[i - 1] + crossSection0[i] * pressure0[i - 1] + crossSection0[i - 1] * pressure0[i] - crossSection0[i + 1] * pressure0[i] - crossSection0[i] * pressure0[i + 1] - crossSection0[i + 1] * pressure0[i + 1])

            # Continuity (we only care about values at n+1, see [2],p.737,eq.(3.16-25))
            res[i + N + 1] = (crossSection0[i] - crossSection1[i]) * dx / tau
            res[i + N + 1] += .25 * theta * (+ crossSection1[i - 1] * velocity1[i - 1] + crossSection1[i] * velocity1[i - 1] + crossSection1[i - 1] * velocity1[i] - crossSection1[i + 1] * velocity1[i] - crossSection1[i] * velocity1[i + 1] - crossSection1[i + 1] * velocity1[i + 1])
            res[i + N + 1] += .25 * (1-theta) * (+ crossSection0[i - 1] * velocity0[i - 1] + crossSection0[i] * velocity0[i - 1] + crossSection0[i - 1] * velocity0[i] - crossSection0[i + 1] * velocity0[i] - crossSection0[i] * velocity0[i + 1] - crossSection0[i + 1] * velocity0[i + 1])
            res[i + N + 1] += alpha * theta * (pressure1[i - 1] - 2 * pressure1[i] + pressure1[i + 1])

        # Boundary

        # Velocity Inlet is prescribed
        res[0] = velocity_in - velocity1[0]

        # Pressure Inlet is lineary interpolated
        res[N + 1] = -pressure1[0] + 2 * pressure1[1] - pressure1[2]

        # Velocity Outlet is lineary interpolated
        res[N] = -velocity1[-1] + 2 * velocity1[-2] - velocity1[-3]

        # Pressure Outlet is "non-reflecting"
        tmp2 = np.sqrt(config.c_mk ** 2 - pressure0[-1] / 2) - (velocity1[-1] - velocity0[-1]) / 4
        res[2 * N + 1] = -pressure1[-1] + 2 * (config.c_mk ** 2 - tmp2 * tmp2)

        k += 1  # Iteration Count

        # compute relative norm of residual
        norm_1 = np.sqrt(res.dot(res))
        norm_2 = np.sqrt(pressure1.dot(pressure1) + velocity1.dot(velocity1))
        norm = norm_1 / norm_2

        if norm < 1e-10 and k > 1:
            break  # Nonlinear Solver success
        elif k > config.k_max_nonlin:
            print "Nonlinear Solver break, iterations: %i, residual norm: %e\n" % (k, norm)
            velocity1[:] = np.nan
            pressure1[:] = np.nan
            crossSection1[:] = np.nan
            success = False
            break
        # else:
        # perform another iteration of newton's method

        # compute Jacobian for Newton's method
        system = np.zeros([N+N+2,N+N+2])

        for i in range(1,N):
            ### Momentum, Velocity see [1] eq. (13b) ###

            # df[i]/du[i-1], f[i] = -res[i]
            system[i][i - 1] += 0.25 * theta * (- crossSection1[i] * velocity1[i] -  crossSection1[i - 1] * velocity1[i])
            system[i][i - 1] += 0.25 * theta * (- 2 * crossSection1[i - 1] * velocity1[i - 1] - 2 * crossSection1[i] * velocity1[i - 1])
            # df[i]/du[i], f[i] = -res[i]
            system[i][i] += crossSection1[i] * dx/tau
            system[i][i] += 0.25 * theta * (+ crossSection1[i + 1] * velocity1[i + 1] + crossSection1[i] * velocity1[i + 1])
            system[i][i] += 0.25 * theta * (+ 2 * crossSection1[i + 1] * velocity1[i] + 2 * crossSection1[i] * velocity1[i] - crossSection1[i] * velocity1[i - 1] - crossSection1[i - 1] * velocity1[i - 1])
            # df[i]/du[i+1], f[i] = -res[i]
            system[i][i + 1] += 0.25 * theta * (crossSection1[i + 1] * velocity1[i] + crossSection1[i] * velocity1[i])

            ### Momentum, Pressure see [1] eq. (13b) ###

            # df[i]/dp[i-1], f[i] = -res[i]
            system[i][N + 1 + i - 1] += 0.25 * theta * (- dcrossSection1[i - 1] * velocity1[i - 1] * velocity1[i])
            system[i][N + 1 + i - 1] += 0.25 * theta * (- dcrossSection1[i - 1] * velocity1[i - 1] * velocity1[i - 1])
            system[i][N + 1 + i - 1] += 0.25 * theta * (- crossSection1[i - 1] - dcrossSection1[i - 1] * pressure1[i - 1] - crossSection1[i] - dcrossSection1[i - 1] * pressure1[i])
            # df[i]/dp[i], f[i] = -res[i]
            system[i][N + 1 + i] += velocity1[i] * dcrossSection1[i] * dx/tau
            system[i][N + 1 + i] += 0.25 * theta * (+ dcrossSection1[i] * velocity1[i] * velocity1[i + 1])
            system[i][N + 1 + i] += 0.25 * theta * (+ dcrossSection1[i] * velocity1[i] * velocity1[i] - dcrossSection1[i] * velocity1[i - 1] * velocity1[i])
            system[i][N + 1 + i] += 0.25 * theta * (- dcrossSection1[i] * velocity1[i - 1] * velocity1[i - 1])
            system[i][N + 1 + i] += 0.25 * theta * (- crossSection1[i - 1] + crossSection1[i + 1] - dcrossSection1[i] * pressure1[i - 1] + dcrossSection1[i] * pressure1[i + 1])
            # df[i]/dp[i+1], f[i] = -res[i]
            system[i][N + 1 + i + 1] += 0.25 * theta * (+ dcrossSection1[i + 1] * velocity1[i] * velocity1[i + 1])
            system[i][N + 1 + i + 1] += 0.25 * theta * (+ dcrossSection1[i + 1] * velocity1[i] * velocity1[i])
            system[i][N + 1 + i + 1] += 0.25 * theta * (+ dcrossSection1[i + 1] * pressure1[i] + crossSection1[i] + dcrossSection1[i + 1] * pressure1[i + 1] + crossSection1[i + 1])

            ### Continuity, Velocity see [1] eq. (13a) ###

            # df[i]/du[i-1], f[i] = -res[i]
            system[i + N + 1][i - 1] += 0.25 * theta * (- crossSection1[i - 1] - crossSection1[i])
            # df[i]/du[i], f[i] = -res[i]
            system[i + N + 1][i] += 0.25 * theta * (- crossSection1[i - 1] + crossSection1[i + 1])
            # df[i]/du[i+1], f[i] = -res[i]
            system[i + N + 1][i + 1] += 0.25 * theta * (+ crossSection1[i] + crossSection1[i + 1])

            # Continuity, Pressure see [1] eq. (13a)

            # dg[i]/dp[i-1], g[i] = -res[i + N + 1]
            system[i + N + 1][N + 1 + i - 1] += 0.25 * theta * (- dcrossSection1[i - 1] * velocity1[i - 1] - dcrossSection1[i - 1] * velocity1[i])
            system[i + N + 1][N + 1 + i - 1] += - alpha * theta
            # dg[i]/dp[i], g[i] = -res[i + N + 1]
            system[i + N + 1][N + 1 + i] += dcrossSection1[i] * dx / tau
            system[i + N + 1][N + 1 + i] += 0.25 * theta * (- dcrossSection1[i] * velocity1[i - 1] + dcrossSection1[i] * velocity1[i + 1])
            system[i + N + 1][N + 1 + i] += 2 * alpha * theta
            # dg[i]/dp[i+1], g[i] = -res[i + N + 1]
            system[i + N + 1][N + 1 + i + 1] += 0.25 * theta * (+ dcrossSection1[i + 1] * velocity1[i] + dcrossSection1[i + 1] * velocity1[i + 1])
            system[i + N + 1][N + 1 + i + 1] += - alpha * theta

        # Velocity Inlet is prescribed
        system[0][0] = 1
        # Pressure Inlet is linearly interpolated [1] eq. (14a)
        system[N + 1][N + 1] = 1
        system[N + 1][N + 2] = -2
        system[N + 1][N + 3] = 1
        # Velocity Outlet is linearly interpolated [1] eq. (14b)
        system[N][N] = 1
        system[N][N - 1] = -2
        system[N][N - 2] = 1
        # Pressure Outlet is Non-Reflecting [1] eq. (15)
        system[2 * N + 1][2 * N + 1] = 1
        system[2 * N + 1][N] = -(np.sqrt(config.c_mk ** 2 - pressure0[-1] / 2) - (velocity1[-1] - velocity0[-1]) / 4)

        try:
            solution = np.linalg.solve(system, res)
        except np.linalg.LinAlgError:
            print "LINALGERROR! SINGULAR MATRIX"
            velocity1[:] = np.nan
            pressure1[:] = np.nan
            crossSection1[:] = np.nan
            success = False
            break

        velocity1 += config.relax_newton * solution[:N + 1]
        pressure1 += config.relax_newton * solution[N + 1:]

    return velocity1, pressure1, crossSection1, success


def perform_partitioned_theta_scheme_step(velocity0, pressure0, crossSection0, crossSection1, dx, tau, velocity_in, custom_coupling, theta=1):

    k = 0

    # initial guess for Newtons method
    pressure1 = np.copy(pressure0)
    velocity1 = np.copy(velocity0)

    crossSection_couple = 2 * [None]
    if custom_coupling:
        # set cross sections corresponding to point in time
        crossSection_couple[0] = crossSection0
        crossSection_couple[1] = crossSection1
    else:
        # set both cross sections equal to input -> depending on input: implicit or explicit coupling
        crossSection_couple[0] = crossSection1
        crossSection_couple[1] = crossSection1

    N = pressure0.shape[0]-1

    alpha = 0#pp.a0 / (pp.u0 + dx/tau)
    success = True

    while success:  # perform Newton iterations to solve nonlinear system of equations

        # compute residual
        res = np.zeros(2 * N + 2)

        for i in range(1,N):
            # Momentum
            res[i] = (velocity0[i] * crossSection0[i] - velocity1[i] * crossSection1[i]) * dx / tau

            res[i] += .25 * theta * (- crossSection_couple[1][i + 1] * velocity1[i] * velocity1[i + 1] - crossSection_couple[1][i] * velocity1[i] * velocity1[i + 1])
            res[i] += .25 * (1-theta) * (- crossSection_couple[0][i + 1] * velocity0[i] * velocity0[i + 1] - crossSection_couple[0][i] * velocity0[i] * velocity0[i + 1])

            res[i] += .25 * theta * (- crossSection_couple[1][i + 1] * velocity1[i] * velocity1[i] - crossSection_couple[1][i] * velocity1[i] * velocity1[i] + crossSection_couple[1][i] * velocity1[i - 1] * velocity1[i] + crossSection_couple[1][i - 1] * velocity1[i - 1] * velocity1[i])
            res[i] += .25 * (1-theta) * (- crossSection_couple[0][i + 1] * velocity0[i] * velocity0[i] - crossSection_couple[0][i] * velocity0[i] * velocity0[i] + crossSection_couple[0][i] * velocity0[i - 1] * velocity0[i] + crossSection_couple[0][i - 1] * velocity0[i - 1] * velocity0[i])

            res[i] += .25 * theta * (+ crossSection_couple[1][i - 1] * velocity1[i - 1] * velocity1[i - 1] + crossSection_couple[1][i] * velocity1[i - 1] * velocity1[i - 1])
            res[i] += .25 * (1-theta) * (+ crossSection_couple[0][i - 1] * velocity0[i - 1] * velocity0[i - 1] + crossSection_couple[0][i] * velocity0[i - 1] * velocity0[i - 1])

            res[i] += .25 * theta * (+ crossSection_couple[1][i - 1] * pressure1[i - 1] + crossSection_couple[1][i] * pressure1[i - 1] + crossSection_couple[1][i - 1] * pressure1[i] - crossSection_couple[1][i + 1] * pressure1[i] - crossSection_couple[1][i] * pressure1[i + 1] - crossSection_couple[1][i + 1] * pressure1[i + 1])
            res[i] += .25 * (1-theta) * (+ crossSection_couple[0][i - 1] * pressure0[i - 1] + crossSection_couple[0][i] * pressure0[i - 1] + crossSection_couple[0][i - 1] * pressure0[i] - crossSection_couple[0][i + 1] * pressure0[i] - crossSection_couple[0][i] * pressure0[i + 1] - crossSection_couple[0][i + 1] * pressure0[i + 1])

            # Continuity (we only care about values at n+1, see [2],p.737,eq.(3.16-25))
            res[i + N + 1] = (crossSection0[i] - crossSection1[i]) * dx / tau
            res[i + N + 1] += .25 * theta * (+ crossSection_couple[1][i - 1] * velocity1[i - 1] + crossSection_couple[1][i] * velocity1[i - 1] + crossSection_couple[1][i - 1] * velocity1[i] - crossSection_couple[1][i + 1] * velocity1[i] - crossSection_couple[1][i] * velocity1[i + 1] - crossSection_couple[1][i + 1] * velocity1[i + 1])
            res[i + N + 1] += .25 * (1-theta) * (+ crossSection_couple[0][i - 1] * velocity0[i - 1] + crossSection_couple[0][i] * velocity0[i - 1] + crossSection_couple[0][i - 1] * velocity0[i] - crossSection_couple[0][i + 1] * velocity0[i] - crossSection_couple[0][i] * velocity0[i + 1] - crossSection_couple[0][i + 1] * velocity0[i + 1])
            res[i + N + 1] += alpha * theta * (pressure1[i - 1] - 2 * pressure1[i] + pressure1[i + 1])

        # Boundary

        # Velocity Inlet is prescribed
        res[0] = velocity_in - velocity1[0]

        # Pressure Inlet is lineary interpolated
        res[N + 1] = -pressure1[0] + 2 * pressure1[1] - pressure1[2]

        # Velocity Outlet is lineary interpolated
        res[N] = -velocity1[-1] + 2 * velocity1[-2] - velocity1[-3]

        # Pressure Outlet is "non-reflecting"
        in_sqrt = config.c_mk ** 2 - pressure0[-1] / 2
        assert in_sqrt > 0
        tmp2 = np.sqrt(in_sqrt) - (velocity1[-1] - velocity0[-1]) / 4
        res[2 * N + 1] = -pressure1[-1] + 2 * (config.c_mk ** 2 - tmp2 * tmp2)

        k += 1  # Iteration Count

        # compute relative norm of residual
        norm_1 = np.sqrt(res.dot(res))
        norm_2 = np.sqrt(pressure1.dot(pressure1) + velocity1.dot(velocity1))
        norm = norm_1 / norm_2

        if norm < 1e-10 and k > 1:
            break  # Nonlinear Solver success
        elif k > config.k_max_nonlin:
            print "Nonlinear Solver break, iterations: %i, residual norm: %e\n" % (k, norm)
            velocity1[:] = np.nan
            pressure1[:] = np.nan
            success = False
            break
        # else:
        # perform another iteration of newton's method

        # compute Jacobian for Newton's method
        system = np.zeros([N+N+2,N+N+2])

        for i in range(1,N):
            # Momentum, Velocity see [1] eq. (13b)
            system[i][i - 1] += .25 * theta * (- 2 * crossSection_couple[1][i - 1] * velocity1[i - 1] - 2 * crossSection_couple[1][i] * velocity1[i - 1] - crossSection_couple[1][i] * velocity1[i] + crossSection_couple[1][i - 1] * velocity1[i])
            system[i][i] += crossSection1[i] * dx/tau
            system[i][i] += .25 * theta * (+ crossSection_couple[1][i + 1] * velocity1[i + 1] + crossSection_couple[1][i] * velocity1[i + 1] + crossSection_couple[1][i + 1] * velocity1[i] * 2 + crossSection_couple[1][i] * velocity1[i] * 2 - crossSection_couple[1][i] * velocity1[i - 1] - crossSection_couple[1][i - 1] * velocity1[i - 1])
            system[i][i + 1] += .25 * theta * (crossSection_couple[1][i + 1] * velocity1[i] + crossSection_couple[1][i] * velocity1[i])

            # Momentum, Pressure see [1] eq. (13b)
            system[i][N + 1 + i - 1] += .25 * theta * (- crossSection_couple[1][i - 1] - crossSection_couple[1][i])
            system[i][N + 1 + i] += .25 * theta * (+ crossSection_couple[1][i - 1] - crossSection_couple[1][i + 1])
            system[i][N + 1 + i + 1] += .25 * theta * (+ crossSection_couple[1][i] + crossSection_couple[1][i + 1])

            # Continuity, Velocity see [1] eq. (13a)
            system[i + N + 1][i - 1] += .25 * theta * (- crossSection_couple[1][i - 1] - crossSection_couple[1][i])
            system[i + N + 1][i] += .25 * theta * (- crossSection_couple[1][i - 1] + crossSection_couple[1][i + 1])
            system[i + N + 1][i + 1] += .25 * theta * (+ crossSection_couple[1][i] + crossSection_couple[1][i + 1])

            # Continuity, Pressure see [1] eq. (13a)
            system[i + N + 1][N + 1 + i - 1] += - alpha * theta
            system[i + N + 1][N + 1 + i] += 2 * alpha * theta
            system[i + N + 1][N + 1 + i + 1] += - alpha * theta

        # Velocity Inlet is prescribed
        system[0][0] = 1
        # Pressure Inlet is linearly interpolated [1] eq. (14a)
        system[N + 1][N + 1] = 1
        system[N + 1][N + 2] = -2
        system[N + 1][N + 3] = 1
        # Velocity Outlet is linearly interpolated [1] eq. (14b)
        system[N][N] = 1
        system[N][N - 1] = -2
        system[N][N - 2] = 1

        # Pressure Outlet is Non-Reflecting [1] eq. (15)
        system[2 * N + 1][2 * N + 1] = 1
        in_sqrt = config.c_mk ** 2 - pressure0[-1] / 2
        assert in_sqrt > 0
        system[2 * N + 1][N] = -1*(np.sqrt(in_sqrt) - (velocity1[-1] - velocity0[-1]) / 4)

        try:
            solution = np.linalg.solve(system, res)
        except np.linalg.LinAlgError:
            print "LINALGERROR! SINGULAR MATRIX"
            velocity1[:] = np.nan
            pressure1[:] = np.nan
            success = False
            break

        velocity1 += config.relax_newton * solution[:N + 1]
        pressure1 += config.relax_newton * solution[N + 1:]

    return velocity1, pressure1, success


def perform_monolithic_implicit_euler_step(velocity0, pressure0, crossSection0, dx, tau, velocity_in):
    return perform_monolithic_theta_scheme_step(velocity0, pressure0, crossSection0, dx, tau, velocity_in, theta=1)


def perform_monolithic_implicit_trapezoidal_rule_step(velocity0, pressure0, crossSection0, dx, tau, velocity_in):
    return perform_monolithic_theta_scheme_step(velocity0, pressure0, crossSection0, dx, tau, velocity_in, theta=.5)


def perform_partitioned_implicit_euler_step(velocity0, pressure0, crossSection0, crossSection1, dx, tau, velocity_in):
    return perform_partitioned_theta_scheme_step(velocity0, pressure0, crossSection0, crossSection1, dx, tau, velocity_in, custom_coupling=True, theta=1)


def perform_partitioned_implicit_trapezoidal_rule_step(velocity0, pressure0, crossSection0, crossSection1, dx, tau, velocity_in, custom_coupling):
    return perform_partitioned_theta_scheme_step(velocity0, pressure0, crossSection0, crossSection1, dx, tau, velocity_in, custom_coupling, theta=.5)


def solve_1DTube(N=config.n_elem, tau=config.tau0, T_max=config.T_max, L=config.L, velocity_in=config.velocity_in, coupling_mode=config.CouplingAlgorithm.Monolitic, time_stepping_scheme=config.TimeStepping.ImplicitEuler, plotting_mode=config.PlottingModes.OFF):
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
        writer.setup(fig, "writer_test_" + config.coupling_mode + ".mp4", 100)

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

        if (coupling_mode is config.CouplingAlgorithm.PartitionedPythonExplicit) or \
                (coupling_mode is config.CouplingAlgorithm.PartitionedPythonImplicit) or \
                (coupling_mode is config.CouplingAlgorithm.PartitionedPythonCustomized):
            if (coupling_mode is config.CouplingAlgorithm.PartitionedPythonImplicit) or \
                    (coupling_mode is config.CouplingAlgorithm.PartitionedPythonCustomized):
                k_max = config.k_max_coupling
            elif coupling_mode is config.CouplingAlgorithm.PartitionedPythonExplicit:
                k_max = 1
            else:
                raise Exception("unknown coupling algorithm!")

            k = 0
            error = np.inf

            crossSection1 = np.copy(crossSection0)  # initial guess for new cross section

            while error > config.e_coupling and k < k_max:  # implicit coupling: iteratively improve crossSection1
                k += 1
                if time_stepping_scheme is config.TimeStepping.ImplicitEuler:
                    velocity1, pressure1, success = perform_partitioned_implicit_euler_step(velocity0, pressure0, crossSection0, crossSection1, dx, tau, velocity_in(t+tau))
                elif (time_stepping_scheme is config.TimeStepping.TrapezoidalRule) or (time_stepping_scheme is config.TimeStepping.TrapezoidalRuleCustom):
                    velocity1, pressure1, success = perform_partitioned_implicit_trapezoidal_rule_step(velocity0, pressure0, crossSection0, crossSection1, dx, tau, velocity_in(t + tau), custom_coupling=coupling_mode is config.CouplingAlgorithm.PartitionedPythonCustomized)
                else:
                    raise Exception("unknown time stepping scheme!")

                crossSection1_tilde = solve_solid(pressure1)  # new cross section corresponding to computed pressure
                if (coupling_mode is config.CouplingAlgorithm.PartitionedPythonImplicit) or (coupling_mode is config.CouplingAlgorithm.PartitionedPythonCustomized):
                    crossSection1, error = fixed_point_solver.iterate(crossSection1, crossSection1_tilde)
                elif coupling_mode is config.CouplingAlgorithm.PartitionedPythonExplicit:
                    crossSection1 = crossSection1_tilde
                else:
                    raise Exception("unknown coupling algorithm!")

            if k == config.k_max_coupling:
                raise Exception("Implicit coupling break! Error: %.4g" % error)
                success = False

        elif coupling_mode is config.CouplingAlgorithm.Monolitic:
            if time_stepping_scheme is config.TimeStepping.ImplicitEuler:
                velocity1, pressure1, crossSection1, success = perform_monolithic_implicit_euler_step(velocity0, pressure0, crossSection0, dx, tau, velocity_in(t+tau))
            elif (time_stepping_scheme is config.TimeStepping.TrapezoidalRule) or (time_stepping_scheme is config.TimeStepping.TrapezoidalRuleCustom):
                velocity1, pressure1, crossSection1, success = perform_monolithic_implicit_trapezoidal_rule_step(velocity0, pressure0, crossSection0, dx, tau, velocity_in(t+tau))
            else:
                raise Exception("unknown time stepping scheme!")
        else:
            raise Exception("unknown coupling algorithm!")

        if plotting_mode is config.PlottingModes.DEBUG:
            tubePlotting.doPlotting(ax, crossSection0, velocity0, pressure0, dx, t)
        elif plotting_mode is config.PlottingModes.VIDEO:
            tubePlotting.doPlotting(ax, crossSection0, velocity0, pressure0, dx, t)
            writer.grab_frame()
            ax.cla()

        # swap at end of timestep
        pressure0 = pressure1
        velocity0 = velocity1
        crossSection0 = crossSection1

    if plotting_mode is config.PlottingModes.VIDEO:
        writer.finish()

    return velocity0, pressure0, crossSection0
