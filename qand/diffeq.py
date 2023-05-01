from diffeqpy import ode
from . trajectory import Trajectory
from . bifurcation import Bifurcation
from numpy import transpose, full, nan, argmax


class DiffEq:
    '''Base class for Ordinary Differential Equations'''

    def __init__(self, diffeq, ndim) -> None:
        self.diffeq = diffeq
        self.ndim = ndim

    def trajectory(self, u0, tspan, p, transient=1, reltol=1E-5, abstol=1E-5, time_res=1E-2, alg=ode.Tsit5()) -> Trajectory:
        prob = ode.ODEProblem(self.diffeq, u0, tspan, p)
        sol = ode.solve(prob, alg, reltol=reltol,
                        abstol=abstol, saveat=time_res)
        solver = {
            'alg': alg, 'reltol': reltol, 'abstol': abstol, 'time_res': time_res}
        if transient < tspan[1]:
            ind = argmax(sol.t >= transient)
        else:
            ind = 0
        t = sol.t[ind:]
        u = transpose(sol.u[ind:])
        traj = Trajectory(self, t, u, p, u0, tspan, solver)
        return traj

    def bifurcation(self, p_ind, p, p_array, u0, tspan, transient=1, mode='max', reltol=1E-5, abstol=1E-5, time_res=1E-2, u0_forward=False, num_points=500, tqdm=False) -> Bifurcation:
        assert p_array.ndim == 1
        p_array_len = p_array.shape[0]
        bif_mat = full((p_array_len, num_points), nan)
        p_array_range = range(p_array_len)
        if tqdm:
            from tqdm import tqdm as tqdm_
            p_array_range = tqdm_(p_array_range)

        for i in p_array_range:
            p[p_ind] = p_array[i]

            traj = DiffEq.trajectory(self, u0, tspan, p, transient=transient, reltol=reltol,
                                     abstol=abstol, time_res=time_res)
            peaks, _ = traj.get_peaks()
            num_points_ = min(num_points, peaks.shape[0])
            bif_mat[i, -num_points_:] = peaks[-num_points_:]

            if u0_forward:
                u0 = traj.u[-1]
        solver = traj.solver
        bif = Bifurcation(self, bif_mat, p_ind, p, p_array,
                          u0, tspan, mode, u0_forward, solver)
        return bif
