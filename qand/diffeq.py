from diffeqpy import ode
from . trajectory import Trajectory

class DiffEq:
    '''Base class for Ordinary Differential Equations'''
    def __init__(self, diffeq, ndim) -> None:
        self.diffeq = diffeq
        self.ndim = ndim
    
    def solve(self, u0, tspan, p, reltol, abstol, saveat) -> Trajectory:
        alg = ode.Tsit5()
        prob = ode.ODEProblem(self.diffeq, u0, tspan, p)
        sol = ode.solve(prob, alg, reltol=reltol, abstol=abstol, saveat=saveat)
        traj = Trajectory(sol.t, sol.u, self, p, u0, tspan, {'alg':alg, 'reltol': reltol, 'abstol': abstol, 'saveat':saveat})
        return traj