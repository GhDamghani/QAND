class Trajectory:
    def __init__(self, t, u, diffeq, p, u0, tspan, solver) -> None:
        self.t = t
        self.u = u
        self.diffeq = diffeq
        self.p = p
        self.u0 = u0
        self.tspan = tspan
        self.solver = solver