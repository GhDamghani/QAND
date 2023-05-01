from scipy.signal import find_peaks


class Trajectory:
    symbols = ('x', 'y', 'z')
    fontname = 'Times New Roman'
    fontsize = 20

    def __init__(self, DiffEq, t, u, p, u0, tspan, solver) -> None:
        self.DiffEq = DiffEq
        self.t = t
        self.u = u
        self.p = p
        self.u0 = u0
        self.tspan = tspan
        self.solver = solver

    def plot_trajectory(self, ax, alpha=0.8, fontname=fontname, fontsize=fontsize):
        assert self.DiffEq.ndim == 3
        ax.plot(self.u[0, :], self.u[1, :], self.u[2, :], alpha=alpha)
        ax.set_xlabel('x', fontname=fontname, fontsize=fontsize)
        ax.set_ylabel('y', fontname=fontname, fontsize=fontsize)
        ax.set_zlabel('z', fontname=fontname, fontsize=fontsize)

    def plot_timeseries(self, axs, symbols=symbols, fontname=fontname, fontsize=fontsize):
        assert self.DiffEq.ndim == len(axs)
        for i, ax in enumerate(axs):
            ax.plot(self.t, self.u[i, :])
            ax.set_xlabel('t', fontname=fontname, fontsize=fontsize)
            ax.set_ylabel(symbols[i], fontname=fontname, fontsize=fontsize)

    def get_peaks(self, u_ind=0, mode='max'):
        signal = self.u[u_ind]
        if mode == 'max':
            peaks, _ = find_peaks(signal)
        if mode == 'min':
            peaks, _ = find_peaks(-signal)
        return signal[peaks], self.t[peaks]
