class Bifurcation:
    symbols = ('x', 'y', 'z')
    fontname = 'Times New Roman'
    fontsize = 20

    def __init__(self, DiffEq, bif_mat, p_ind, p, p_array, u0, tspan, mode, u0_forward, solver) -> None:
        self.DiffEq = DiffEq
        self.bif_mat = bif_mat
        self.p_ind = p_ind
        self.p = p
        self.p_array = p_array
        self.u0 = u0
        self.tspan = tspan
        self.mode = mode
        self.u0_forward = u0_forward
        self.solver = solver

    def plot(self, ax, p_sym, alpha=0.1, ylabel=None, fontname=fontname, fontsize=fontsize):
        num_points = self.bif_mat.shape[1]
        for i in range(num_points):
            ax.plot(self.p_array, self.bif_mat[:, i], 'k,', alpha=alpha)

        ax.set_xlabel(p_sym, fontname=fontname, fontsize=fontsize)
        if ylabel == None:
            ylabel = f'${Bifurcation.symbols[self.p_ind]}_{{{self.mode}}}$'
        ax.set_ylabel(ylabel, fontname=fontname, fontsize=fontsize)
