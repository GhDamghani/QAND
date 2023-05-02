from diffeqpy import ode
from numpy import transpose, full, nan, argmax
from julia import Main
from numpy import save as np_save, load as np_load
from zipfile import ZipFile
from functools import partial
from os import remove as os_remove
from scipy.signal import find_peaks
from numpy import diff
np_save = partial(np_save, allow_pickle = True)
np_load = partial(np_load, allow_pickle = True)


class DiffEq:
    '''Base class for Ordinary Differential Equations'''
    _io_file_names = ('diffeq', 'ndim.npy')
    def __init__(self, diffeq, ndim) -> None:
        self.diffeq = diffeq
        self.ndim = ndim
    
    def get_julia_func_name(self):
        import re
        regex = r"<PyCall.jlwrap (\S+)>"
        test_str = self.diffeq.__repr__()
        matches = re.search(regex, test_str)
        if matches:
            func_name = matches.groups()[0]
            return func_name
        else:
            raise TypeError("It's not a Julia Function (PyCall.jlwrap)")
    
    def save(self, file) -> None:
        func_name = DiffEq.get_julia_func_name(self)
        Main.eval(f'''using Serialization: serialize;
        serialize("{DiffEq._io_file_names[0]}", {func_name})''')
        np_save(DiffEq._io_file_names[1], self.ndim)
        with ZipFile(file, mode="w") as archive:
            for filename in DiffEq._io_file_names:
                archive.write(filename)
        for x in DiffEq._io_file_names:
            os_remove(x)
        
    
    @staticmethod
    def load(file):
        with ZipFile(file, mode="r") as archive:
            archive.extractall()
        diffeq = Main.eval(f'''using Serialization: deserialize;
        deserialize("{DiffEq._io_file_names[0]}")''')
        ndim = np_load(DiffEq._io_file_names[1])
        
        for x in DiffEq._io_file_names:
            os_remove(x)
        return DiffEq(diffeq, ndim)
        


    def trajectory(self, u0, tspan, p, transient=1, reltol=1E-5, abstol=1E-5, time_res=1E-2, alg_name = 'Tsit5'):
        if alg_name == 'Tsit5':
            alg = ode.Tsit5()
        prob = ode.ODEProblem(self.diffeq, u0, tspan, p)
        sol = ode.solve(prob, alg, reltol=reltol,
                        abstol=abstol, saveat=time_res)
        solver = {
            'alg_name': alg_name, 'reltol': reltol, 'abstol': abstol, 'time_res': time_res}
        if transient < tspan[1]:
            ind = argmax(sol.t >= transient)
        else:
            ind = 0
        t = sol.t[ind:]
        u = transpose(sol.u[ind:])
        traj = Trajectory(self, t, u, p, u0, tspan, solver)
        return traj

    def bifurcation(self, p_ind, p, p_array, u0, tspan, transient=1, mode='max', reltol=1E-5, abstol=1E-5, time_res=1E-2, u0_forward=False, num_points=500, tqdm=False):
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
            peaks, _ = traj.get_nonlinear_feature(mode=mode)
            num_points_ = min(num_points, peaks.shape[0])
            bif_mat[i, -num_points_:] = peaks[-num_points_:]

            if u0_forward:
                u0 = traj.u[:, -1]
        solver = traj.solver
        bif = Bifurcation(self, bif_mat, p_ind, p, p_array,
                          u0, tspan, mode, u0_forward, solver)
        return bif


class Trajectory:
    _symbols = ('x', 'y', 'z')
    _fontname = 'Times New Roman'
    _fontsize = 20
    _io_file_names = ('diffeq.qand', 't.npy', 'u.npy', 'p.npy', 'u0.npy', 'tspan.npy', 'solver.npy')

    def __init__(self, diffeq, t, u, p, u0, tspan, solver) -> None:
        self.diffeq = diffeq
        self.t = t
        self.u = u
        self.p = p
        self.u0 = u0
        self.tspan = tspan
        self.solver = solver
    
    def save(self, file) -> None:
        self.diffeq.save(Trajectory._io_file_names[0])
        np_save(Trajectory._io_file_names[1], self.t)
        np_save(Trajectory._io_file_names[2], self.u)
        np_save(Trajectory._io_file_names[3], self.p)
        np_save(Trajectory._io_file_names[4], self.u0)
        np_save(Trajectory._io_file_names[5], self.tspan)
        np_save(Trajectory._io_file_names[6], self.solver)
        with ZipFile(file, mode="w") as archive:
            for filename in Trajectory._io_file_names:
                archive.write(filename)
        for x in Trajectory._io_file_names:
            os_remove(x)
    
    @staticmethod
    def load(file):
        with ZipFile(file, mode="r") as archive:
            archive.extractall()
        diffeq = DiffEq.load(Trajectory._io_file_names[0])
        t = np_load(Trajectory._io_file_names[1])
        u = np_load(Trajectory._io_file_names[2])
        p = np_load(Trajectory._io_file_names[3])
        u0 = np_load(Trajectory._io_file_names[4])
        tspan = np_load(Trajectory._io_file_names[5])
        solver = np_load(Trajectory._io_file_names[6])
        
        for x in Trajectory._io_file_names:
            os_remove(x)
        return Trajectory(diffeq, t, u, p, u0, tspan, solver)

    def plot_trajectory(self, ax, alpha=0.8, fontname=_fontname, fontsize=_fontsize):
        assert self.diffeq.ndim == 3
        ax.plot(self.u[0, :], self.u[1, :], self.u[2, :], alpha=alpha)
        ax.set_xlabel('x', fontname=fontname, fontsize=fontsize)
        ax.set_ylabel('y', fontname=fontname, fontsize=fontsize)
        ax.set_zlabel('z', fontname=fontname, fontsize=fontsize)

    def plot_timeseries(self, axs, symbols=_symbols, fontname=_fontname, fontsize=_fontsize):
        assert self.diffeq.ndim == len(axs)
        for i, ax in enumerate(axs):
            ax.plot(self.t, self.u[i, :])
            ax.set_xlabel('t', fontname=fontname, fontsize=fontsize)
            ax.set_ylabel(symbols[i], fontname=fontname, fontsize=fontsize)

    def get_nonlinear_feature(self, u_ind=0, mode='max'):
        signal = self.u[u_ind]
        if mode == 'max':
            ind, _ = find_peaks(signal)
        if mode == 'min':
            ind, _ = find_peaks(-signal)
        if mode == 'max_interval':
            ind, _ = find_peaks(signal)
        if mode == 'min_interval':
            ind, _ = find_peaks(-signal)
        
        u_out = signal[ind]
        t_out = self.t[ind]

        if mode in ('max_interval', 'min_interval'):
            u_out = diff(u_out)
            t_out = diff(t_out)
        
        return u_out, t_out

class Bifurcation:
    _symbols = ('x', 'y', 'z')
    _fontname = 'Times New Roman'
    _fontsize = 20
    _io_file_names = ('diffeq.qand', 'bif_mat.npy', 'p_ind.npy', 'p.npy', 'p_array.npy', 'u0.npy', 'tspan.npy', 'mode.npy', 'u0_forward.npy', 'solver.npy')

    def __init__(self, diffeq, bif_mat, p_ind, p, p_array, u0, tspan, mode, u0_forward, solver) -> None:
        self.diffeq = diffeq
        self.bif_mat = bif_mat
        self.p_ind = p_ind
        self.p = p
        self.p_array = p_array
        self.u0 = u0
        self.tspan = tspan
        self.mode = mode
        self.u0_forward = u0_forward
        self.solver = solver
    
    def save(self, file) -> None:
        self.diffeq.save(Bifurcation._io_file_names[0])
        np_save(Bifurcation._io_file_names[1], self.bif_mat)
        np_save(Bifurcation._io_file_names[2], self.p_ind)
        np_save(Bifurcation._io_file_names[3], self.p)
        np_save(Bifurcation._io_file_names[4], self.p_array)
        np_save(Bifurcation._io_file_names[5], self.u0)
        np_save(Bifurcation._io_file_names[6], self.tspan)
        np_save(Bifurcation._io_file_names[7], self.mode)
        np_save(Bifurcation._io_file_names[8], self.u0_forward)
        np_save(Bifurcation._io_file_names[9], self.solver)
        with ZipFile(file, mode="w") as archive:
            for filename in Bifurcation._io_file_names:
                archive.write(filename)
        for x in Bifurcation._io_file_names:
            os_remove(x)
    
    @staticmethod
    def load(file):
        with ZipFile(file, mode="r") as archive:
            archive.extractall()
        diffeq = DiffEq.load(Bifurcation._io_file_names[0])
        bif_mat = np_load(Bifurcation._io_file_names[1])
        p_ind = np_load(Bifurcation._io_file_names[2])
        p = np_load(Bifurcation._io_file_names[3])
        p_array = np_load(Bifurcation._io_file_names[4])
        u0 = np_load(Bifurcation._io_file_names[5])
        tspan = np_load(Bifurcation._io_file_names[6])
        mode = np_load(Bifurcation._io_file_names[7])
        u0_forward = np_load(Bifurcation._io_file_names[8])
        solver = np_load(Bifurcation._io_file_names[9])
        
        for x in Bifurcation._io_file_names:
            os_remove(x)
        return Bifurcation(diffeq, bif_mat, p_ind, p, p_array, u0, tspan, mode, u0_forward, solver)

    def plot(self, ax, p_sym, alpha=0.1, ylabel=None, fontname=_fontname, fontsize=_fontsize):
        num_points = self.bif_mat.shape[1]
        for i in range(num_points):
            ax.plot(self.p_array, self.bif_mat[:, i], 'k,', alpha=alpha)

        ax.set_xlabel(p_sym, fontname=fontname, fontsize=fontsize)
        if ylabel == None:
            ylabel = f'${Bifurcation._symbols[self.p_ind]}_{{{self.mode}}}$'
        ax.set_ylabel(ylabel, fontname=fontname, fontsize=fontsize)

