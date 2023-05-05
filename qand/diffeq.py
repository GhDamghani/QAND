from diffeqpy import ode
from numpy import transpose, full, nan, argmax
from julia import Main
from numpy import save as np_save, load as np_load, savez_compressed as np_savez
from zipfile import ZipFile
from functools import partial
from os import getcwd, chdir, listdir, remove
from scipy.signal import find_peaks
from numpy import diff
import warnings
import tempfile
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
        with tempfile.TemporaryDirectory(suffix='DiffEq') as tmpdirname:
            func_name = self.get_julia_func_name()
            curdirname = getcwd()
            chdir(tmpdirname)
            Main.eval(f'''using Serialization: serialize;
            serialize("{DiffEq._io_file_names[0]}", {func_name})''')
            np_save(DiffEq._io_file_names[1], self.ndim)
            chdir(curdirname)
            with ZipFile(file, mode="w") as archive:
                chdir(tmpdirname)
                for filename in DiffEq._io_file_names:
                    archive.write(filename)
            chdir(curdirname)
        
        
    
    @classmethod
    def load(cls, file):
        with tempfile.TemporaryDirectory(suffix='DiffEq') as tmpdirname:
            curdirname = getcwd()
            with ZipFile(file, mode="r") as archive:
                chdir(tmpdirname)
                archive.extractall()
            diffeq = Main.eval(f'''using Serialization: deserialize;
            deserialize("{DiffEq._io_file_names[0]}")''')
            ndim = np_load(DiffEq._io_file_names[1])
            chdir(curdirname)
        
        return cls(diffeq, ndim)
        


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
        if t[-1] != tspan[-1]:
            warnings.warn(f'System {self.get_julia_func_name()} with this configuration is unbounded, list timepoint saved: {sol.t[-1]}', UserWarning)
        traj = Trajectory(self, t=t, u=u, p=p, u0=u0, tspan=tspan, solver=solver)
        return traj

    @staticmethod
    def _bif_nextstep_u0(i, p_array_len, u0_forward, u0_array, u):
        if i != p_array_len-1:
                if u0_forward:
                    u0_array[i+1] = u[:, -1]
                else:
                    u0_array[i+1] = u0_array[i]

    
    def bifurcation(self, p_ind, p, p_array, u0, tspan, transient=1, mode='max', reltol=1E-4, abstol=1E-4, time_res=1E-2, u0_forward=False, num_points=500, tqdm=False, print_max_num_points=False):
        assert p_array.ndim == 1
        p_array_len = p_array.shape[0]
        bif_mat = full((p_array_len, num_points), nan)
        u0_array = full((p_array_len, self.ndim), nan)
        u0_array[0] = u0.copy()
        p_array_range = range(p_array_len)
        if print_max_num_points:
            max_num_points = 0
        if tqdm:
            from tqdm import tqdm as tqdm_
            p_array_range = tqdm_(p_array_range)

        for i in p_array_range:
            p[p_ind] = p_array[i]

            traj = DiffEq.trajectory(self, u0_array[i], tspan, p, transient=transient, reltol=reltol,
                                     abstol=abstol, time_res=time_res)
            if traj.t[-1] != tspan[-1]:
                break
                # self._bif_nextstep_u0(i, p_array_len, u0_forward, u0_array, traj.u)
                # continue
            points, _ = traj.get_nonlinear_feature(mode=mode)
            points_n = points.shape[0]
            if print_max_num_points:
                if points_n>max_num_points:
                    max_num_points = points_n
                    print_str = f'Current max num points: {max_num_points} at p={p[p_ind]}, ind={i}/{p_array_len}'
                    if tqdm:
                        tqdm_.write(print_str)
                    else:
                        print(print_str)
            num_points_ = min(num_points, points_n)
            bif_mat[i, -num_points_:] = points[-num_points_:]

            self._bif_nextstep_u0(i, p_array_len, u0_forward, u0_array, traj.u)
        solver = traj.solver
        bif = Bifurcation(self, bif_mat=bif_mat, p_ind=p_ind, p=p, p_array=p_array,
                          u0_array=u0_array, tspan=tspan, mode=mode, u0_forward=u0_forward, solver=solver)
        return bif


class Trajectory:
    _symbols = ('x', 'y', 'z')
    _fontname = 'Times New Roman'
    _fontsize = 20
    _io_file_names = ('diffeq.qand', 'parameters.npz')

    def __init__(self, diffeq, t, u, p, u0, tspan, solver) -> None:
        self.diffeq = diffeq
        self.t = t
        self.u = u
        self.p = p
        self.u0 = u0
        self.tspan = tspan
        self.solver = solver
    
    def save(self, file) -> None:
        with tempfile.TemporaryDirectory(suffix='Trajectory') as tmpdirname:
            curdirname = getcwd()
            chdir(tmpdirname)
            self.diffeq.save(Trajectory._io_file_names[0])
            np_savez(Trajectory._io_file_names[1], t=self.t, u=self.u, p=self.p, u0=self.u0, tspan=self.tspan, solver=self.solver)
            chdir(curdirname)
            with ZipFile(file, mode="w") as archive:
                chdir(tmpdirname)
                for filename in Trajectory._io_file_names:
                    archive.write(filename)
            chdir(curdirname)


    @classmethod
    def load(cls, file):
        with tempfile.TemporaryDirectory(suffix='Trajectory') as tmpdirname:
            curdirname = getcwd()
            with ZipFile(file, mode="r") as archive:
                chdir(tmpdirname)
                archive.extractall()
                
            diffeq = DiffEq.load(Trajectory._io_file_names[0])
            parameters = np_load(Trajectory._io_file_names[1])
            chdir(curdirname)

            t = parameters['t']
            u = parameters['u']
            p = parameters['p']
            u0 = parameters['u0']
            tspan = parameters['tspan']
            solver = parameters['solver']
            parameters.close()
            
        return cls(diffeq=diffeq, t=t, u=u, p=p, u0=u0, tspan=tspan, solver=solver)

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
    _io_file_names = ('diffeq.qand', 'parameters.npz')

    def __init__(self, diffeq, bif_mat, p_ind, p, p_array, u0_array, tspan, mode, u0_forward, solver) -> None:
        self.diffeq = diffeq
        self.bif_mat = bif_mat
        self.p_ind = p_ind
        self.p = p
        self.p_array = p_array
        self.u0_array = u0_array
        self.tspan = tspan
        self.mode = mode
        self.u0_forward = u0_forward
        self.solver = solver
    
    def save(self, file) -> None:
        with tempfile.TemporaryDirectory(suffix='Bifurcation') as tmpdirname:
            curdirname = getcwd()
            chdir(tmpdirname)
            self.diffeq.save(Bifurcation._io_file_names[0])
            np_savez(Bifurcation._io_file_names[1], bif_mat=self.bif_mat, p_ind=self.p_ind, p=self.p, p_array=self.p_array, u0_array=self.u0_array, tspan=self.tspan, mode=self.mode, u0_forward=self.u0_forward, solver=self.solver)
            chdir(curdirname)
            with ZipFile(file, mode="w") as archive:
                chdir(tmpdirname)
                for filename in Bifurcation._io_file_names:
                    archive.write(filename)
            chdir(curdirname)
        
    
    @staticmethod
    def load(file):
        with tempfile.TemporaryDirectory(suffix='Bifurcation') as tmpdirname:
            curdirname = getcwd()
            with ZipFile(file, mode="r") as archive:
                chdir(tmpdirname)
                archive.extractall()

            diffeq = DiffEq.load(Bifurcation._io_file_names[0])
            parameters = np_load(Bifurcation._io_file_names[1])
            chdir(curdirname)
        

            bif_mat = parameters['bif_mat']
            p_ind = parameters['p_ind']
            p = parameters['p']
            p_array = parameters['p_array']
            u0_array = parameters['u0_array']
            tspan = parameters['tspan']
            mode = parameters['mode']
            u0_forward = parameters['u0_forward']
            solver = parameters['solver']
            parameters.close()
        
        return Bifurcation(diffeq, bif_mat, p_ind, p, p_array, u0_array, tspan, mode, u0_forward, solver)

    def plot(self, ax, p_sym, alpha=0.1, ylabel=None, fontname=_fontname, fontsize=_fontsize):
        num_points = self.bif_mat.shape[1]
        for i in range(num_points):
            ax.plot(self.p_array, self.bif_mat[:, i], 'k,', alpha=alpha)

        ax.set_xlabel(p_sym, fontname=fontname, fontsize=fontsize)
        if ylabel == None:
            ylabel = f'${Bifurcation._symbols[self.p_ind]}_{{{self.mode}}}$'
        ax.set_ylabel(ylabel, fontname=fontname, fontsize=fontsize)

