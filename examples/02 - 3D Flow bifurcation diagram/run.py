import qand
import matplotlib.pyplot as plt
from numpy import linspace
from julia import Main

f = Main.eval(open('ES1.jl', 'r').read())

sys1 = qand.DiffEq(f, 3)
u0 = [6., 0., -1.0]
tspan = (0., 2000.)
p = [1.54]
n = 500
p_array = linspace(1.52, 1.56, n, 1)

from os import listdir
# Check if the bifurcation has already been saved so we don't have to compute it again.
# Also we make sure to save it for later use
bif_file = 'Bifurcation.qand'
current_files = listdir()
if bif_file not in current_files:
    # Computing trajectory
    bif = sys1.bifurcation(0, p, p_array, u0, tspan, mode='min', transient=30, u0_forward=True, tqdm=True, num_points=100, print_max_num_points=True)
    bif.save(bif_file)
else:
    print('Already computed, loading...')
    bif = qand.Bifurcation.load(bif_file)

fig = plt.figure(figsize=(12, 8))
plt.suptitle('ES1 (Jafari et al, 2016) Bifurcation Diagram Sample', fontsize=24)
ax = fig.add_subplot(111)
bif.plot(ax, p_sym='a', alpha=1)
ax.set_xlim([p_array[0], 1.56])
fig.savefig('Bifurcation.png')

plt.show()
