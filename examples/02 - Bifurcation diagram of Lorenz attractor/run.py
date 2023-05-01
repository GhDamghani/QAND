import qand
import matplotlib.pyplot as plt
from numpy import linspace
from julia import Main
from math import pi

f = Main.eval(open('sprott.jl', 'r').read())

sys1 = qand.DiffEq(f, 3)
u0 = [pi, 0.0]
tspan = (0., 5000.)
p = [0.33, 0, 0.73]
n = 500
p_array = linspace(0, 2, n)

bif = sys1.bifurcation(1, p, p_array, u0, tspan, transient=10, tqdm=True, num_points=2000)
ylabel = r'$v_{max}$'
fig = plt.figure(figsize=(16, 12))
plt.suptitle('(Sprott, Julien C., et al, 2017) Bifurcation Diagram Sample', fontsize=24)
ax = fig.add_subplot(111)
bif.plot(ax, p_sym='A', ylabel=ylabel, alpha=1)
fig.savefig('Bifurcation.png')

# plt.show()
