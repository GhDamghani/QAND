import qand
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from julia import Main

f = Main.eval(open('lorenz.jl', 'r').read())

sys1 = qand.DiffEq(f, 3)
u0 = [1.0,0.0,0.0]
tspan = (0., 100.)
p = [10.0,28.0,8/3]
reltol = abstol = 1E-5
saveat = 0.01

sol = sys1.solve(u0, tspan, p, reltol, abstol, saveat)

ut = np.transpose(sol.u)

fig = plt.figure(figsize=(12, 12))
ax = fig.add_subplot(111, projection='3d')
ax.plot(ut[0,:],ut[1,:],ut[2,:], alpha=0.6)
plt.show()