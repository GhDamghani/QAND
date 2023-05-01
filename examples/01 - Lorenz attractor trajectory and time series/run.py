import qand
import matplotlib.pyplot as plt
from julia import Main

f = Main.eval(open('lorenz.jl', 'r').read())

sys1 = qand.DiffEq(f, 3)
u0 = [1.0, 0.0, 0.0]
tspan = (0., 100.)
p = [10.0, 28.0, 8/3]

sol = sys1.trajectory(u0, tspan, p)

fig = plt.figure(figsize=(8, 8))
plt.suptitle('Lorenz Attractor Trajectory Sample')
ax = fig.add_subplot(111, projection='3d')
sol.plot_trajectory(ax)
fig.savefig('Trajectory.png')

fig, axs = plt.subplots(3, figsize=(8, 8))
plt.suptitle('Lorenz Attractor Time Series Sample')
sol.plot_timeseries(axs)
fig.savefig('Time Series.png')

plt.show()
