import qand
import matplotlib.pyplot as plt
from os import listdir

sys1 = qand.DiffEq('lorenz.jl', 3)
u0 = [1.0, 0.0, 0.0]
tspan = (0., 100.)
p = [10.0, 28.0, 8/3]

# Check if the trajectory has already been saved so we don't have to compute it again.
# Also we make sure to save it for later use
traj_file = 'Trajectory.qand'
current_files = listdir()
if traj_file not in current_files:
    # Computing trajectory
    sol = sys1.trajectory(u0, tspan, p)
    sol.save(traj_file)
else:
    print('Already computed, loading...')
    sol = qand.Trajectory.load(traj_file)

# Visualizing and saving the figure
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
