import numpy as np
import matplotlib.pyplot as plt

GM = 4.0 * np.pi**2

def derivatives(state):
    x, y, vx, vy = state
    r = np.sqrt(x**2 + y**2)
    ax = -GM * x / (r**3)
    ay = -GM * y / (r**3)
    return np.array([vx, vy, ax, ay])

def rk4_step(state, dt):
    k1 = derivatives(state)
    k2 = derivatives(state + 0.5 * dt * k1)
    k3 = derivatives(state + 0.5 * dt * k2)
    k4 = derivatives(state + dt * k3)
    return state + (dt / 6.0) * (k1 + 2*k2 + 2*k3 + k4)

r0 = 1.0
v0 = 5.5
state = np.array([r0, 0.0, 0.0, v0])

dt = 0.001
total_time = 2.0
steps = int(total_time / dt)

trajectory = np.zeros((steps, 2))
energy = np.zeros(steps)

for i in range(steps):
    trajectory[i] = [state[0], state[1]]
    r = np.sqrt(state[0]**2 + state[1]**2)
    v_sq = state[2]**2 + state[3]**2
    energy[i] = 0.5 * v_sq - GM / r
    state = rk4_step(state, dt)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 5), dpi=300)

ax1.plot(trajectory[:, 0], trajectory[:, 1], color='navy', label='Orbital Path (RK4)')
ax1.plot(0, 0, marker='o', markersize=12, color='gold', markeredgecolor='darkorange', label='Host Star')
ax1.scatter(trajectory[0, 0], trajectory[0, 1], color='forestgreen', s=40, zorder=5, label='Perihelion')
ax1.set_xlabel('X Position (AU)', fontsize=10)
ax1.set_ylabel('Y Position (AU)', fontsize=10)
ax1.set_title('Planetary Trajectory', fontsize=12)
ax1.axis('equal')
ax1.grid(True, linestyle='--', alpha=0.5)
ax1.legend(loc='upper right', fontsize=8)

time_line = np.linspace(0, total_time, steps)
fractional_energy_error = (energy - energy[0]) / np.abs(energy[0])
ax2.plot(time_line, fractional_energy_error, color='darkred')
ax2.set_xlabel('Time (Years)', fontsize=10)
ax2.set_ylabel(r'$\Delta E / |E_0|$', fontsize=10)
ax2.set_title('Energy Conservation Stability', fontsize=12)
ax2.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax2.grid(True, linestyle='--', alpha=0.5)

plt.tight_layout()
plt.savefig('orbit_simulation.png')
plt.show()
print("Project 2 Complete! Image saved.")
