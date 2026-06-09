import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

T0 = 1
I0 = 0
V0 = .01

B = float(input("Input the frequency of infection per hour: "))
O = float(input("Input the death rate of infected cells: "))
P = float(input("Input the virus production rate: "))
C = float(input("Input the clearing rate: "))

time0 = float(input("Input the time at start of sim: "))
time  = float(input("Input the time at the end of the sim: "))
num   = int(input("Input the number of time plot points: "))

def step(y, t, B, O, P, C):
    T, I, V = y
    dTdt = -B * T * V
    dIdt = B * T * V - O * I
    dVdt = P * I - C * V
    return [dTdt, dIdt, dVdt]

y0 = [T0, I0, V0]
t = np.linspace(time0, time, num)
sol = odeint(step, y0, t, args=(B, O, P, C))

T = sol[:, 0]
I = sol[:, 1]
V = sol[:, 2]

fig, ax1 = plt.subplots(figsize=(10, 6))

ax1.plot(t, T, color='blue',   label="Target Cells (T)")
ax1.plot(t, I, color='orange', label="Infected Cells (I)")
ax1.set_xlabel("Time (hours)")
ax1.set_ylabel("Cells (linear scale)", color='black')
ax1.tick_params(axis='y')

ax2 = ax1.twinx()
ax2.semilogy(t, V, color='green', label="Virus (V)")
ax2.set_ylabel("Virus Load (log scale)", color='green')
ax2.tick_params(axis='y', labelcolor='green')

lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper right')

plt.title("Viral Infection Model")
plt.grid(True)
plt.tight_layout()
plt.show()