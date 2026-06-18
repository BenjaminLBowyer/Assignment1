import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

T0 = 1000000
I0 = 1
V0 = 1

Y = float(input("Y Input the Rate at which target cells generate : "))
U = float(input("U Input the natural death rate of target cells : "))
B = float(input("B Input the target cell infection rate : "))
P = float(input("P Input the virus production rate : "))
C = float(input("C Input the rate at which viruses are cleared from the system : "))
S = float(input("S Input the infected cells death rate : "))

time0 = 0
time  = float(input("How long is the sim run for: "))
num   = int(input("Input the number of time plot points: "))

def step(y, t, Y, U, B, P, C, S):
    T, I, V = y
    dTdt = Y - U * T - B * V * T
    dIdt = B * V * T - S * I
    dVdt = P * I - C * V
    return [dTdt, dIdt, dVdt]

y0 = [T0, I0, V0]
t = np.linspace(time0, time, num)
sol = odeint(step, y0, t, args=(Y, U, B, P, C, S))

T = sol[:, 0]
I = sol[:, 1]
V = sol[:, 2]

fig, ax1 = plt.subplots()

ax1.plot(t, T, label="Target Cells (T)")
ax1.plot(t, I, label="Infected Cells (I)")
ax1.set_xlabel("Time (hours)")
ax1.set_ylabel("Cells")

ax2 = ax1.twinx()
ax2.semilogy(t, V, color='green', label="Virus (V)")
ax2.set_ylabel("Virus Load (log scale)")

lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax1.legend(lines1 + lines2, labels1 + labels2)

plt.title("Viral Infection Model")
plt.show()