import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint


Lambda = 0.272
mu =1.36e-3
Beta_s = 2.7e-4
Beta_ns = 2.7e-4
p_s = 100
p_ns = 100
Gamma = 1.0e-4 
Delta_s = .33
Delta_ns = .33
c_s = 2
c_ns = 2
rp = 1
rdelta = 1

Vns0 = 4.0e-7
T0 = 1000
Ins0 = 1
Is0 = 1
S0 = 0
Vs0 = 4.0e-7

Starttime = 0
Endtime = 1000
points = 10000

y0 = [T0, Ins0, Vns0, Is0, S0, Vs0]


def sim ( y, t, Lambda, mu, Beta_s, Beta_ns, Delta_s, Delta_ns, p_s, p_ns, Gamma, c_s, c_ns, rp, rdelta):

    T, Ins, Vns, Is, S, Vs, = y

    dTdt = ( Lambda - mu*T - Beta_ns*Vns*T - Beta_s*Vs*T - Gamma * (S+Is) *T )
    dInsdt = ( Beta_ns * Vns * T - Delta_ns * Ins )
    dVnsdt = ( p_ns * Ins - c_ns * Vns )
    dIsdt = ( Beta_s * Vs * T - Gamma * (T + 2 * Is + S) * Is - Delta_s * Is )
    dSdt = ( Gamma * T * (2 * Is + S) + Gamma * Is * (2 * Is + S) - rdelta * Delta_s * S )
    dVsdt = ( p_s * Is + rp * p_s * S - c_s * Vs )

    return [dTdt, dInsdt, dVnsdt, dIsdt, dSdt, dVsdt]


t = np.linspace(Starttime,Endtime,points)

solution = odeint(
    sim,
    y0,
    t,
    args=(
        Lambda, mu,
        Beta_s, Beta_ns,
        Delta_s, Delta_ns,
        p_s, p_ns,
        Gamma,
        c_s, c_ns,
        rp, rdelta
    )
)

T = solution[:, 0]
Ins = solution[:, 1]
Vns = solution[:, 2]
Is = solution[:, 3]
S = solution[:, 4]
Vs = solution[:, 5]



peak_vns = np.max(Vns)
peak_vs = np.max(Vs)

peak_time_vns = t[np.argmax(Vns)]
peak_time_vs = t[np.argmax(Vs)]

chronic_vns = Vns[-1]
chronic_vs = Vs[-1]

print("\n================ RESULTS ================\n")

print(f"Non-Syncytia Peak Viral Load : {peak_vns:.6e}")
print(f"Non-Syncytia Peak Time       : {peak_time_vns:.2f} days")
print(f"Non-Syncytia Chronic Load    : {chronic_vns:.6e}\n")

print(f"Syncytia Peak Viral Load     : {peak_vs:.6e}")
print(f"Syncytia Peak Time           : {peak_time_vs:.2f} days")
print(f"Syncytia Chronic Load        : {chronic_vs:.6e}")


fig, ax1 = plt.subplots(figsize=(12, 7))

# Cell populations (left axis)
ax1.plot(t, T, label='T (Healthy Cells)')
ax1.plot(t, Ins, label='Ins (Non-Syncytia Infected)')
ax1.plot(t, Is, label='Is (Syncytia Infected)')
ax1.plot(t, S, label='S (Syncytia)')

ax1.set_xlabel('Time (days)')
ax1.set_ylabel('Cell Population')
ax1.grid(True)

# Virus populations (right axis)
ax2 = ax1.twinx()

ax2.semilogy(t, Vns, label='Vns (Virus)')
ax2.semilogy(t, Vs, label='Vs (Virus)')

ax2.set_ylabel('Virus Concentration (log scale)')

# Combine legends from both axes
lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()

ax1.legend(
    lines1 + lines2,
    labels1 + labels2,
    loc='best'
)

plt.title('Competition Between Syncytia and Non-Syncytia HIV Strains')

plt.tight_layout()
plt.show()