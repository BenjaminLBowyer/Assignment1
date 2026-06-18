import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

def model(y, t,
          Lambda, mu, beta, gamma,
          delta, delta_s, r,
          p, p_s, c):

    T, I, S, V = y

    dTdt = (
        Lambda
        - mu*T
        - beta*V*T
        - gamma*(S + I)*T
    )

    dIdt = (
        beta*V*T
        - gamma*(T + 2*I + S)*I
        - delta*I
    )

    dSdt = (
        gamma*(T + I)*(2*I + S)
        - r*delta_s*S
    )

    dVdt = (
        p*I
        + p_s*S
        - c*V
    )

    return [dTdt, dIdt, dSdt, dVdt]


Lambda = float(input("Lambda: "))
mu = float(input("mu: "))
beta = float(input("beta: "))
gamma = float(input("gamma: "))

delta = float(input("delta: "))
delta_s = float(input("delta_s: "))
r = float(input("r: "))

p = float(input("p: "))
p_s = float(input("p_s: "))
c = float(input("c: "))

T0 = float(input("Initial T: "))
I0 = float(input("Initial I: "))
S0 = float(input("Initial S: "))
V0 = float(input("Initial V: "))

start = float(input("Start time: "))
end = float(input("End time: "))
points = int(input("Number of time points: "))

y0 = [T0, I0, S0, V0]

t = np.linspace(start, end, points)

solution = odeint(
    model,
    y0,
    t,
    args=(
        Lambda, mu, beta, gamma,
        delta, delta_s, r,
        p, p_s, c
    )
)

T = solution[:,0]
I = solution[:,1]
S = solution[:,2]
V = solution[:,3]

plt.plot(t, T, label='Target Cells')
plt.plot(t, I, label='Infected Cells')
plt.plot(t, S, label='Syncytia')
plt.semilogy(t, V, label='Virus')

plt.xlabel('Time')
plt.ylabel('Population')
plt.legend()
plt.grid()
plt.show()