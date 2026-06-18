import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

def model(y, t, Lambda, mu, beta, delta, p, c):

    T, I, V = y

    dTdt = Lambda - mu*T - beta*V*T
    dIdt = beta*V*T - delta*I
    dVdt = p*I - c*V

    return [dTdt, dIdt, dVdt]


Lambda = float(input("Lambda: "))
mu = float(input("mu: "))
beta = float(input("beta: "))
delta = float(input("delta: "))
p = float(input("p: "))
c = float(input("c: "))

T0 = float(input("Initial T: "))
I0 = float(input("Initial I: "))
V0 = float(input("Initial V: "))

start = float(input("Start time: "))
end = float(input("End time: "))
points = int(input("Number of time points: "))

y0 = [T0, I0, V0]

t = np.linspace(start, end, points)

solution = odeint(
    model,
    y0,
    t,
    args=(Lambda, mu, beta, delta, p, c)
)

T = solution[:,0]
I = solution[:,1]
V = solution[:,2]

plt.plot(t, T, label='Target Cells')
plt.plot(t, I, label='Infected Cells')
plt.semilogy(t, V, label='Virus')

plt.xlabel('Time')
plt.ylabel('Population')
plt.legend()
plt.grid()
plt.show()