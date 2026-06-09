import numpy
import matplotlib.pyplot
from scipy.integrate import odeint

T0 = 1
I0 = 0
V0 = .01

B = input("Input the frequency of infection per hour: ")
O = input("Input the death rate of infected cells: ")
P = input("Input the virus production rate: ")
C = input("Input the clearing rate: ")

time0 = input("Input the time at start of sim: ")
time = input("Input the time at the end of the sim: ")

def step( y, B, O, P, C ):
    T, I, V = y

    dTdt = -B*T*V
    dIdt = B*T*V - O*I
    dVdt = P*I - C*V

    return [ dTdt, dIdt, dVdt ]

y0 = [T0, I0, V0]
t = numpy.linspace(time0, time, 100)
sol = odeint(step, y0, t, args =(B, O, P, C))

T = step[:, 0]
I = step[:, 1]
V = step[:, 2]


