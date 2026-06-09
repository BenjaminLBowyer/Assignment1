import numpy
import matplotlib.pyplot
from scipy.integrate import odeint

T = 1
I = 0
V = .01

B = input("Input the frequency of infection per hour: ")
O = input("Input the death rate of infected cells: ")
P = input("Input the virus production rate: ")
C = input("Input the clearing rate: ")

