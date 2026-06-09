#T
targetCells = 0

#I
infectedCells = 0

#V
virusParticles = 0

#how fast is the number of target cells changing? ( dT/ dt )
#how fast is the number of infected cells changing? ( dI/ dt )
#how fast is the number of virus population changing? ( dV/ dt )

#B / hour
frequencyOfInfection = 1e-5

#O / delta / hour
deathRateOfInfectedCells = 4

#p / hour
virusProductionRate = 2e6

#c / hour
clearanceRate = 4

#function dy / dt = func(y, t, ...)
#this would mean dT / dt ---> func(T / time)