"Aleksi Bossart, 27.3.2019"
"Plot the recursive decay"

import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit

def rational(x, p, q):
    return np.polyval(p, x) / np.polyval(q, x)

def rational3_3(x, p0, p1, p2, q0, q1, q2):
    return rational(x, [p0, p1, p2], [q0, q1, q2])


"Number of iterations and initial conditions"

n = 200

delta = np.zeros(n)
delta[0] = -0.00004

hourglass = np.zeros(n)

chi = np.zeros(n)

position = np.array(range(n))

cutoff = n

#for i in range(1,n):

	#delta[i] = delta[0]

	#for k in range(i):

		#delta[i] *= (np.cos(chi[k-1]) - np.sin(chi[k-1])) / (np.cos(chi[k-1]) + np.sin(chi[k-1]))

	#hourglass[i-1] = delta[i] * ( 2 *np.cos(chi[i-1])) / (np.cos(chi[i-1]) + np.sin(chi[i-1]))
		
	#chi[i] = chi[i-1] + (np.cos(chi[i-1]) + np.sin(chi[i-1]) + 1) / (np.cos(chi[i-1]) + np.sin(chi[i-1])) * delta[i-1]

#hourglass[-1] = delta[n-1] * ( 2 *np.cos(chi[n-1])) / (np.cos(chi[n-1]) + np.sin(chi[n-1]))

for i in range(1,n):

	delta[i] = delta[0]

	for k in range(i):

		delta[i] *= np.tan(np.pi/4 - chi[k-1])

	hourglass[i-1] = delta[i] * ( 2 *np.cos(chi[i-1])) / (np.sin(chi[i-1] + np.pi/4))
		
	chi[i] = chi[i-1] + (1 + 1/(np.sin(chi[i-1] + np.pi/4))) * delta[i-1]

hourglass[-1] = delta[n-1] * ( 2 *np.cos(chi[n-1])) / (np.sin(chi[n-1] + np.pi/4))

plt.plot(position, chi / chi[n-1])
plt.plot(position, hourglass / hourglass[0], 'r')


"Fit with a rational curve"
#guess = (1, 400, 250, 1, 1, 100)
#popt, pcov = curve_fit(rational3_3, position, chi / chi[n-1], p0 = guess)
#print(popt.astype(int))
#plt.plot(position, rational3_3(position, *popt), label='fit', color='y')


"Compare with an exponential decay"
length = 113.1
#plt.plot(position, (1-np.exp(-position/length))/(1-np.exp(-position[n-1]/length)), label='fit', color='y')

plt.show()
