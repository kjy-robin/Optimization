import matplotlib
import numpy as np
import matplotlib.pyplot as plt
from apm import *

server = 'http://byu.apmonitor.com'

app = 'nlp'
apm(server,app,'clear all')
apm_load(server,app,app + '.apm')
apm_option(server,app,'nlc.solver',3)
solver_output = apm(server,app,'solve')
print solver_output
url = apm_web_var(server,app)

## Generate a contour plot
x = np.arange(-13, 13, 0.1)
y = np.arange(-6.0, 1.1, 0.1)
x1, x2 = np.meshgrid(x, y)

# objective function:
obj = np.power(x1,2) - 2.0*x1*x2 + 4.0 * np.power(x2,2)

# constraint1: 
c1 = 0.1 * x1 - x2

n = 50
mu = np.logspace(2,-2,n)

for i in range(0,n):
    # barrier problem
    psi = obj - mu[i] * np.log(c1-1.0)

    # create contour plot
    plt.figure()
    #plt.set_size_inches(3.5,2.5)
    plt.title('Contour Plot with log10(mu) = ' + str(float(int(10.0*np.log10(mu[i])))/10.0))
    lines = [2,5,10,15,20,30,40,50,60,70,80,90,100]
    CS = plt.contour(x1,x2,obj,lines)
    plt.clabel(CS,inline=1,fontsize=10)

    # add x and y plot labels
    plt.xlabel('x1')
    plt.ylabel('x2')

    # constraint #1
    CS = plt.contour(x1,x2,c1, \
    	[1.0], \
	colors='r', \
	linewidths=[4.0])
    plt.clabel(CS, inline=1, fontsize=10)

    # barrier problem
    psi_min = np.floor(np.nanmin(psi))
    #lines = np.linspace(psi_min,100,10)
    CS = plt.contour(x1,x2,psi,lines)
    plt.clabel(CS,inline=1,fontsize=10)

    # save figures
    if i<=9:
        plt.savefig('contour0' + str(i) + '.png')
    else:
        plt.savefig('contour' + str(i) + '.png')
    
# show figures
#plt.show()
