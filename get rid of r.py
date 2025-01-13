
# currently working on 

import time
from math import sin
from numpy import array, arange, zeros
from pylab import plot, xlabel, show, legend, yscale, figure

threshold = 800000
r1 = 0.5
r2 = 0.49
mu = 10**-5
d = 0.1
k = 10**4
eta = 0.01
h = 0.1
Patches = 110

xpoints =  [[0] for _ in range(Patches)] 
xpoints[0] = [1]# add a wildtype to the first patch
# xpoints = [[1], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0]]

ypoints = [[0] for _ in range(Patches)] 
# ypoints = [[0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0]]

xsum = []
ysum = []
timepoints = []
iteration = 0
t = 0
current_size = 0
xtemp = zeros(Patches)
# xtemp = array([0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.])

ytemp = zeros(Patches)
# xtemp = array([0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.])

while current_size < threshold:
    
    xs = 0 # each iteration start with a brand new sum for summing across populations 
    ys = 0
    
    for n in range(Patches):
        
        # new iteration
        x = xpoints[n][-1] # current x value
        y = ypoints[n][-1] # current y value

        sum_x = sum_y = 0.0
        for j in range(Patches): # sum across all other patches except current one for migration term
            if j != n:
                sum_x += xpoints[j][-1]  
                sum_y += ypoints[j][-1]

        fx_element = (1.0 - mu) * r1 * x * (1.0 - (x + y) / k) - d * x - eta * (x - sum_x / (Patches - 1))
        fy_element = mu * r1 * x * (1.0 - (x + y) / k) + r2 * y * (1.0 - (x + y) / k) - d * y - eta * (y - sum_y / (Patches - 1))
        
        xtemp[n] = x + h*fx_element # store xnext as temporary variable so i can update all patches simultaneously at the end
        ytemp[n] = y + h*fy_element

        xs += xpoints[n][iteration] # im adding a patch per iteration to the total sum
        ys += ypoints[n][iteration]
        
    for n in range(Patches): # adding xnext and ynext to the list of timepoints
        xpoints[n].append(xtemp[n])
        ypoints[n].append(ytemp[n])        
    
    current_size = xs + ys
    #print(current_size)

    xsum.append(xs)
    ysum.append(ys)
    iteration += 1
    t += h
    timepoints.append(t)

# Plot
figure()
plot(timepoints, xsum, label='total number of wildtype')
plot(timepoints, ysum, label='total number of mutants')
#yscale('log')
legend()
xlabel("t")
show()
print(f"The number of mutants at size {threshold} is ", ysum[-1])
