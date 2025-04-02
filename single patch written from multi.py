import time
from math import sin
from numpy import array, arange, zeros
from pylab import plot, xlabel, show, legend, yscale, figure, ylabel

threshold = 799
r1 = 0.5
r2 = 0.49
mu = 0.001
d = 0.1
k = 1000

h = 0.001
xpoints =  [] 
ypoints = [] 
xsum = []
ysum = []
timepoints = [0]
iteration = 0
t = 0
current_size = 0

xnext = 1
ynext = 0

while timepoints[-1] < 5000:
    
    xpoints.append(xnext)
    ypoints.append(ynext)
    
    x = xpoints[-1]
    y = ypoints[-1]
    
    
    fx_element = (1.0 - mu) * r1 * x * (1.0 - (x + y) / k) - d * x 
    fy_element = mu * r1 * x * (1.0 - (x + y) / k) + r2 * y * (1.0 - (x + y) / k) - d * y 
    
    xnext = x + h*fx_element 
    ynext = y + h*fy_element

    current_size = xnext + ynext
    
    iteration += 1
    t += h
    timepoints.append(t)

# Plot
figure()
plot(timepoints[:-1], xpoints, label='total number of wildtype')
plot(timepoints[:-1], ypoints, label='total number of mutants')
#yscale('log')
legend()
xlabel("t")
#ylabel("logscale of population abundances")
show()
print(f"The number of mutants at size {threshold} is ", ypoints[-1]) 


'''
threshold = 799
r1 = 0.5
r2 = 0.49
mu = 0.001
d = 0.1
k = 1000# why is k off by a factor of 10, maybe it's not - we have 100 patches, not 10 

h = 0.0001
xpoints =  [] 
ypoints = [] 
xsum = []
ysum = []
timepoints = []
iteration = 0
t = 0
current_size = 0


xnext = 1
ynext = 0

while current_size < threshold:
    print( threshold - current_size)
    
    xpoints.append(xnext)
    ypoints.append(ynext)
    
    x = xpoints[-1]
    y = ypoints[-1]
    
    
    fx_element = (1.0 - mu) * r1 * x * (1.0 - (x + y) / k) - d * x 
    fy_element = mu * r1 * x * (1.0 - (x + y) / k) + r2 * y * (1.0 - (x + y) / k) - d * y 
    
    xnext = x + h*fx_element 
    ynext = y + h*fy_element

    current_size = xnext + ynext
    
    iteration += 1
    t += h
    timepoints.append(t)

# Plot
figure()
plot(timepoints, xpoints, label='total number of wildtype')
plot(timepoints, ypoints, label='total number of mutants')
#yscale('log')
legend()
xlabel("t")
show()
print(f"The number of mutants at size {threshold} is ", ypoints[-1]) 

# CODE
# x = 791.8762411848188
# y = 7.121390458538885

# for 110 patches, this becomes 8  mutants


# ANALYTICALLY 
u = 0.001
xequilibrium = k *(d - r1 + r1 *u) * (-r1 + r2 + r1 *u)   /   (r1 *(-1 + u)*(-r1 + r2 + r1 *u**2))
print("xequilibrium", xequilibrium) # 760
# number of mutants is 40 analytically 
'''