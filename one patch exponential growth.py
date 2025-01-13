# Deme Model - one deme
# simulating ODE's

# TO DO - TURN INTO A WHILE LOOP INSTEAD - JUST LIKE MULTIDEME ONE?

# so far best version - ONE DEME VERSION

# apply text book to my model and then compare with chatgtp dominik's code

# fourth order Runge-Kutta

# equations
# dx/dt = r1 * x * (1 - mu) * (1 - ((x + y) / k)) - d * x
# dy/dt = mu * r1 * x * (1 - mu) * (1 - ((x + y) / k)) + r2 * y * (1 - ((x + y) / k)) - d * y

threshold = 800000#795000#10**4

# define function parameters
r1 = 0.5  # 0.1
r2 = 0.49
mu = 10**-5
d = 0.1  # 0.01
#if you go to 1E5 cells, the mutation rate should be round 1E-5 or 2E-5!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
k = 10**5  # for comparison to 10 deme model, equilibrium will be 800, run untill total pop is 800
# equilibrium 76
# import statements
from math import sin
from numpy import array, arange
from pylab import plot, xlabel, show, legend, yscale, xscale

# define funtion
def f(r, t):  # r is a vector with the first element as the first function and the second element as the second function
    x = r[0]
    y = r[1]
    fx = r1 * x * (1 - mu)  - d * x  # never calls t
    fy = mu * r1 * x + r2 * y - d * y
    return array([fx, fy], float)



# define parameters
a = 0  # start of interval
b = 10000  # 10000  # 10000  # 10 end of interval
# N =  # 1000000  # number of steps
h = 0.01  # (b - a) / N  # step size # 0.01

# iterate through
xpoints = []
ypoints = []
tpoints=[]
r = array([1, 0], float)  # x and y start at 1 at time 0. so at time 0, there is one wildtype cell and 1 mutant cell
# print(r)
#t   = 2
#print('f=', f(r, t))
t=0
current_size = 0
while current_size < threshold:
    

    xpoints.append(r[0])
    ypoints.append(r[1])
    k1 = h * f(r, t)            # remember r is a vector, so is f and so is k1,k2,k3,k4
    k2 = h * f(r + 0.5 * k1, t + 0.5 * h)
    k3 = h * f(r + 0.5 * k2, t + 0.5 * h)
    k4 = h * f(r + k3, t + h)
    r += (k1 + 2 * k2 + 2 * k3 + k4) / 6
    
    # also - problem there is no t in my function
    
    # if threshold population size is reached then stop simulation
    t+= h 
    tpoints.append(t)
    current_size = (xpoints[-1] + ypoints[-1])
    print(current_size)

# plot
plot(tpoints, xpoints, label='x(t) number of wildtype cells')
plot(tpoints, ypoints, label='y(t) number of mutant cells')
#yscale('log')
xlabel("t")
legend()
show()

'''
equilibrium values
xpoints[-10]
75.98478021996047
ypoints[-1]
3.9951997600176385
'''


