'''
Goal of the simulation:
stochastically simulate the ODEs

dx/dt = rx(1-(x+y)/k)
= rx     -     rx(x+y)/k)
dy/dt = urx(1-(x+y)/k)
= urx    -    urx(x+y)/k
using a Gillespie algorithm
'''
# import statements
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import random
import sys

# initialize lists
X = [1] 
t = [0]
Y = [0]

# end of simulation
tend = 100 

# parameters
k = 10**5
r1 = 0.5
r2=0.49
u = 0#10**-5
d = 0.1


# HOW COME I KEEP GETTING A DIVISION BY ZERO ERROR? Sometimes rate_sum = 0 - WHY?

# loop
while t[-1] < tend:
        
        # define the current x value ------------------------------------------------------|
        current_X = X[-1]
        current_Y = Y[-1]
        
        # define the rates of the reaction -----------------------------------------------|
        rates = [r1*current_X*(1 - ((current_X+current_Y)/k)),                      (d*current_X),                            u*r1*current_X*(1 - ((current_X+current_Y)/k))+(r2*current_Y)*(1 - ((current_X+current_Y)/k)),                                                                  (d*current_Y)]

        # find the total sum of the rates-------------------------------------------------|
        rate_sum = sum(rates) # propensity that any event will occur
        if rate_sum == 0:
                print("Population crashed")
                sys.exit(0)
        # randomly pick the time untill the next event ------------------------------|
        tau = np.random.exponential(scale=1/rate_sum)
        # add to time list
        t.append(t[-1] + tau)
        
        # randomly pick which event happens next----------------------------------|
        # pick a random number
        rand = random.uniform(0,1)
        
        # make sure the carrying capacity is not exceeded -----------------------|
        # need to make (1 - ((current_X+current_Y)/k) ) probability 0 if x+y > k so it's not negative, let another event occur in that time step
        if (current_X+current_Y)>k:
                rates[0] = 0
                rates[2] = 0

        # production event -----------------------------------------------------------------|
        # if this random number is between 0 and the first event
        if rand * rate_sum > 0 and rand * rate_sum <= rates[0]: # in the case of current_X+current_Y)>k .... rates[0] = 0 so the condition becomes
                # if rand*rate_sum > 0 and rand * rate_sum <0 which will never happen
                # but then do we need to reassign the probabilities?
                # no because the next one will be between 0 and then next probability 
                
                # the first event occurs
                # population increases by 1
                X.append(X[-1] + 1)
                Y.append(Y[-1])

        # decay event -----------------------------------------------------------------------|
        # if this random number is instead between the first rate and then second rate
        elif rand * rate_sum > rates[0] and rand * rate_sum <= rates[0] + rates[1]:
                # the second event occurs
                # population decreases by 1
                X.append(X[-1] - 1)
                Y.append(Y[-1])
        
        # producrtion event for y --------------------------------------------------------|
        elif rates[0] + rates[1] < rand * rate_sum and rand * rate_sum <= rates[0] + rates[1] + rates[2] :
                X.append(X[-1])
                Y.append(Y[-1]+1)
                
        # decay event for y  --------------------------------------------------------------|
        elif rates[0] + rates[1]+rates[2] < rand * rate_sum and rand * rate_sum <= rates[0] + rates[1] + rates[2]+rates[3]:
                X.append(X[-1])
                Y.append(Y[-1]-1)
                
equilibria_list = []
equilibria = 80000
for i in range(len(t)):
        equilibria_list.append(equilibria)
# plot
plt.plot(t,X)
plt.plot(t,Y)
plt.plot(t,equilibria_list)
plt.xlabel("time")
plt.ylabel(" Wildtype Population ")
plt.show()


'''

# import statements
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import random

# initialize lists
X = [1] 
t = [0]
Y = [0]

# end of simulation
tend = 100 

# parameters
k = 10**5
r1 = 0.5
r2=0.49
u = 10**-5
d = 0.1

# loop
while t[-1] < tend:
        
        # define the current x value
        current_X = X[-1]
        current_Y = Y[-1]
        
        # define the rates of the reaction
        rates = [r1*current_X*(1 - ((current_X+current_Y)/k)),       (d*current_X),                            u*r1*current_X*(1 - ((current_X+current_Y)/k))+(r2*current_Y)*(1 - ((current_X+current_Y)/k)),              (d*current_Y)]
        
        # need to make (1 - ((current_X+current_Y)/k) ) probability 0 if x+y > k so it's not negative, let another event occur in that time step
        
        
        
        # find the total sum of the rates
        rate_sum = sum(rates) # propensity that any event will occur
        
        # randomly pick the time untill the next event ------------------------------|
        tau = np.random.exponential(scale=1/rate_sum)
        # add to time list
        t.append(t[-1] + tau)
        
        # randomly pick which event happens next----------------------------------|
        # pick a random number
        rand = random.uniform(0,1)

        # production event
        # if this random number is between 0 and the first event
        if rand * rate_sum > 0 and rand * rate_sum <= rates[0]:
                # the first event occurs
                # population increases by 1
                X.append(X[-1] + 1)
                Y.append(Y[-1])

        # decay event
        # if this random number is instead between the first rate and then second rate
        elif rand * rate_sum > rates[0] and rand * rate_sum <= rates[0] + rates[1]:
                # the second event occurs
                # population decreases by 1
                X.append(X[-1] - 1)
                Y.append(Y[-1])
        
        # producrtion event for y
        elif rates[0] + rates[1] < rand * rate_sum and rand * rate_sum <= rates[0] + rates[1] + rates[2] :
                X.append(X[-1])
                Y.append(Y[-1]+1)
                
        # decay event for y        
        elif rates[0] + rates[1]+rates[2] < rand * rate_sum and rand * rate_sum <= rates[0] + rates[1] + rates[2]+rates[3]:
                X.append(X[-1])
                Y.append(Y[-1]-1)

# plot
plt.plot(t,X)
plt.plot(t,Y)
plt.xlabel("time")
plt.ylabel(" Wildtype Population ")
plt.show()

'''