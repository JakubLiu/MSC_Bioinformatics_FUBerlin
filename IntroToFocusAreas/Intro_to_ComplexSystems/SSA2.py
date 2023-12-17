import numpy as np
from matplotlib import pyplot as plt
import random
import statistics

In = np.loadtxt("C:/Users/Qba Liu/Documents/STUDIA/BIOINF_MASTER_BERLIN/SEMESTER_I/INTRODUCTION_TO_FOCUS_AREAS/COMPLEX_SYSTEMS/LAB/Input.txt")
NrSimulations = int(In[1])
np.random.seed(seed=int(In[0]))
#random.seed(int(In[0]))

# function to calculate the propensities based on the current state and the reaction constants
def propensities(state_vector, lam, delta, beta, k_r):
    R = np.zeros((6,1))
    R[0] = lam
    R[1] = state_vector[0]*delta
    R[2] = state_vector[0]*state_vector[1]*beta
    R[3] = delta*state_vector[1]*3e7
    R[4] = state_vector[1]*k_r
    R[5] = state_vector[2]*delta

    return R

# function to calculate the time to the next reaction
def time_to_next_reaction(lambd):
    r = np.random.rand()
    while r == 0:
        r = np.random.rand()
    return (1.0/lambd)*np.log(1.0/r)

def which_reaction(propensities):
    reaction = random.choices(population=range(len(list(propensities))),
                             weights = list(propensities),
                             k = 1)
    return reaction

def choose_reaction(propensities):
    r = np.random.rand()
    while r == 0:
        r = np.random.rand()
    return np.sum(np.cumsum(propensities) < r*np.sum(propensities))

    
    

def SSA(S, prop_fun, init_state, simulation_time, lam, delta, beta, k_r):
    # S --> stoichioetrix matrix
    # init_state --> values at the start of the simulation
    # simulation time --> how long we run the simulation
    # lam, delta, beta, k_r --> reaction specific parameters

    t = 0.0
    ts = [0.0]  # a list to keep track of the times
    current_state = np.copy(init_state)
    res_all = [list(current_state)] # a list to keep track of the states (for example the concentrations)
    #res_x2 = [list(current_state)[1]]   # a list to keep track of the states only of X2
    #res_x4 = [list(current_state)[3]]   # a list to keep track of the states only of X4

    while True:

        prop = prop_fun(current_state, lam, delta, beta, k_r)
        prop_sum = np.sum(prop)

        if prop_sum == 0:
            stop = "0 prop"
            break
        
        dt = time_to_next_reaction(prop_sum) # sample random time to next event

        if t + dt > simulation_time:
            stop = "timeout"
            break

        idx = choose_reaction(propensities=prop) # sample random event from event space
        
        change_to_apply = S[:, idx]
        change_to_apply.shape = len(change_to_apply)

        # update state vector and time
        t = t + dt
        current_state = current_state + change_to_apply
        ts.append(t)
        res_all.append(list(current_state))  # counts for all elements
        #res_x2.append(list(current_state)[1])  # counts for only X2
        #res_x4.append(list(current_state)[3])  # counts for only X4


    
    return (ts, np.array(res_all), stop)
    


    





# stoichiometric matrix
stoi_matrix = np.array([[1, -1, -1, 0, 0, 0],
              [0, 0, 1, -1, -1, 0],
              [0, 0, 0, 0, 1, -1],
              [0, 0, 0, 1, 0, 0]])

# reaction rates k
lam = 1e-4
delta = 1e-8
beta = 5e-5
k_r = 0.3

# starting values
x1_init = lam/delta
x2_init = 2
x3_init = 0
x4_init = 0
initial_state_vector = np.array([x1_init, x2_init, x3_init, x4_init])


x2_counts = []
x4_counts = []

for i in range(1000):
    (times, counts, stop) = SSA(S = stoi_matrix, prop_fun=propensities, init_state=initial_state_vector, simulation_time=10, lam = lam, delta=delta, beta = beta, k_r = k_r)
    x2_counts.append(list(counts[-1])[1])
    x4_counts.append(list(counts[-1])[3])

infection_ongoing = 0
for i in x2_counts:
    if i != 0.0:
        infection_ongoing = infection_ongoing + 1
result = str(infection_ongoing/len(x2_counts)*100) + "%"

fig, ax = plt.subplots()
ax.hist(x2_counts, bins = 30)
textstr = "Probability of ongoing infection\nat time t = 10: {res}".format(res = result)
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
ax.text(0.4, 0.95, textstr, transform=ax.transAxes, fontsize=10,
        verticalalignment='top', bbox=props)
plt.grid()
plt.title("number of infected individuals at time t = 10 for 1000 simulations")
plt.xlabel("number of infected at time t = 10")
plt.ylabel("frequency")
plt.show()


fig, ax = plt.subplots()
avg = str(round(statistics.mean(x4_counts),2))
std = str(round(statistics.stdev(x4_counts),2))
textstr = "Mean: {mean}\nStdev: {sd}".format(mean = avg, sd = std)

ax.hist(x4_counts, bins = 30)
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)

ax.text(0.6, 0.95, textstr, transform=ax.transAxes, fontsize=10,
        verticalalignment='top', bbox=props)
plt.grid()
plt.xlabel("number of dead at time t = 10")
plt.ylabel("frequency")
plt.title("Number of dead individuals at time t = 10 for 1000 simulations")

plt.show()
