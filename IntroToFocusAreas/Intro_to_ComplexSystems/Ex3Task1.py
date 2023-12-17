import numpy as np
from matplotlib import pyplot as plt

path_inputtxt = "C:/Users/Qba Liu/Documents/STUDIA/BIOINF_MASTER_BERLIN/SEMESTER_I/INTRODUCTION_TO_FOCUS_AREAS/COMPLEX_SYSTEMS/LAB/Input_ex3.txt"
In = np.loadtxt(path_inputtxt)
NrSimulations = int(In[1])
np.random.seed(seed=int(In[0]))

# function to calculate the propensities based on the current state and the reaction constants
def propensities(state_scalar, k1, k2, k3, k4):
    R = np.zeros((4,1))
    R[0] = k1*state_scalar*(state_scalar-1)
    R[1] = k2*state_scalar*(state_scalar-1)*(state_scalar-2)
    R[2] = k3
    R[3] = k4*state_scalar

    return R

# function to calculate the time to the next reaction
def time_to_next_reaction(lambd):
    r = np.random.rand()
    while r == 0:
        r = np.random.rand()
    return (1.0/lambd)*np.log(1.0/r)

def choose_reaction(propensities):
    r = np.random.rand()
    while r == 0:
        r = np.random.rand()
    return np.sum(np.cumsum(propensities) < r*np.sum(propensities))

    

def SSA(S, prop_fun, init_state, simulation_time, k1, k2, k3, k4):
    # S --> stoichioetrix matrix
    # init_state --> values at the start of the simulation
    # simulation time --> how long we run the simulation
    # lam, delta, beta, k_r --> reaction specific parameters

    t = 0.0
    ts = [0.0]  # a list to keep track of the times
    current_state = np.copy(init_state)
    res_all = [current_state] # a list to keep track of the states (for example the concentrations)

    while True:

        prop = prop_fun(current_state, k1, k2, k3, k4)
        prop_sum = np.sum(prop)

        if prop_sum == 0:
            stop = "0 prop"
            break
        
        dt = time_to_next_reaction(prop_sum) # sample random time to next event

        if t + dt > simulation_time:
            stop = "timeout"
            break

        idx = choose_reaction(propensities=prop) # sample random event from event space
        
        change_to_apply = S[idx]
        #change_to_apply.shape = len(change_to_apply)

        # update state vector and time
        t = t + dt
        current_state = current_state + change_to_apply
        ts.append(t)
        res_all.append(current_state)  # counts for all elements
    
    return (ts, np.array(res_all), stop)
    

# stoichiometric matrix
#stoi_matrix = np.array([[1, -1, 1, -1]])
stoi_matrix = [1, -1, 1, -1]

# reaction rates k
k1 = 0.1500
k2 = 0.0015
k3 = 20.0000
k4 = 3.5000

# starting values
initial_state = 40

#Excercise 1a *********************************************************************************************************************
"""

for i in range(0,NrSimulations+1):
    (times, counts, stop) = SSA(S = stoi_matrix, prop_fun=propensities, init_state=initial_state, simulation_time=5, k1 = k1, k2 = k2, k3 = k3, k4 = k4)

    Output = np.concatenate((np.array(times,ndmin=2),np.array(counts,ndmin=2)), axis=0)
    np.savetxt("C:/Users/Qba Liu/Documents/STUDIA/BIOINF_MASTER_BERLIN/SEMESTER_I/INTRODUCTION_TO_FOCUS_AREAS/COMPLEX_SYSTEMS/LAB/Ex3Task1Traj"+str(i+1)+".txt",Output,delimiter = ",",fmt="%1.3f")

"""

# Excercise 1b ***********************************************************************************************************************
"""
final_state = []
for i in range(0,200+1):
    (times, counts, stop) = SSA(S = stoi_matrix, prop_fun=propensities, init_state=initial_state, simulation_time=5, k1 = k1, k2 = k2, k3 = k3, k4 = k4)
    final_state.append(counts[-1])
    print(i)

sum = 0
for i in final_state:
    sum = sum + i
sample_mean = round(sum/len(final_state),2)


fig, ax = plt.subplots()
textstr = "Sample mean: {mean}".format(mean = sample_mean)
ax.hist(final_state, bins = 30)
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)

ax.text(0.5, 0.95, textstr, transform=ax.transAxes, fontsize=10,
        verticalalignment='top', bbox=props)
plt.grid()
plt.xlabel("States of x1 at time t=5")
plt.ylabel("Frequency")
plt.title("Histogram of the state of x1 at time t=5")
plt.show()
"""
