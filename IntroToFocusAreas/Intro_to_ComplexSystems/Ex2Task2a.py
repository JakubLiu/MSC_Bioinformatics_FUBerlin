"""
S --> (r1) --> O
S --> (r2) --> I

S --> suscpetible particiapnts that are still taking part in the clinical trail
I --> individuals that got infected durign the clinical trail
O --> individuals that dropped out from the clinical trail

Initial values based on paper:
    S_test: 21720 --> the number of patients that received the vaccine
    S_placebo: 21728 --> number of patients that received the placebo 
    I: 0

r1 --> the drop-out reaction
r2 --> observed infection reaction

k_dropout --> reaction constant of r1
k_inf --> reaction constant of r2
k_inf --> given in file (first value used)
k_drop_out --> can be dervied from k_inf and the averge follow up time (assumed follow up time is 21 days)

propensities:
    r1 = k_dropoout*S
    r2 = k_inf*S
"""

import numpy as np
from matplotlib import pyplot as plt
import statistics

# set random seed
path_inputtxt = "C:/Users/Qba Liu/Documents/STUDIA/BIOINF_MASTER_BERLIN/SEMESTER_I/INTRODUCTION_TO_FOCUS_AREAS/COMPLEX_SYSTEMS/LAB/Input.txt"
In = np.loadtxt(path_inputtxt)
NrSimulations = int(In[1])
np.random.seed(seed=int(In[0]))

# set reaction rate constants
path_kinftxt = "C:/Users/Qba Liu/Documents/STUDIA/BIOINF_MASTER_BERLIN/SEMESTER_I/INTRODUCTION_TO_FOCUS_AREAS/COMPLEX_SYSTEMS/LAB/kinf.txt"
k_inf_file = np.loadtxt(path_kinftxt)
k_inf = k_inf_file[0]   # choose the first value in the kinf file (this is an assumption, probably wrong)
avg_follow_up_days = 21  # assumed average follow up time in days (this is another assumption, probably wrong)
"""
formula for k_dropout:
    k_dropout = (1-avg_follow_up_days*k_inf)/avg_follow_up_days
"""
k_dropout = (1-avg_follow_up_days*k_inf)/avg_follow_up_days

# setting stoichiometry matrix
stoi_matrix = np.array([[-1, -1],
                        [0, 1]])

# setting the initial states
S_init_placebo = 21728
I_init = 0
initial_state_vector = np.array([S_init_placebo, I_init])

# function to calculate the propensities based on the current state and the reaction constants
def propensities(state_vector, k_inf, k_dropout):
    R = np.zeros((2,1))
    R[0] = state_vector[0]*k_dropout
    R[1] = state_vector[0]*k_inf
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


def SSA(S, prop_fun, init_state, simulation_time, k_inf, k_dropout):
    # S --> stoichioetrix matrix
    # init_state --> values at the start of the simulation
    # simulation time --> how long we run the simulation
    # lam, delta, beta, k_r --> reaction specific parameters

    t = 0.0
    ts = [0.0]  # a list to keep track of the times
    current_state = np.copy(init_state)
    res_all = [list(current_state)] # a list to keep track of the states (for example the concentrations)

    while True:

        prop = prop_fun(current_state, k_inf, k_dropout)
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
    
    return (ts, np.array(res_all), stop)



# do 1 simulation and write the results to a file
(times, counts, stop) = (times, counts, stop) = SSA(S = stoi_matrix, prop_fun=propensities, init_state=initial_state_vector, simulation_time=100, k_inf=k_inf, k_dropout=k_dropout)
I_change = []

for count in counts:
    I_change.append(count[-1])

Output = np.concatenate((np.array(times,ndmin=2),np.array(I_change,ndmin=2)), axis=0)
np.savetxt("C:/Users/Qba Liu/Documents/STUDIA/BIOINF_MASTER_BERLIN/SEMESTER_I/INTRODUCTION_TO_FOCUS_AREAS/COMPLEX_SYSTEMS/LAB/Task2Infected.txt",Output,delimiter = ",",fmt="%1.3f")