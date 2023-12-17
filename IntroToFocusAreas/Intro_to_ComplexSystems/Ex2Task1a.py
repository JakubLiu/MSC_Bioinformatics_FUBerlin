import numpy as np
from matplotlib import pyplot as plt

path_inputtxt = "C:/Users/Qba Liu/Documents/STUDIA/BIOINF_MASTER_BERLIN/SEMESTER_I/INTRODUCTION_TO_FOCUS_AREAS/COMPLEX_SYSTEMS/LAB/Input.txt"
In = np.loadtxt(path_inputtxt)
NrSimulations = int(In[1])
np.random.seed(seed=int(In[0]))

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
x2_init = 20
x3_init = 0
x4_init = 0
initial_state_vector = np.array([x1_init, x2_init, x3_init, x4_init])


for i in range(0,NrSimulations+1):
    (times, counts, stop) = SSA(S = stoi_matrix, prop_fun=propensities, init_state=initial_state_vector, simulation_time=10, lam = lam, delta=delta, beta = beta, k_r = k_r)

    x1_change = []
    x2_change = []
    x3_change = []
    x4_change = []

    for count in counts:
        x1_change.append(count[0])
        x2_change.append(count[1])
        x3_change.append(count[2])
        x4_change.append(count[3])

    Output = np.concatenate((np.array(times,ndmin=2),np.array(x2_change,ndmin=2)), axis=0)
    np.savetxt("C:/Users/Qba Liu/Documents/STUDIA/BIOINF_MASTER_BERLIN/SEMESTER_I/INTRODUCTION_TO_FOCUS_AREAS/COMPLEX_SYSTEMS/LAB/traj"+str(i+1)+".txt",Output,delimiter = ",",fmt="%1.3f")


