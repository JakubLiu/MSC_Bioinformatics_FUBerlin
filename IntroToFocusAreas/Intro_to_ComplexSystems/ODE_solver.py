import numpy as np
from matplotlib import pyplot as plt

# set values ------------------------------------------------------------------------------------------------------------------
ka = 0.5
ke = 0.3
reaction_rates = [ka, ke]
S = np.array([[-1, 0],
              [1, -1]])
x0_init = 200
x1_init = 0
initial_state_vector = [x0_init, x1_init]
# --------------------------------------------------------------------------------------------------------------------------------------

# function to calculate the propensities (R)
# as an input it takes the state vector and the reaction rates
# it outputs the propesity vector (R)
def propensities(state_vector, reaction_rates):
    R = np.zeros(2)
    R[0] = state_vector[0]*reaction_rates[0]
    R[1] = state_vector[1]*reaction_rates[1]
    return R

# function to calculate the ordinary differential equation (dot product of the stoichiometrix matrix and the reaction rate vector R)
# as an input it takes the stoichiometric matrix,
#, the function to calculate the propensity vector (R),
#, the current state and the reaction rates
# it outputs the ODEs
def ODE(stoi_matrix, prop_fun, state_vector, rates):
   return np.dot(stoi_matrix, prop_fun(state_vector, rates))


# this is the function to approximate the ODE with the explicit Euler method
# as an input it takes the function to calculate the propenisty vector (R),
#, the function to calculate the ODE, aswell as all the arguments of these two functions
# it also takes as an input the initial state vector, the time step and the final time
def ExplicitEuler(init_state_vector, time_step, time_final, prop_fun, ODE_fun, stoi_matrix, reaction_rates):
    NumSteps = int(time_final/time_step) # number of time points
    times = np.zeros(NumSteps+1)  # 1D array to store the times

    # 2D array to store the states of x0, x1 at each given time
    # for every time step there is a row
    # one column for x0 and one column for x1
    states = np.zeros((np.size(times), len(init_state_vector)))
    states[0] = init_state_vector

    for step in range(0, NumSteps):
        times[step+1] = times[step] + time_step
        states[step+1] = states[step] + time_step*ODE_fun(stoi_matrix = stoi_matrix, prop_fun = prop_fun, state_vector = states[step], rates = reaction_rates)
    
    return times, states

(times, states) = ExplicitEuler(init_state_vector=[200,0], time_step=1, time_final=10, prop_fun=propensities, ODE_fun=ODE, stoi_matrix=S, reaction_rates=reaction_rates)


x0_states = []
x1_states = []
for i in states:
    x0_states.append(i[0])
    x1_states.append(i[1])

#Output = np.concatenate((np.array(times,ndmin=2),np.array(x1_states,ndmin=2)), axis=0)
#np.savetxt("C:/Users/Qba Liu/Documents/STUDIA/BIOINF_MASTER_BERLIN/SEMESTER_I/INTRODUCTION_TO_FOCUS_AREAS/COMPLEX_SYSTEMS/LAB/Ex3Task2aTraj.txt",Output,delimiter = ",",fmt="%1.3f")


plt.plot(times, x0_states, label="x0")
plt.plot(times,x1_states, label = "x1")
plt.show()
