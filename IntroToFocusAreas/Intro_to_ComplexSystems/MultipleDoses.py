import numpy as np
from matplotlib import pyplot as plt


def propensities(state_vector, reaction_rates):
    R = np.zeros(2)
    R[0] = state_vector[0]*reaction_rates[0]
    R[1] = state_vector[1]*reaction_rates[1]
    return R

def ODE(stoi_matrix, prop_fun, state_vector, rates):
   return np.dot(stoi_matrix, prop_fun(state_vector, rates))


def ExplicitEuler(init_state_vector, time_step, time_final, prop_fun, ODE_fun, stoi_matrix, reaction_rates):
    NumSteps = int(time_final/time_step) # number of time points
    times = np.zeros(NumSteps+1)  # 1D array to store the times

    states = np.zeros((np.size(times), len(init_state_vector)))
    states[0] = init_state_vector

    for step in range(0, NumSteps):
        times[step+1] = times[step] + time_step
        states[step+1] = states[step] + time_step*ODE_fun(stoi_matrix = stoi_matrix, prop_fun = prop_fun, state_vector = states[step], rates = reaction_rates)
    
    return times, states


S = np.array([[-1, 0],
              [1, -1]])


NumDoses = 28  # 1 dose every 6h for 7 days
times_total = []
states_total = []
state = [200, 0]
simulation_time = 6.00 # because of the given dosing frequency we have to simulate for 6 hours

for dose in range(1, NumDoses):
  (times, states) = ExplicitEuler(init_state_vector=state, time_step=0.01, time_final=simulation_time, prop_fun=propensities, ODE_fun=ODE, stoi_matrix=S, reaction_rates=[0.5, 0.3])
  times_total.append(times)
  states_total.append(states)
  x0_final = states[-1][0]
  x1_final = states[-1][1]
  state = [200+x0_final, 0+x1_final]

x0_states = []
x1_states = []

for row in states_total:
   for i in row:
      x0_states.append(i[0])
      x1_states.append(i[1])


change_x0 = []
for i in range(0, len(x0_states)-1):
   change_x0.append([i, abs(x0_states[i]-x0_states[i+1])])

change_x1 = []
for i in range(0, len(x1_states)-1):
   change_x1.append([i, abs(x1_states[i]-x1_states[i+1])])


add = 0.00
times_cumulative = []

for row in times_total:
   for i in row:
      times_cumulative.append(i+add)
   add = add + simulation_time
   
plt.plot(times_cumulative, x1_states, label = "x1")
plt.plot(times_cumulative, x0_states, label = "x0")
plt.legend()
plt.grid()
plt.show()


def mean(x):
   sum = 0.00
   for i in x:
      sum = sum + i
   return sum/len(x)

x1_steady = x1_states[int(len(x1_states)/2):]
max_x1_steady = round(max(x1_steady),4)
min_x1_steady = round(min(x1_steady),4)
mean_x1_steady = round(mean(x1_steady),4)
text_x1 = "Mean of x1 steady state: {var1}\nMaximum of x1 steady state: {var2}\nMinimum of x1 steady state: {var3}"
print(text_x1.format(var1 = mean_x1_steady, var2 = max_x1_steady, var3 = min_x1_steady))
print("\n")
x0_steady = x0_states[int(len(x0_states)/2):]
max_x0_steady = round(max(x0_steady),4)
min_x0_steady = round(min(x0_steady),4)
mean_x0_steady = round(mean(x0_steady),4)
text_x0 = "Mean of x0 steady state: {var1}\nMaximum of x0 steady state: {var2}\nMinimum of x0 steady state: {var3}"
print(text_x0.format(var1 = mean_x0_steady, var2 = max_x0_steady, var3 = min_x0_steady))

"""
-After how many doses does the concentrations in the bloodstream stabilize? (this is called
"pharmacokinetic steady state")

  When we add 1 dose every 12h for 7 days with a simulation time of 12h the steady state is never reached. 
  This is because a simulation time of 12h is so long, that the concentrations of x0 and x1 fall to 0.
  Since the concentrations of x0 and x1 are 0, then multiple doses have no cumulative effect.
  If we increase the dosing frequency to 1 dose every 6h (and offcourse decrease the simulation time to 6h) we achieve
  a pharmacokinetic steady state after approximately 15h.


-What are the min (trough), max and average concentration levels in the blood at the "pharmacokinetic steady state"?

   - min concentration: 72.7053
   - max concentration: 135.3499
   - mean concentration: 110.6269
"""