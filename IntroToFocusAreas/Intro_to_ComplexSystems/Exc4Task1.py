import numpy as np
import scipy
import math
from matplotlib import pyplot as plt

# stoichiometrix matrix ___________________________________________________________________________________________________________
S = np.array([[-1, 0],
              [1, -1]])


# ODE function for solve_ivp()_____________________________________________________________________________________________________
def ODE(t,state):

  rates = [0.5, 0.3]
  S = np.array([[-1, 0],
       [1, -1]])

  def prop(state, rates):
    R = np.zeros(2)
    R[0] = state[0]*rates[0]
    R[1] = state[1]*rates[1]
    return R
  
  return np.dot(S, prop(state, rates))



# exact ODE function ____________________________________________________________________________________________________
def exact_ODE(t, reaction_rates, init):
    ka = reaction_rates[0]
    ke = reaction_rates[1]
    e1 = (-1)*ke*t
    e2 = (-1)*ka*t
    return init * (ka/(ka-ke)) * (math.exp(e1) - math.exp(e2))



# solve_ivp() usage __________________________________________________________________________________________________________________
initial_state = np.array([[200],[0]])
sol = scipy.integrate.solve_ivp(ODE, [1.00,24.00], initial_state.flatten('F'), method = 'RK45')
# sol.y --> states
# sol.t --> times

# generating exact trajectory________________________________________________________________________________________________________
exact_states = []
for time in sol.t:
   state = exact_ODE(t=time, reaction_rates=[0.5, 0.3], init = 200)
   exact_states.append(state)


RK45_states = sol.y[1,:]
plt.plot(sol.t, RK45_states, label = "RK45 approximation")
plt.plot(sol.t, exact_states, label = "exact solution")
plt.legend()
plt.xlabel("times")
plt.ylabel("states")
plt.grid()
plt.show()