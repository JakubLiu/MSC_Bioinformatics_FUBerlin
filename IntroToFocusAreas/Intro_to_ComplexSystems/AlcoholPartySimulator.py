import numpy as np
from matplotlib import pyplot as plt

# set values ------------------------------------------------------------------------------------------------------------------
ka = 0.150
ke = 0.154
reaction_rates = [ka, ke]
S = np.array([[-1, 0],
              [1, -1]])
x0_init = 200
x1_init = 0
initial_state_vector = [x0_init, x1_init]
# --------------------------------------------------------------------------------------------------------------------------------------

# necessary function definitions ------------------------------------------------------------------------------------------------------
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

# ---------------------------------------------------------------------------------------------------------------------------------------

# define how much mg of ethanol does a drink have
long_drink = 16
wine_glass = 24
beer = 25

# define what and in which order you drink
water = 0 # for coding reasons, the last intake has to have 0mg of ethanol
drinking_sequence = [long_drink*2, wine_glass*2, beer, beer*2, water]

# define how often have you drunk these drinks
time_from_last_drink_to_party_end = 80  # the last time interval corresponds to the time from the last alcohol intake to the leave
drinking_intervals = [10,10,20,time_from_last_drink_to_party_end]


# set lists to hold the values and set initial states
state = [0,0]
times_total = []
states_total = []
state = [drinking_sequence[0], 0]

for what, time_interval in zip(drinking_sequence, drinking_intervals):
    (times, states) = ExplicitEuler(init_state_vector = state, time_step=0.01, time_final=time_interval, prop_fun=propensities, ODE_fun=ODE, stoi_matrix=S, reaction_rates=reaction_rates)
    
    alcohol_conc_in_blood = []
    alcohol_conc_in_digestive_tract = []
    
    for row in states:
        alcohol_conc_in_blood.append(row[1])
        alcohol_conc_in_digestive_tract.append(row[0])

    states_total.append(alcohol_conc_in_blood)   # nested list will get created, deal with this later
    times_total.append(times)
    state = [alcohol_conc_in_digestive_tract[-1] + what, alcohol_conc_in_blood[-1]]


# unlist the nested states_total list
states_total = [state for states_nested in states_total for state in states_nested]


add = 0.00
times_cumulative = []

for row, time_interval in zip(times_total, drinking_intervals):
   for i in row:
      times_cumulative.append(i+add)
   add = add + time_interval


weight = 74
constant = (520*weight)
promiles = []

for state in states_total:
    promiles.append((state/constant)*1000)


# just for us to help
times_bicycle = []
times_car = []
threshold = 0.0009

for conc, time in zip(promiles, times_cumulative):
    if abs(conc - 0.5) <= threshold:
        times_bicycle.append(time)
    elif abs(conc - 1.6) <= threshold:
        times_car.append(time)

print("Times when the promiles are equal to the bicycle threshold: {x}".format(x = times_bicycle))
print("Times when the promiles are equal to the car threshold: {y}".format(y = times_car))

plt.subplot(1,3,1)
plt.plot(times_cumulative, states_total, label = "Alcohol concentration in the bloodstream")
plt.xlabel("time")
plt.ylabel("alcohol concentration in bloodstream")
plt.legend()
plt.grid()

plt.subplot(1,3,2)
plt.plot(times_cumulative, promiles, label = "promiles", color = "green")
plt.axhline(y=0.5, color='blue', linestyle='-', label = "car driving")
plt.axhline(y=1.5, color='r', linestyle='-', label = "bicycle driving")
plt.xlabel("time")
plt.ylabel("promiles")
plt.legend()


plt.grid()
plt.show()
