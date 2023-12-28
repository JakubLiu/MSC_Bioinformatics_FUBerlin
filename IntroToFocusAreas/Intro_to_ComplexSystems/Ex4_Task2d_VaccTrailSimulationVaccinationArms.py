import numpy as np
from matplotlib import pyplot as plt
import statistics

# function definitons ===============================================================================================================

# function to calculate the propensities based on the current state and the reaction constants
# k_inf --> reaction rate constant
# k_dr --> drop out rate constant
def propensities(state_vector, k_inf, k_dr):
    R = np.zeros((2,1))
    R[0] = state_vector[0]*k_dr
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

    

def SSA(S, prop_fun, init_state, simulation_time, k_inf, k_dr):

    t = 0.0
    ts = [0.0]  # a list to keep track of the times
    current_state = np.copy(init_state)
    res_all = [list(current_state)] # a list to keep track of the states (for example the concentrations)

    while True:

        prop = prop_fun(current_state, k_inf, k_dr)
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



def SI_simulate(vacc_name, S, S_init, I_init, avg_follow_up_time, k_inf_file_path, actual_infected_num):

    initial_state_vector = np.array([S_init, I_init])

    # reaction rates k
    #k_inf = k_inf will be set in a loop
    T_hat = avg_follow_up_time / S_init  # average follow up per person
    #k_dr =  this will also be set inside a loop

    # opening one of the k_inf file
    file = open(k_inf_file_path, 'r')
    k_inf_all = file.readlines()
    final_state_I = []
    count = 0

    for k_inf in k_inf_all:  # loop over the k_inf values in the k_inf file
        while count <= 100:  # this specifies the number of simulations
            k_dr = 1/T_hat - float(k_inf) # calculate drop out rate based on formula (from assignment 2)
            (times, counts, stop) = SSA(S = S, prop_fun=propensities, init_state=initial_state_vector, simulation_time=10, k_inf = float(k_inf), k_dr = k_dr)
            final_state_I.append(counts[-1][1])   # final state
            count = count + 1

    stdev = round(statistics.stdev(final_state_I),2)
    fig, ax = plt.subplots()
    plt.hist(final_state_I)
    plt.axvline(x = actual_infected_num, color = "red")
    plt.grid()
    plt.title(vacc_name + "(vaccination arm)")
    plt.xlabel("Number of infected at end time")
    plt.ylabel("Frequency")
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    textstr = "-Red line shows actual value from study\n-Standard deviation (uncertainity): {}"
    ax.text(0.05, 0.95, textstr.format(stdev), transform=ax.transAxes, fontsize=8,verticalalignment='top', bbox=props)
    plt.show()

    return(final_state_I)
    





# study parameter definitions ==============================================================================================================

# Biontech
vacc_name_biontech = "Biontech"
S_biontech = np.array([[-1, -1],
            [0,   1]])
S_init_biontech = 18860.00
I_init_biontech = 0.00
avg_follow_up_time_biontech = 2.214 * 1000.00
k_inf_file_path_biontech = "C:/Users/Qba Liu/Documents/STUDIA/BIOINF_MASTER_BERLIN/SEMESTER_I/INTRODUCTION_TO_FOCUS_AREAS/COMPLEX_SYSTEMS/LAB/kinf_Biontech_Vacc.txt"
actual_infected_num_biontech = 8.00

# Moderna
vacc_name_moderna = "Moderna"
S_moderna = np.array([[-1, -1],
            [0,   1]])
S_init_moderna = 14134.00
I_init_moderna = 0.00
avg_follow_up_time_moderna = 3.39 * 1000.00
k_inf_file_path_moderna = "C:/Users/Qba Liu/Documents/STUDIA/BIOINF_MASTER_BERLIN/SEMESTER_I/INTRODUCTION_TO_FOCUS_AREAS/COMPLEX_SYSTEMS/LAB/kinf_Moderna_Vacc.txt"
actual_infected_num_moderna = 19

# Astra
vacc_name_astra = "Astra Zeneca"
S_astra = np.array([[-1, -1],
            [0,   1]])
S_init_astra = 21587.00
I_init_astra = 0.00
avg_follow_up_time_astra = 64.8/365.00 * 1000.00
k_inf_file_path_astra = "C:/Users/Qba Liu/Documents/STUDIA/BIOINF_MASTER_BERLIN/SEMESTER_I/INTRODUCTION_TO_FOCUS_AREAS/COMPLEX_SYSTEMS/LAB/kinf_Astra_Vacc.txt"
actual_infected_num_astra = 73





# executiuon ===========================================================================================================================

vaccination = "Astra"

if vaccination == "Biontech":
    print(SI_simulate(vacc_name_biontech, S_biontech, S_init_biontech, I_init_biontech, avg_follow_up_time_biontech, k_inf_file_path_biontech, actual_infected_num_biontech))

elif vaccination == "Moderna":
    print(SI_simulate(vacc_name_moderna, S_moderna, S_init_moderna, I_init_moderna, avg_follow_up_time_moderna, k_inf_file_path_moderna, actual_infected_num_moderna))

else:  # vaccination == "Astra"
    print(SI_simulate(vacc_name_astra, S_astra, S_init_astra, I_init_astra, avg_follow_up_time_astra, k_inf_file_path_astra, actual_infected_num_astra))