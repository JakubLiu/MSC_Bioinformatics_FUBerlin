import numpy as np
from matplotlib import pyplot
from itertools import chain
import collections
from functools import reduce
from operator import concat

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



def SI_simulate(StoiMatrix, S_init, I_init, avg_follow_up_time, k_inf, NumSim):

    initial_state_vector = np.array([S_init, I_init])

    T_hat = avg_follow_up_time / S_init  # average follow up per person
    final_state_I = []
    count = 0
    iter = 1

    while count <= NumSim:  # this specifies the number of simulations
        print("---->",iter/NumSim*100, "%")
        k_dr = 1/T_hat - float(k_inf) # calculate drop out rate based on formula (from assignment 2)
        (times, counts, stop) = SSA(S = StoiMatrix, prop_fun=propensities, init_state=initial_state_vector, simulation_time=10, k_inf = float(k_inf), k_dr = 1/T_hat - float(k_inf))
        final_state_I.append(counts[-1][1])   # final state
        count = count + 1
        iter = iter + 1

    return(final_state_I)


"""

MAKE THIS INTO A FUNCTION SO THAT THIS FUNCTION CAN BE USED IN THE OTHER EXERCISE, WHERE THE UNCERTAINITY MUST BE CALCULATED.

P(phi) = sum( P(Itend) + P(phi | Itend) )


A.) create a phi vector [0, 0.01, . . . , 0.99]

    P(phi | Itend)
    1.) use the created phi vector [0, 0.01, . . . , 0.99] (100 phi values)
    2.) for i in phi vector calculate k_inf_vacc = (1 - i) * k_inf_placebo
    3.) for every value of k_inf_vacc simulate the vaccination arm 100 times (so 100x100 total simulations)
    P(Itend)
    4.) for each simulation save the number of infected individuals at the end time (Itend)
    5.) for each Itend calculate its probability and save this probability somewhere
        - P(Itend) = number of simulations that gave that Itend / total number of simulations**
            ** the total number of simulations is 100x100
    6.)* calculate how many phi values have generated a given Itend value and calculate the probability

    * for example:
        - we have 100 phi values
        - from these 100 phi values we obtain 100 k_inf_vacc values
        - for each of the 100 k_inf_vacc values we do 100 simulations (10 000 total simulations)
        - for each of the 10 000 simulations we obtain a Itend
        - the number of distinct Itend values is less then the number of simulations...
          ... and each 100 simulations have a corresponding k_inf_vacc value that each have a corresponding phi value
        - so lets say we consider the 100 simulations done for phi = 0.45
        - 2 simulations with phi = 0.45 yielded a Itend of 16
        - over the total 100 000 simulations, 120 simulations yielded a Itend of 16
        - so P(phi | Itend) = 2/120

    7.) - for each distinct phi value we ran 100 simulations that resulted in < 100 distinct Itend values
        - so for each phi we have a vector of Itend values (length is not known a priori)
        - for each distinct Itend value in the Itend vector calculate the probability of this Itend and store it (vector X)
        - for each distinct Itend value in the Itend vector calculate P(phi | Itend) (as discussed above) and store it (vector Y)
        - for i,j in vector X, vector Y calculate ixj and store it
        - sum it and this is our P(phi) for one of the 100 phi values
        - so, now we need to repeat this for the other 99 phi values
"""

import numpy as np

def Uncertainity(file_path_in):
    phi_vector = np.linspace(0,1,100) #1.)
    file = open(file_path_in, 'r')
    k_inf_placebo = file.readlines()
    k_inf_vacc = []
    NumberOfSimulations = 100

    for i,j in zip(phi_vector, k_inf_placebo):  #2.)
        k_inf_vacc.append((1 - i) * float(j))


    Itend_only = []  # this list will hold only the Itend values

    # 3.), 4.)
    SimResSorted = np.zeros((NumberOfSimulations+1))
    counter = 1
    for kinf in k_inf_vacc:
        print(counter, " out of ", len(k_inf_vacc))
        res = np.array(SI_simulate(StoiMatrix = np.array([[-1, -1], [0,   1]]), S_init = 18860.00, I_init = 0.00, avg_follow_up_time = 2.214 * 1000.00, k_inf=kinf, NumSim=NumberOfSimulations))
        Itend_only.append(res)
        SimResSorted = np.vstack([SimResSorted, res])
        counter = counter + 1

    SimResSorted = np.delete(SimResSorted, (0), axis=0)  # delete redundant first row
    Itend_only = SimResSorted.tolist()
    Itend_only = np.array(reduce(concat, Itend_only))
    # 5.)
    # calculate P(Itend)
    # this is a nested list
    # the first element of the nested list is the Itend value
    # the second element of the nested list is it's probability
    Prop_Itend = {}
    for key, value in collections.Counter(Itend_only).items():
        Prop_Itend[key] = value/len(Itend_only)

    # calculate P(phi | Itend)
    SimResSorted_NumRows, SimResSorted_NumCols = SimResSorted.shape
    unique, counts = np.unique(SimResSorted, return_counts=True)
    CountsOverall = dict(zip(unique, counts))  # dictionary of counts of all Itend values

  
    PhiProb = np.zeros((1,NumberOfSimulations+1))
    Itend_identifier_array = np.zeros((1,NumberOfSimulations+1))
    for phi in range(0, SimResSorted_NumRows):
        unique, counts = np.unique(SimResSorted[phi], return_counts=True)
        CountsSmall = dict(zip(unique, counts)) # dictionary of counts of Itend values for a given phi value
        bucket_list = []
        Itend_identifier = []
        for Itend in CountsSmall.keys():
            bucket_list.append(CountsSmall[Itend]/CountsOverall[Itend])
            Itend_identifier.append(Itend)
        bucket_list = bucket_list + [0.00] * (NumberOfSimulations+1-len(bucket_list))
        Itend_identifier = Itend_identifier + [0.00] * (NumberOfSimulations+1-len(Itend_identifier))
        PhiProb = np.vstack([PhiProb, np.array(bucket_list)])
        Itend_identifier_array = np.vstack([Itend_identifier_array, np.array(Itend_identifier)])

    PhiProb = np.delete(PhiProb, (0), axis=0)  # delete redundant first row
    Itend_identifier_array = np.delete(Itend_identifier_array, (0), axis=0)  # delete redundant first row

    Uncertainity_out = np.zeros((len(k_inf_vacc),1))

    numrow, numcol = PhiProb.shape

    for i in range(0, numrow):
        sum = 0
        for phi_given_Itend, Itend in zip(PhiProb[i], Itend_identifier_array[i]):
            if phi_given_Itend != 0.0 and Itend != 0.0:
                prod = phi_given_Itend * Prop_Itend[Itend]
                sum = sum + prod
        Uncertainity_out[i] = sum

    Uncertainity_out = Uncertainity_out.reshape((-1,1))
    return(Uncertainity_out)


Biontech_Uncertainity = Uncertainity("C:/Users/Qba Liu/Documents/STUDIA/BIOINF_MASTER_BERLIN/SEMESTER_I/INTRODUCTION_TO_FOCUS_AREAS/COMPLEX_SYSTEMS/LAB/kinf_Biontech.txt")
np.savetxt('C:/Users/Qba Liu/Documents/STUDIA/BIOINF_MASTER_BERLIN/SEMESTER_I/INTRODUCTION_TO_FOCUS_AREAS/COMPLEX_SYSTEMS/LAB/EfficiacyUncertainityBiontech.txt',Biontech_Uncertainity,delimiter=',',fmt='%1.4f')


Moderna_Uncertainity = Uncertainity("C:/Users/Qba Liu/Documents/STUDIA/BIOINF_MASTER_BERLIN/SEMESTER_I/INTRODUCTION_TO_FOCUS_AREAS/COMPLEX_SYSTEMS/LAB/kinf_Moderna.txt")
np.savetxt('C:/Users/Qba Liu/Documents/STUDIA/BIOINF_MASTER_BERLIN/SEMESTER_I/INTRODUCTION_TO_FOCUS_AREAS/COMPLEX_SYSTEMS/LAB/EfficiacyUncertainityModerna.txt',Moderna_Uncertainity,delimiter=',',fmt='%1.4f')

Astra_Uncertainity = Uncertainity("C:/Users/Qba Liu/Documents/STUDIA/BIOINF_MASTER_BERLIN/SEMESTER_I/INTRODUCTION_TO_FOCUS_AREAS/COMPLEX_SYSTEMS/LAB/kinf_Astra.txt")
np.savetxt('C:/Users/Qba Liu/Documents/STUDIA/BIOINF_MASTER_BERLIN/SEMESTER_I/INTRODUCTION_TO_FOCUS_AREAS/COMPLEX_SYSTEMS/LAB/EfficiacyUncertainityAstra.txt',Astra_Uncertainity,delimiter=',',fmt='%1.4f')
