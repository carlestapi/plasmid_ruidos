import math
import numpy as np
import cellpylib as cpl

cellular_automaton = cpl.init_simple(25, dtype=np.float64)

def apply_rule(n, c, t):
    result = (sum(n) / len(n)) * (3 / 2)
    frac, whole = math.modf(result)
    return frac

cellular_automaton = cpl.evolve(cellular_automaton, timesteps=50,
                                apply_rule=apply_rule)

# Creating an array
List = cellular_automaton
Array = np.asarray(List, dtype='float')

# Displaying the array
print(Array)

# Saving the array in a text file
np.savetxt("contCA_50.txt", Array)

cpl.plot(cellular_automaton)
