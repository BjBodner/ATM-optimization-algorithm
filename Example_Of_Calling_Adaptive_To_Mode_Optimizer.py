
# coding: utf-8

# In[1]:


import numpy as np
import Adaptive_Two_Mode_Optimizer


# In[3]:


def fun(x): 
    return np.sum((1000000**(np.arange(D)/D))*((x-X_Optimal)**2) )


D = 2
remaining_evals = 10000*D
Number_Of_Tests = 3
test_Optimization = 1


SOLVER = Adaptive_Two_Mode_Optimizer.Adaptive_Two_Mode

print("Now running example optimization test on f2 - on a Budget of 10000*D:")
print()
for test in range (Number_Of_Tests):
    print("round ",test+1,":")
    x0 = 10*(np.random.rand(D)-0.5)
    X_Optimal = 10*(np.random.rand(D)-0.5)
    print("Initial Cost",fun(x0))

    ## This is the first method of calling the optimizer
    Best_Parameters_After_Optimization = Adaptive_Two_Mode_Optimizer.optimize(fun, x0, remaining_evals )
    ## This is the second method of calling the optimizer
    Best_Parameters_After_Optimization = SOLVER.Optimize(fun, x0, remaining_evals )
    
    print("Final Cost",fun(Best_Parameters_After_Optimization))
    print()
    
print("All Done")

