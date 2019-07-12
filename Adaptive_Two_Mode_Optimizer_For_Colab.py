
# coding: utf-8

# In[1]:



# coding: utf-8

# In[1]:



# coding: utf-8

# In[1]:


## Two Mode Section

import numpy as np
import matplotlib.pyplot as plt

np.seterr(invalid='ignore')


# In[2]:



def Two_Mode_Optimizer(Seed_Parameters,Second_Best_Seed_Parameters,NumberOfSamples,Parameter_ChangeVector,Cost_Change,Random_Search_exponential_Growth_Factor,Targeted_Search_exponential_Growth_Factor,Targeted_Search_Decay_Rate,Random_Search_Growth_Rate,Targeted_Search_Growth_Rate,Random_Search_Period,Significant_Change_Value,Maximal_Random_Search,Adaptive_Amplitude):

    Targeted_Search_Mechanism = 2

    NumberOfSamples = np.maximum(NumberOfSamples,1)## This is to avoid allocation errors
    Number_Of_Parameters_11 = len(Seed_Parameters)
    TargetedSearch = np.zeros((Number_Of_Parameters_11,NumberOfSamples) )
    Random_Search_exponential_Growth_Factor1 = Random_Search_exponential_Growth_Factor
    

    ################################# Testing if the improvement since last iteration was big enough #############
    if np.sum(np.abs(Parameter_ChangeVector)) > Significant_Change_Value:
        Random_Search_exponential_Growth_Factor = 0
        Targeted_Search_exponential_Growth_Factor = Targeted_Search_exponential_Growth_Factor + 1

    if np.sum(np.abs(Parameter_ChangeVector)) < Significant_Change_Value:
        Targeted_Search_exponential_Growth_Factor = 0
        
        
    if Print_Two_Mode_Analytics == 1:
        print("Two_Mode_Analytics:")
        print("---------------------------------------------------------------------")
        print("The np.sum(np.abs(Parameter_ChangeVector)) is ",np.sum(np.abs(Parameter_ChangeVector)), )
        print("As Compared To the Significant_Change_Value of ", Significant_Change_Value)
    ###############################################################################################################        
        
        
        
        
    ################################# This the Random Search Section #############################################
    ## Applying a Genetic Algorithm
    # This is to broadcast them to a len(InitialParameters),NumberOfSamples size matrix, for element wise multiplication
    Seed_Parameters1 = np.zeros((Number_Of_Parameters_11,NumberOfSamples))
    Seed_Parameters1[:,:] = np.array([Seed_Parameters]).T
    Second_Best_Seed_Parameters1 = np.zeros((Number_Of_Parameters_11,NumberOfSamples))    
    Second_Best_Seed_Parameters1[:,:] = np.array([Second_Best_Seed_Parameters]).T    
    
    
    Random_Search_Scalar_Amplitude = Maximal_Random_Search*np.exp(Random_Search_Growth_Rate*(np.sin(np.mod(Random_Search_exponential_Growth_Factor*np.pi/(2*Random_Search_Period),np.pi/2) + 0.01)-1))
    Random_Search_Scalar_Amplitude = np.minimum(Random_Search_Scalar_Amplitude,10)
    Random_Samples = 2*(np.random.rand(Number_Of_Parameters_11,NumberOfSamples)-0.5)
    
    ## Take half of the parameter each of the best two individuals from last iteration + add a smart 
    ## RMS addaptive amplitude search around the children of the next generation (Mutation) + at the end 
    ## subtract the initial seed parameters so that we are only left with change vectors
    
    Genes = np.random.rand(Number_Of_Parameters_11,NumberOfSamples)


    RandomSearch = (Genes > 0.5)*Seed_Parameters1 +     (Genes <= 0.5)*Second_Best_Seed_Parameters1 +     Random_Samples*Random_Search_Scalar_Amplitude*Adaptive_Amplitude - Seed_Parameters1



    
    if Print_Two_Mode_Analytics == 1:
        print("")
        print("Random_Search_Scalar_Amplitude = ",Random_Search_Scalar_Amplitude)
    ###############################################################################################################

        
    ################################# This the Targeted Search Section #############################################
    if np.max(Parameter_ChangeVector) > 0:
        RandomNumbers = np.random.rand(np.maximum(NumberOfSamples,1),1) ## To avoid allocation errors

        Random_Search_Scalar_Amplitude_1 = Maximal_Random_Search*np.exp(Random_Search_Growth_Rate*(np.sin(np.mod(Random_Search_exponential_Growth_Factor1*np.pi/(2*Random_Search_Period),np.pi/2) + 0.01)-1))
        TargetedSearch_Amplitude = 0


        if Targeted_Search_exponential_Growth_Factor*Targeted_Search_Growth_Rate -Random_Search_exponential_Growth_Factor*Targeted_Search_Decay_Rate < 100:
            TargetedSearch_Amplitude = Random_Search_Scalar_Amplitude_1*RandomNumbers*            np.exp(Targeted_Search_exponential_Growth_Factor*Targeted_Search_Growth_Rate                   -Random_Search_exponential_Growth_Factor*Targeted_Search_Decay_Rate)
        else: TargetedSearch_Amplitude = 100*RandomNumbers

        if Print_Two_Mode_Analytics == 1:
            print("Mean Targeted_Search_Scalar_Amplitude =  = ",np.mean(TargetedSearch_Amplitude))



        if np.max(np.abs(TargetedSearch_Amplitude)) > 1000: ## This is if the change vector is too big
                TargetedSearch_Amplitude = 100*RandomNumbers

        if np.isnan(np.max(TargetedSearch_Amplitude)) == 1: ## This is if the change vector is too big
                TargetedSearch_Amplitude = 100*RandomNumbers



        if Print_Two_Mode_Analytics == 1:
                print("np.mean(TargetedSearch_Amplitude) after suppression is = ",np.mean(TargetedSearch_Amplitude))


        for sample in range (1,NumberOfSamples):

            TargetedSearch[:,sample] = TargetedSearch_Amplitude[sample]*Parameter_ChangeVector

    ###############################################################################################################
    
    
    ################################ This is the Final Readout and generating Suggestions ########################
    if Print_Two_Mode_Analytics == 1:   
        print("")
        print(" The Adaptive_Amplitude Vector is :",Adaptive_Amplitude[:,0])
        print("---------------------------------------------------------------------")
        print("")
    

    Suggested_Parameter_Samples_From_Two_Mode = np.array([Seed_Parameters]).T + TargetedSearch + RandomSearch
    Suggested_Parameter_Samples_From_Two_Mode[:,0] = Seed_Parameters
    
    Random_Search_exponential_Growth_Factor = Random_Search_exponential_Growth_Factor +1
    ##############################################################################################################

    return Suggested_Parameter_Samples_From_Two_Mode, Random_Search_exponential_Growth_Factor, Targeted_Search_exponential_Growth_Factor







def Two_Mode_Optimizer2(Seed_Parameters,Best_Quarter_Of_Individuals,NumberOfSamples,Parameter_ChangeVector,Cost_Change,Random_Search_exponential_Growth_Factor,Targeted_Search_exponential_Growth_Factor,Targeted_Search_Decay_Rate,Random_Search_Growth_Rate,Targeted_Search_Growth_Rate,Random_Search_Period,Significant_Change_Value,Maximal_Random_Search,Adaptive_Amplitude):


    Targeted_Search_Mechanism = 2

    NumberOfSamples = np.maximum(NumberOfSamples,1)## This is to avoid allocation errors
    Number_Of_Parameters_11 = len(Seed_Parameters)
    TargetedSearch = np.zeros((Number_Of_Parameters_11,NumberOfSamples))
    Random_Search_exponential_Growth_Factor1 = Random_Search_exponential_Growth_Factor
    

    ################################# Testing if the improvement since last iteration was big enough #############
    if np.sum(np.abs(Parameter_ChangeVector)) > Significant_Change_Value:
        Random_Search_exponential_Growth_Factor = 0
        Targeted_Search_exponential_Growth_Factor = Targeted_Search_exponential_Growth_Factor + 1

    if np.sum(np.abs(Parameter_ChangeVector)) < Significant_Change_Value:
        Targeted_Search_exponential_Growth_Factor = 0
        
        
    if Print_Two_Mode_Analytics == 1:
        print("Two_Mode_Analytics:")
        print("---------------------------------------------------------------------")
        print("The np.sum(np.abs(Parameter_ChangeVector)) is ",np.sum(np.abs(Parameter_ChangeVector)), )
        print("As Compared To the Significant_Change_Value of ", Significant_Change_Value)
    ###############################################################################################################        
        
        
        
        
    ################################# This the Random Search Section #############################################
    ## Applying a Genetic Algorithm

    ## Crossover
    RandomSearch = np.zeros((Number_Of_Parameters_11,NumberOfSamples))
    if Best_Quarter_Of_Individuals.shape[1] > 0:
        Seed_Parameters1 = np.zeros((Number_Of_Parameters_11,NumberOfSamples))
        Seed_Parameters1[:,:] = np.array([Seed_Parameters]).T
        Crossover_Indexes = np.random.randint(0,Best_Quarter_Of_Individuals.shape[1],(Number_Of_Parameters_11,NumberOfSamples))
        for i in range (0,Number_Of_Parameters_11):
            RandomSearch[i,:] = Best_Quarter_Of_Individuals[i,Crossover_Indexes[i,:]]

        ## Mutation
        Random_Search_Scalar_Amplitude = Maximal_Random_Search*np.exp(Random_Search_Growth_Rate*(np.sin(np.mod(Random_Search_exponential_Growth_Factor*np.pi/(2*Random_Search_Period),np.pi/2) + 0.01)-1))    
        Random_Search_Scalar_Amplitude = np.minimum(Random_Search_Scalar_Amplitude,10)
        RandomSearch += 2*(np.random.rand(Number_Of_Parameters_11,NumberOfSamples)-0.5)*Random_Search_Scalar_Amplitude*Adaptive_Amplitude - Seed_Parameters1

    
    
    if Print_Two_Mode_Analytics == 1:
        print("")
        print("Random_Search_Scalar_Amplitude = ",Random_Search_Scalar_Amplitude)
    ###############################################################################################################

        
    ################################# This the Targeted Search Section #############################################
    if np.max(Parameter_ChangeVector) > 0:
        RandomNumbers = np.random.rand(np.maximum(NumberOfSamples,1),1) ## To avoid allocation errors

        Random_Search_Scalar_Amplitude_1 = Maximal_Random_Search*np.exp(Random_Search_Growth_Rate*(np.sin(np.mod(Random_Search_exponential_Growth_Factor1*np.pi/(2*Random_Search_Period),np.pi/2) + 0.01)-1))
        TargetedSearch_Amplitude = 0


        if Targeted_Search_exponential_Growth_Factor*Targeted_Search_Growth_Rate -Random_Search_exponential_Growth_Factor*Targeted_Search_Decay_Rate < 100:
            TargetedSearch_Amplitude = Random_Search_Scalar_Amplitude_1*RandomNumbers*            np.exp(Targeted_Search_exponential_Growth_Factor*Targeted_Search_Growth_Rate                   -Random_Search_exponential_Growth_Factor*Targeted_Search_Decay_Rate)
        else: TargetedSearch_Amplitude = 100*RandomNumbers


        if Print_Two_Mode_Analytics == 1:
            print("Mean Targeted_Search_Scalar_Amplitude =  = ",np.mean(TargetedSearch_Amplitude))



        if np.max(np.abs(TargetedSearch_Amplitude)) > 1000: ## This is if the change vector is too big
                TargetedSearch_Amplitude = 100*RandomNumbers

        if np.isnan(np.max(TargetedSearch_Amplitude)) == 1: ## This is if the change vector is too big
                TargetedSearch_Amplitude = 100*RandomNumbers



        if Print_Two_Mode_Analytics == 1:
                print("np.mean(TargetedSearch_Amplitude) after suppression is = ",np.mean(TargetedSearch_Amplitude))


        for sample in range (1,NumberOfSamples):

            TargetedSearch[:,sample] = TargetedSearch_Amplitude[sample]*Parameter_ChangeVector

    ###############################################################################################################
    
    
    ################################ This is the Final Readout and generating Suggestions ########################
    if Print_Two_Mode_Analytics == 1:   
        print("")
        print(" The Adaptive_Amplitude Vector is :",Adaptive_Amplitude[:,0])
        print("---------------------------------------------------------------------")
        print("")
    

    Suggested_Parameter_Samples_From_Two_Mode = np.array([Seed_Parameters]).T + TargetedSearch + RandomSearch
    Suggested_Parameter_Samples_From_Two_Mode[:,0] = Seed_Parameters
    
    Random_Search_exponential_Growth_Factor = Random_Search_exponential_Growth_Factor +1
    ##############################################################################################################

    return Suggested_Parameter_Samples_From_Two_Mode, Random_Search_exponential_Growth_Factor, Targeted_Search_exponential_Growth_Factor


# In[3]:


## Conductor Algorithm

def Symmetric_Section_Conductor(Current_Resource_Allocation,Changevector,Hyper_Parameters):
    ## The difference is that this function can allow for different masses and interactions between the algorithms
    
    
    Initial_Resource_Allocation = Hyper_Parameters[:,0]
    MassVector = Hyper_Parameters[:,1]
    Self_Interaction_Spring_Constants = Hyper_Parameters[:,2]
    InteractionMatrix = Hyper_Parameters[:,3:3+Number_Of_Algorithms]

    
    
    Norm = np.sqrt(np.matmul(Changevector,Changevector))
    if Norm != 0:
        Changevector = Changevector/Norm
    
    R0 = Initial_Resource_Allocation
    R = Current_Resource_Allocation
    K0 = Self_Interaction_Spring_Constants

    ## This updates the allocated resources
    H_Self = -(K0/MassVector)*(R-R0) ## This slowly restores the Resouce allocation to the original values
    H_Interaction = np.matmul(InteractionMatrix,(Changevector/MassVector).T) ## this is a resouce conserving interaction
    R = R + H_Self + H_Interaction.T
        

    RecomendedResourceAllocation1 = np.round(R)

    ## This makes sure that if the total resources aren't conserved then add or subtract from the best Member.
    I = np.argmax(RecomendedResourceAllocation1)
    RecomendedResourceAllocation1[I] = RecomendedResourceAllocation1[I] + (np.sum(Initial_Resource_Allocation) - np.sum(RecomendedResourceAllocation1))

    ## This convertes the array to an int array
    RecomendedResourceAllocation = np.arange(len(RecomendedResourceAllocation1))
    for i in range(len(RecomendedResourceAllocation1)):
        if np.isnan(RecomendedResourceAllocation1[i]):
            Changevector = np.zeros(4)
            RecomendedResourceAllocation1 = 16*np.ones(4)
        RecomendedResourceAllocation[i] = int(RecomendedResourceAllocation1[i])
    
    return RecomendedResourceAllocation



def Set_Hyperparameters_For_Symmetric_Conductor(Total_Resources,Number_Of_Algorithms,Mass,Self_Spring_Constants,Interaction_Spring_Constants):
    

    Initial_Resource_Allocation =int(Total_Resources/Number_Of_Algorithms)*np.ones(Number_Of_Algorithms)

    ## This fixes Rounding issues
    Initial_Resource_Allocation[0] = Initial_Resource_Allocation[0] + (Total_Resources-np.sum(Initial_Resource_Allocation))
    MassVector = Mass*np.ones(Number_Of_Algorithms)
    Self_Interaction_Spring_Constants = Self_Spring_Constants*np.ones(Number_Of_Algorithms)
    InteractionMatrix = Interaction_Spring_Constants*(-np.ones((Number_Of_Algorithms,Number_Of_Algorithms)) + Number_Of_Algorithms*np.diag(np.ones(Number_Of_Algorithms)))

    
    
    Hyperparameters_Array = np.zeros((Number_Of_Algorithms,3+Number_Of_Algorithms))
    Hyperparameters_Array[:,0] = Initial_Resource_Allocation
    Hyperparameters_Array[:,1] = MassVector
    Hyperparameters_Array[:,2] = Self_Interaction_Spring_Constants
    Hyperparameters_Array[:,3:3+Number_Of_Algorithms] = InteractionMatrix
    
    
    


    return Hyperparameters_Array



### Test of the conductor algorithm

Total_Resources = 60
Number_Of_Algorithms = 4
Mass =10
Self_Spring_Constants = 10
Interaction_Spring_Constants = 40
    
    
Hyperparameters_For_Symmetric_Conductor = Set_Hyperparameters_For_Symmetric_Conductor(Total_Resources,Number_Of_Algorithms,Mass,Self_Spring_Constants,Interaction_Spring_Constants)
Initial_Resource_Allocation = Hyperparameters_For_Symmetric_Conductor[:,0]
Current_Resource_Allocation = Hyperparameters_For_Symmetric_Conductor[:,0]

test_resource_allocation = 0
if test_resource_allocation == 1:
    print("Initial Total_Resources",np.sum(Initial_Resource_Allocation))
    for i in range (1,50):
        Changevector = np.array([-1, 0 ,0,0]) ## Negative is bad
        Recomended_Resource_Allocation = Symmetric_Section_Conductor(Current_Resource_Allocation,Changevector,Hyperparameters_For_Symmetric_Conductor) 
        Current_Resource_Allocation = Recomended_Resource_Allocation


    Minimal_Allocation = np.min(Current_Resource_Allocation)
    print("")
    print("Minimum Resouces test")
    print("------------------------------------------------------------------")
    print("Minimal_Allocation is:",Minimal_Allocation)
    print("(Maximal_Allocation is:",np.max(Current_Resource_Allocation),")")
    print("Total_Resources After minumum resource allocation Test",np.sum(Current_Resource_Allocation))

    for i in range (1,50):
        Changevector = np.array([1,0 ,0,0]) ## Positive is good
        Recomended_Resource_Allocation = Symmetric_Section_Conductor(Current_Resource_Allocation,Changevector,Hyperparameters_For_Symmetric_Conductor)
        Current_Resource_Allocation = Recomended_Resource_Allocation

    Maximal_Allocation = np.max(Current_Resource_Allocation)
    print("")
    print("Maximum Resouces test")
    print("------------------------------------------------------------------")
    print("Maximal_Allocation is:",Maximal_Allocation)
    print("(Minimal_Allocation is:",np.min(Current_Resource_Allocation),")")
    print("Total_Resources After maximum resource allocation Test",np.sum(Current_Resource_Allocation))


# In[4]:


class Two_Mode_Suggestion_Tools(object):

    # The class "constructor" - It's actually an initializer 
    def __init__(self, Best_Individual, Second_Best_Individual, Parameter_ChangeVector, Number_Of_Samples,Cost_Change,Random_Search_exponential_Growth_Factor,Targeted_Search_exponential_Growth_Factor,Adaptive_Amplitude_Vector,Squared_Gradient):
        self.Best_Individual = Best_Individual
        self.Second_Best_Individual = Second_Best_Individual
        self.Parameter_ChangeVector = Parameter_ChangeVector
        self.Number_Of_Samples = Number_Of_Samples
        self.Cost_Change = Cost_Change
        self.Random_Search_exponential_Growth_Factor = Random_Search_exponential_Growth_Factor
        self.Targeted_Search_exponential_Growth_Factor = Targeted_Search_exponential_Growth_Factor
        self.Adaptive_Amplitude_Vector = Adaptive_Amplitude_Vector
        self.Squared_Gradient = Squared_Gradient

        
def Create_New_Two_Mode_Suggestion_Tools_Object(Number_Of_Parameters): ## this is where the directional pool is initialized
    N = Number_Of_Parameters
    Best_Individual = np.random.rand(N) - 0.5
    Second_Best_Individual = np.random.rand(N) - 0.5
    Parameter_ChangeVector = np.random.rand(N) - 0.5
    Number_Of_Samples = 10
    Cost_Change = -1
    Random_Search_exponential_Growth_Factor = 1.0
    Targeted_Search_exponential_Growth_Factor = 1.0
    Adaptive_Amplitude_Vector = np.random.rand(N) - 0.5
    Squared_Gradient = np.random.rand(N)
    ## Generating and instance For The object the Two mode Suggestions tools object
    return  Two_Mode_Suggestion_Tools(Best_Individual, Second_Best_Individual, Parameter_ChangeVector, Number_Of_Samples,Cost_Change,Random_Search_exponential_Growth_Factor,Targeted_Search_exponential_Growth_Factor,Adaptive_Amplitude_Vector,Squared_Gradient)


Test_Object_Functionality = 0
if Test_Object_Functionality == 1:
    ## Setting the inputs for the Two mode Suggestions tools object
    N = 30
    Best_Individual = np.random.rand(N) - 0.5
    Second_Best_Individual = np.random.rand(N) - 0.5
    Parameter_ChangeVector = np.random.rand(N) - 0.5
    Number_Of_Samples = 10
    Cost_Change = -1
    Random_Search_exponential_Growth_Factor = 2.0
    Targeted_Search_exponential_Growth_Factor = 8.0
    Adaptive_Amplitude_Vector = np.random.rand(N) - 0.5
    Squared_Gradient = np.zeros(N)

    ## Generating and instance For The object the Two mode Suggestions tools object
    Suggestion_Tools =  Two_Mode_Suggestion_Tools(Best_Individual, Second_Best_Individual, Parameter_ChangeVector, Number_Of_Samples,Cost_Change,Random_Search_exponential_Growth_Factor,Targeted_Search_exponential_Growth_Factor,Adaptive_Amplitude_Vector,Squared_Gradient)

    Best_Individual2 = Suggestion_Tools.Best_Individual
    Second_Best_Individual2 = Suggestion_Tools.Second_Best_Individual
    Parameter_ChangeVector2 = Suggestion_Tools.Parameter_ChangeVector
    Number_Of_Samples2 = Suggestion_Tools.Number_Of_Samples
    Cost_Change2 = Suggestion_Tools.Cost_Change
    Random_Search_exponential_Growth_Factor2 = Suggestion_Tools.Random_Search_exponential_Growth_Factor
    Targeted_Search_exponential_Growth_Factor2 = Suggestion_Tools.Targeted_Search_exponential_Growth_Factor
    Adaptive_Amplitude_Vector2 = Suggestion_Tools.Adaptive_Amplitude_Vector
    Adaptive_Amplitude2 = np.ones((len(Best_Individual),Number_Of_Samples))
    Adaptive_Amplitude2[:,:] = np.array([Adaptive_Amplitude_Vector2]).T
    Squared_Gradient = np.zeros(len(Best_Individual))




    print("Test If the Items were loaded Properly into the object:")
    print("------------------------------------------------------------")
    print("")
    print(Best_Individual == Best_Individual2)
    print("")
    print(Second_Best_Individual == Second_Best_Individual2)
    print("")
    print(Parameter_ChangeVector == Parameter_ChangeVector2)
    print("")
    print(Number_Of_Samples == Number_Of_Samples2)
    print("")
    print(Cost_Change == Cost_Change2)
    print("")
    print(Random_Search_exponential_Growth_Factor == Random_Search_exponential_Growth_Factor2)
    print("")
    print(Targeted_Search_exponential_Growth_Factor == Targeted_Search_exponential_Growth_Factor2)
    print("")
    print(Adaptive_Amplitude_Vector == Adaptive_Amplitude_Vector2)
    print("------------------------------------------------------------")

    print("")
    print("")
    print("Test if we can update the values in the object propely")
    print("------------------------------------------------------------")
    print("Suggestion_Tools.Number_Of_Samples = ",Suggestion_Tools.Number_Of_Samples)
    print("Suggestion_Tools.Number_Of_Samples = ",2)
    Suggestion_Tools.Number_Of_Samples = 2
    print("Suggestion_Tools.Number_Of_Samples = ",Suggestion_Tools.Number_Of_Samples)


    print("Suggestion_Tools.Best_Individual = ",Suggestion_Tools.Best_Individual)
    print("Suggestion_Tools.Best_Individual = 10* Suggestion_Tools.Best_Individual")
    Suggestion_Tools.Best_Individual = 10*Suggestion_Tools.Best_Individual
    print("Suggestion_Tools.Best_Individual = ",Suggestion_Tools.Best_Individual)



    print("------------------------------------------------------------")


# In[5]:



def Get_TWO_Mode_Suggestions(fun,Hyperparameters,Suggestion_Tools,Print_Two_Mode_Analytics):
   

   

   Random_Search_Growth_Rate  = Hyperparameters[0]
   Random_Search_Period = Hyperparameters[1]
   Maximal_Random_Search = Hyperparameters[2]
   Targeted_Search_Growth_Rate = Hyperparameters[3]
   Targeted_Search_Decay_Rate = Hyperparameters[4]
   Significant_Change_Value = Hyperparameters[5]
   alpha = Hyperparameters[6]
   beta = Hyperparameters[7]

   
   ######################################## There can also be speed up here by just inputing it directly into the function
   Best_Individual = Suggestion_Tools.Best_Individual
   Second_Best_Individual = Suggestion_Tools.Second_Best_Individual
   Parameter_ChangeVector = Suggestion_Tools.Parameter_ChangeVector
   Number_Of_Samples = Suggestion_Tools.Number_Of_Samples
   Cost_Change = Suggestion_Tools.Cost_Change
   Random_Search_exponential_Growth_Factor = Suggestion_Tools.Random_Search_exponential_Growth_Factor
   Targeted_Search_exponential_Growth_Factor = Suggestion_Tools.Targeted_Search_exponential_Growth_Factor
   Adaptive_Amplitude1 = Suggestion_Tools.Adaptive_Amplitude_Vector
   Adaptive_Amplitude = np.ones((len(Best_Individual),Number_Of_Samples))
   Adaptive_Amplitude[:,:] = np.array([Adaptive_Amplitude1]).T

   
   Suggested_Parameter_Samples_From_Two_Mode, Random_Search_exponential_Growth_Factor, Targeted_Search_exponential_Growth_Factor =     Two_Mode_Optimizer(Best_Individual,Second_Best_Individual,Number_Of_Samples,Parameter_ChangeVector,Cost_Change, Random_Search_exponential_Growth_Factor,Targeted_Search_exponential_Growth_Factor,Targeted_Search_Decay_Rate,Random_Search_Growth_Rate,Targeted_Search_Growth_Rate,Random_Search_Period,Significant_Change_Value,Maximal_Random_Search,Adaptive_Amplitude)

   Suggestion_Tools.Random_Search_exponential_Growth_Factor = Random_Search_exponential_Growth_Factor
   Suggestion_Tools.Targeted_Search_exponential_Growth_Factor = Targeted_Search_exponential_Growth_Factor
   
   return Suggested_Parameter_Samples_From_Two_Mode,Suggestion_Tools





def Get_TWO_Mode_Suggestions2(fun,Hyperparameters,Suggestion_Tools,Print_Two_Mode_Analytics,Best_Quarter_Of_Individuals):
   

   

   Random_Search_Growth_Rate  = Hyperparameters[0]
   Random_Search_Period = Hyperparameters[1]
   Maximal_Random_Search = Hyperparameters[2]
   Targeted_Search_Growth_Rate = Hyperparameters[3]
   Targeted_Search_Decay_Rate = Hyperparameters[4]
   Significant_Change_Value = Hyperparameters[5]
   alpha = Hyperparameters[6]
   beta = Hyperparameters[7]

   
   ######################################## There can also be speed up here by just inputing it directly into the function
   Best_Individual = Suggestion_Tools.Best_Individual
   Second_Best_Individual = Suggestion_Tools.Second_Best_Individual
   Parameter_ChangeVector = Suggestion_Tools.Parameter_ChangeVector
   Number_Of_Samples = Suggestion_Tools.Number_Of_Samples
   Cost_Change = Suggestion_Tools.Cost_Change
   Random_Search_exponential_Growth_Factor = Suggestion_Tools.Random_Search_exponential_Growth_Factor
   Targeted_Search_exponential_Growth_Factor = Suggestion_Tools.Targeted_Search_exponential_Growth_Factor
   Adaptive_Amplitude1 = Suggestion_Tools.Adaptive_Amplitude_Vector
   Adaptive_Amplitude = np.ones((len(Best_Individual),Number_Of_Samples))
   Adaptive_Amplitude[:,:] = np.array([Adaptive_Amplitude1]).T

   
   Suggested_Parameter_Samples_From_Two_Mode, Random_Search_exponential_Growth_Factor, Targeted_Search_exponential_Growth_Factor =     Two_Mode_Optimizer2(Best_Individual,Best_Quarter_Of_Individuals,Number_Of_Samples,Parameter_ChangeVector,Cost_Change, Random_Search_exponential_Growth_Factor,Targeted_Search_exponential_Growth_Factor,Targeted_Search_Decay_Rate,Random_Search_Growth_Rate,Targeted_Search_Growth_Rate,Random_Search_Period,Significant_Change_Value,Maximal_Random_Search,Adaptive_Amplitude)

   Suggestion_Tools.Random_Search_exponential_Growth_Factor = Random_Search_exponential_Growth_Factor
   Suggestion_Tools.Targeted_Search_exponential_Growth_Factor = Targeted_Search_exponential_Growth_Factor
   
   return Suggested_Parameter_Samples_From_Two_Mode,Suggestion_Tools


# In[6]:



def fun(x):
    return np.sum((x-2)**2)


def Get_Suggestions_From_Two_Mode_Section(fun,Hyperparameters_1,Hyperparameters_2,Hyperparameters_3,Hyperparameters_4,Suggestion_Tools_1,Suggestion_Tools_2,Suggestion_Tools_3,Suggestion_Tools_4,Hyperparameters_For_Symmetric_Conductor):

    ## This gets The Suggestions from The Different Sets Of hyperparameters Set 1

    Suggested_Parameter_Samples_From_Two_Mode_1,Suggestion_Tools_1  = Get_TWO_Mode_Suggestions(fun,Hyperparameters_1,Suggestion_Tools_1 ,Print_Two_Mode_Analytics)   
    Suggested_Parameter_Samples_From_Two_Mode_2,Suggestion_Tools_2  = Get_TWO_Mode_Suggestions(fun,Hyperparameters_2,Suggestion_Tools_2 ,Print_Two_Mode_Analytics)   
    Suggested_Parameter_Samples_From_Two_Mode_3,Suggestion_Tools_3  = Get_TWO_Mode_Suggestions(fun,Hyperparameters_3,Suggestion_Tools_3 ,Print_Two_Mode_Analytics)   
    Suggested_Parameter_Samples_From_Two_Mode_4,Suggestion_Tools_4  = Get_TWO_Mode_Suggestions(fun,Hyperparameters_4,Suggestion_Tools_4 ,Print_Two_Mode_Analytics) 


    ## This is to orginize the different indexes of the suggestions for later processing
    N1 = len(Suggested_Parameter_Samples_From_Two_Mode_1[0,:])
    N2 = len(Suggested_Parameter_Samples_From_Two_Mode_2[0,:])
    N3 = len(Suggested_Parameter_Samples_From_Two_Mode_3[0,:])
    N4 = len(Suggested_Parameter_Samples_From_Two_Mode_4[0,:])

    ## Save the indexes in the Suggestions array to know which samples came from which indvidual
    Indexes_Of_The_Suggestions_From_The_Different_Algorithms = np.array([0,N1,N1+N2,N1+N2+N3,N1+N2+N3+N4])
    ## Combine all the suggested samples into one array
    All_Suggestions_From_Two_Mode_Section = np.concatenate((Suggested_Parameter_Samples_From_Two_Mode_1,Suggested_Parameter_Samples_From_Two_Mode_2,Suggested_Parameter_Samples_From_Two_Mode_3,Suggested_Parameter_Samples_From_Two_Mode_4),1)


    return All_Suggestions_From_Two_Mode_Section,Indexes_Of_The_Suggestions_From_The_Different_Algorithms








def Get_Suggestions_From_Two_Mode_Section2(fun,Hyperparameters_1,Hyperparameters_2,Hyperparameters_3,Hyperparameters_4,Suggestion_Tools_1,Suggestion_Tools_2,Suggestion_Tools_3,Suggestion_Tools_4,Hyperparameters_For_Symmetric_Conductor,Best_Quarter_Of_Individuals):

    ## This gets The Suggestions from The Different Sets Of hyperparameters Set 1

    Suggested_Parameter_Samples_From_Two_Mode_1,Suggestion_Tools_1  = Get_TWO_Mode_Suggestions2(fun,Hyperparameters_1,Suggestion_Tools_1 ,Print_Two_Mode_Analytics,Best_Quarter_Of_Individuals)   
    Suggested_Parameter_Samples_From_Two_Mode_2,Suggestion_Tools_2  = Get_TWO_Mode_Suggestions2(fun,Hyperparameters_2,Suggestion_Tools_2 ,Print_Two_Mode_Analytics,Best_Quarter_Of_Individuals)   
    Suggested_Parameter_Samples_From_Two_Mode_3,Suggestion_Tools_3  = Get_TWO_Mode_Suggestions2(fun,Hyperparameters_3,Suggestion_Tools_3 ,Print_Two_Mode_Analytics,Best_Quarter_Of_Individuals)   
    Suggested_Parameter_Samples_From_Two_Mode_4,Suggestion_Tools_4  = Get_TWO_Mode_Suggestions2(fun,Hyperparameters_4,Suggestion_Tools_4 ,Print_Two_Mode_Analytics,Best_Quarter_Of_Individuals) 


    ## This is to orginize the different indexes of the suggestions for later processing
    N1 = len(Suggested_Parameter_Samples_From_Two_Mode_1[0,:])
    N2 = len(Suggested_Parameter_Samples_From_Two_Mode_2[0,:])
    N3 = len(Suggested_Parameter_Samples_From_Two_Mode_3[0,:])
    N4 = len(Suggested_Parameter_Samples_From_Two_Mode_4[0,:])

    ## Save the indexes in the Suggestions array to know which samples came from which indvidual
    Indexes_Of_The_Suggestions_From_The_Different_Algorithms = np.array([0,N1,N1+N2,N1+N2+N3,N1+N2+N3+N4])
    ## Combine all the suggested samples into one array
    All_Suggestions_From_Two_Mode_Section = np.concatenate((Suggested_Parameter_Samples_From_Two_Mode_1,Suggested_Parameter_Samples_From_Two_Mode_2,Suggested_Parameter_Samples_From_Two_Mode_3,Suggested_Parameter_Samples_From_Two_Mode_4),1)


    return All_Suggestions_From_Two_Mode_Section,Indexes_Of_The_Suggestions_From_The_Different_Algorithms


# In[7]:



     
def Update_Addaptive_Amplitudes(Gradients,Hyperparameters_1,Hyperparameters_2,Hyperparameters_3,Hyperparameters_4,Suggestion_Tools_1,Suggestion_Tools_2,Suggestion_Tools_3,Suggestion_Tools_4):
         

## This is to do these two lines for all individuals
#   alpha = Hyperparameters_1[6]
#   beta = Hyperparameters_1[7]
#   S = beta*S + (1-beta)*(np.mean(np.abs(Gradients),1))**2
#   Adaptive_Amplitude1 = alpha/np.sqrt(S + (alpha**2))
#   a = Hyperparameters_1[7]*Suggestion_Tools_1.Squared_Gradient + (1-Hyperparameters_1[7])*(np.mean(np.abs(Gradients,1))**2)
 if np.max(np.abs(Gradients)) < 10**20:
     RMS = np.mean(Gradients - (np.array([np.mean(Gradients ,1)]).T)**2,1) ## actually its just the square mean
 else: RMS = np.zeros(Gradients.shape[0])
 ## Update The First Individaul
 Suggestion_Tools_1.Squared_Gradient = Hyperparameters_1[7]*Suggestion_Tools_1.Squared_Gradient + (1-Hyperparameters_1[7])*RMS
 Suggestion_Tools_1.Squared_Gradient = Suggestion_Tools_1.Squared_Gradient*np.sign(np.mean(Gradients,1))
 Suggestion_Tools_1.Adaptive_Amplitude = Hyperparameters_1[6]/np.sqrt(np.maximum(Suggestion_Tools_1.Squared_Gradient + (Hyperparameters_1[6]**2),0.1*(Hyperparameters_1[6]**2)))


 ## Update The Second Individaul
 Suggestion_Tools_2.Squared_Gradient = Hyperparameters_2[7]*Suggestion_Tools_2.Squared_Gradient + (1-Hyperparameters_2[7])*RMS
 Suggestion_Tools_2.Squared_Gradient = Suggestion_Tools_2.Squared_Gradient*np.sign(np.mean(Gradients,1))
 Suggestion_Tools_2.Adaptive_Amplitude = Hyperparameters_2[6]/np.sqrt(np.maximum(Suggestion_Tools_2.Squared_Gradient + (Hyperparameters_2[6]**2),0.1*(Hyperparameters_2[6]**2)))

 ## Update The Third  Individaul
 Suggestion_Tools_3.Squared_Gradient = Hyperparameters_3[7]*Suggestion_Tools_3.Squared_Gradient + (1-Hyperparameters_3[7])*RMS
 Suggestion_Tools_3.Squared_Gradient = Suggestion_Tools_3.Squared_Gradient*np.sign(np.mean(Gradients,1))
 Suggestion_Tools_3.Adaptive_Amplitude = Hyperparameters_3[6]/np.sqrt(np.maximum(Suggestion_Tools_3.Squared_Gradient + (Hyperparameters_3[6]**2),0.1*(Hyperparameters_3[6]**2)))

 ## Update The Third  Individaul
 Suggestion_Tools_4.Squared_Gradient = Hyperparameters_4[7]*Suggestion_Tools_4.Squared_Gradient + (1-Hyperparameters_4[7])*RMS
 Suggestion_Tools_4.Squared_Gradient = Suggestion_Tools_4.Squared_Gradient*np.sign(np.mean(Gradients,1))
 Suggestion_Tools_4.Adaptive_Amplitude = Hyperparameters_4[6]/np.sqrt(np.maximum(Suggestion_Tools_4.Squared_Gradient + (Hyperparameters_4[6]**2),0.1*(Hyperparameters_4[6]**2)))


                                                                                                                                
                                                                                                                                

def Get_Best_Cost_Change_From_Individuals(Cost_Change_From_Probe_Samples,Indexes_Of_The_Suggestions_From_The_Different_Algorithms1):

 Best_Cost_Change_From_Individuals  = np.zeros(4)
 Indexes = Indexes_Of_The_Suggestions_From_The_Different_Algorithms1

 
 ### There is potential for speed up by writing this out explicitely
 for i in range (0,4):                                                                                                                       
     index_min = np.argmin(Cost_Change_From_Probe_Samples[Indexes[i]:Indexes[i+1]])
     Best_Cost_Change_From_Individuals[i] = Cost_Change_From_Probe_Samples[Indexes[i]:Indexes[i+1]][index_min]

 return  Best_Cost_Change_From_Individuals                                                                                

                                                                                                                                
                                                                                                                                
                                                                                                                     
      ## Update the Resouce Allocation
def  Update_Resource_Allocation(Cost_Change_From_Probe_Samples,Indexes_Of_The_Suggestions_From_The_Different_Algorithms1,Suggestion_Tools_1,Suggestion_Tools_2,Suggestion_Tools_3,Suggestion_Tools_4,Hyperparameters_For_Symmetric_Conductor):
                               
 # Normalize Cost
 if np.sum(np.abs(Cost_Change_From_Probe_Samples)) != 0:
     Cost_Change_From_Probe_Samples = Cost_Change_From_Probe_Samples/(np.sqrt(np.matmul(Cost_Change_From_Probe_Samples,Cost_Change_From_Probe_Samples)))
                                                                                                                           
 Best_Cost_Change_From_Individuals  = Get_Best_Cost_Change_From_Individuals(Cost_Change_From_Probe_Samples,Indexes_Of_The_Suggestions_From_The_Different_Algorithms1)                                                                                                                     
 Current_Resource_Allocation = np.array([Suggestion_Tools_1.Number_Of_Samples                                          ,Suggestion_Tools_2.Number_Of_Samples                                           ,Suggestion_Tools_3.Number_Of_Samples                                           ,Suggestion_Tools_4.Number_Of_Samples ])

 if print_Cost == 1:
     print("Current_Resource_Allocation",Current_Resource_Allocation)

 Recomended_Resource_Allocation = Symmetric_Section_Conductor(Current_Resource_Allocation,-Best_Cost_Change_From_Individuals,Hyperparameters_For_Symmetric_Conductor) 
 Current_Resource_Allocation = Recomended_Resource_Allocation

 Suggestion_Tools_1.Number_Of_Samples =      Current_Resource_Allocation[0]                                                                                                                     
 Suggestion_Tools_2.Number_Of_Samples =      Current_Resource_Allocation[1] 
 Suggestion_Tools_3.Number_Of_Samples =      Current_Resource_Allocation[2]                                                                                                                                  
 Suggestion_Tools_4.Number_Of_Samples =      Current_Resource_Allocation[3]    

 return Best_Cost_Change_From_Individuals



def Update_Best_Individuals_And_Costs_From_This_Round(Best_Individual,Second_Best_Individual,Parameter_ChangeVector,Cost_Change,Suggestion_Tools_1,Suggestion_Tools_2,Suggestion_Tools_3,Suggestion_Tools_4):                                                                                                                           


 Suggestion_Tools_1.Best_Individual = Best_Individual                                                                                                                    
 Suggestion_Tools_1.Second_Best_Individual = Second_Best_Individual
 Suggestion_Tools_1.Parameter_ChangeVector = Parameter_ChangeVector[:,0]                                                                                                                    
 Suggestion_Tools_1.Cost_Change = Cost_Change                                                                                                                     

 Suggestion_Tools_2.Best_Individual = Best_Individual                                                                                                                    
 Suggestion_Tools_2.Second_Best_Individual = Second_Best_Individual
 Suggestion_Tools_2.Parameter_ChangeVector = Parameter_ChangeVector [:,0]                                                                                                                   
 Suggestion_Tools_2.Cost_Change = Cost_Change                                                                                                                                        

 Suggestion_Tools_3.Best_Individual = Best_Individual                                                                                                                    
 Suggestion_Tools_3.Second_Best_Individual = Second_Best_Individual
 Suggestion_Tools_3.Parameter_ChangeVector = Parameter_ChangeVector[:,0]                                                                                                                    
 Suggestion_Tools_3.Cost_Change = Cost_Change                                                                                                                                        

 Suggestion_Tools_4.Best_Individual = Best_Individual                                                                                                                    
 Suggestion_Tools_4.Second_Best_Individual = Second_Best_Individual
 Suggestion_Tools_4.Parameter_ChangeVector = Parameter_ChangeVector[:,0]                                                                                                                    
 Suggestion_Tools_4.Cost_Change = Cost_Change  


 
def Optimize_With_TWO_Mode_Section(Best_Individual,Second_Best_Individual,PreviousCost,fun,current_iteration,Total_Number_Of_Iterations,Number_Of_Iterations_With_These_Hyperparameters,Hyperparameters_1,Hyperparameters_2,Hyperparameters_3,Hyperparameters_4,Suggestion_Tools_1,Suggestion_Tools_2,Suggestion_Tools_3,Suggestion_Tools_4,Number_Of_Evaluations,Allowed_Number_Of_Function_Evaluations,ForceStop,Cost_Tracker ,Best_Quarter_Of_Individuals):

 
 Restart = 0
 Best_Cost_Change_From_Individuals1 = np.zeros((Number_Of_Algorithms,Number_Of_Iterations_With_These_Hyperparameters))
      
 if current_iteration == 0:
     Parameter_ChangeVector = np.zeros((len(Best_Individual),1))
     Cost_Change = 0
     Update_Best_Individuals_And_Costs_From_This_Round(Best_Individual,Second_Best_Individual,Parameter_ChangeVector,Cost_Change,Suggestion_Tools_1,Suggestion_Tools_2,Suggestion_Tools_3,Suggestion_Tools_4)                                                                                                                           

                                                                                                                       
 for i in range(0,Number_Of_Iterations_With_These_Hyperparameters):

         
         
     current_iteration= current_iteration + 1
     ## Readout ony if wanted                                                                                                                           
     if print_Cost == 1:
         print("____________________________________________________")
         print("Iteration",current_iteration,"\\",Total_Number_Of_Iterations)

     
                                                                                                                             
     ##  Get Sample Suggestions from the Section    
     if len(Cost_Tracker) <=1:
         All_Suggestions_From_Two_Mode_Section1,Indexes_Of_The_Suggestions_From_The_Different_Algorithms1 =  Get_Suggestions_From_Two_Mode_Section(fun,Hyperparameters_1,Hyperparameters_2,Hyperparameters_3,Hyperparameters_4,Suggestion_Tools_1,Suggestion_Tools_2,Suggestion_Tools_3,Suggestion_Tools_4,Hyperparameters_For_Symmetric_Conductor)
     
     
     if len(Cost_Tracker) > 1:
         All_Suggestions_From_Two_Mode_Section1,Indexes_Of_The_Suggestions_From_The_Different_Algorithms1 =  Get_Suggestions_From_Two_Mode_Section2(fun,Hyperparameters_1,Hyperparameters_2,Hyperparameters_3,Hyperparameters_4,Suggestion_Tools_1,Suggestion_Tools_2,Suggestion_Tools_3,Suggestion_Tools_4,Hyperparameters_For_Symmetric_Conductor,Best_Quarter_Of_Individuals)

     Cost = np.zeros(len(All_Suggestions_From_Two_Mode_Section1[1,:]))                                                                                                                         
     ## Sample the Function at the suggested sample points
     for j in range(0,len(All_Suggestions_From_Two_Mode_Section1[1,:])):
         
         if Number_Of_Evaluations > Allowed_Number_Of_Function_Evaluations:

             CurrentCost = PreviousCost
             ForceStop = 1
             break
         Cost[j] =fun(All_Suggestions_From_Two_Mode_Section1[:,j])
         Number_Of_Evaluations +=1
     if ForceStop == 1:
         break


     ## Find the Optimzal Sample and calculate the change vectors                                                                                                                          
     index_min = np.argmin(Cost)
     CurrentCost = Cost[index_min ]
     Cost_Change = CurrentCost - PreviousCost
     Parameter_ChangeVector = np.array([All_Suggestions_From_Two_Mode_Section1[:,index_min] - Best_Individual]).T


     # Update the Addaptive Amplitudes for the next search
     # This should minimize the amplitude of ones that have large variance, and increase those which have high variance.
     Cost_Change_From_Probe_Samples  = Cost - PreviousCost

     Change_Vectors = All_Suggestions_From_Two_Mode_Section1 - np.array([Best_Individual]).T

     Gradients =   Cost_Change_From_Probe_Samples/((np.maximum(Change_Vectors,-0.0001) + np.minimum(Change_Vectors,-0.0001)) + 0.001)

     Update_Addaptive_Amplitudes(Gradients,Hyperparameters_1,Hyperparameters_2,Hyperparameters_3,Hyperparameters_4,Suggestion_Tools_1,Suggestion_Tools_2,Suggestion_Tools_3,Suggestion_Tools_4)
 


      ## Update Resource Allocation          
     Best_Cost_Change_From_Individuals1[:,i] = Update_Resource_Allocation(Cost_Change_From_Probe_Samples,Indexes_Of_The_Suggestions_From_The_Different_Algorithms1,Suggestion_Tools_1,Suggestion_Tools_2,Suggestion_Tools_3,Suggestion_Tools_4,Hyperparameters_For_Symmetric_Conductor)
                                                                                                                   
                                                                                                                                                                                                                                                                                                                                                     
     ## Update Best_Individuals_And_Costs_From_This_Round                                                                                                                                       
     Best_Individual = All_Suggestions_From_Two_Mode_Section1[:,index_min]

     ## Get the Best Quarter of the population
     # This is 1/8th actually
     Best_Quarter_Of_Individuals =  All_Suggestions_From_Two_Mode_Section1[:,np.array(np.nonzero(Cost <= np.median(Cost[np.nonzero(Cost <= np.median(Cost[np.nonzero(Cost <= np.median(Cost))]))])))[0,:]]

     
     Cost[index_min ] = 10000*Cost[index_min ]
     index_min = np.argmin(Cost)
     Second_Best_Individual = All_Suggestions_From_Two_Mode_Section1[:,index_min]
     PreviousCost = CurrentCost                                                                                                                
     Update_Best_Individuals_And_Costs_From_This_Round(Best_Individual,Second_Best_Individual,Parameter_ChangeVector,Cost_Change,Suggestion_Tools_1,Suggestion_Tools_2,Suggestion_Tools_3,Suggestion_Tools_4)                                                                                                                           

     
     Cost_Tracker.append(CurrentCost) 


     ############### Restart Conditions ##########################
     if np.min(np.array([Suggestion_Tools_1.Random_Search_exponential_Growth_Factor,               Suggestion_Tools_2.Random_Search_exponential_Growth_Factor,Suggestion_Tools_3.Random_Search_exponential_Growth_Factor,               Suggestion_Tools_4.Random_Search_exponential_Growth_Factor])) > 30:
         Restart = 1
#             print("Restart Triggered by Random Search Growth Factor")
         break
     
     
     if len(Cost_Tracker) > 100: ## Not Made Much Improvement lately - on absolute scale
         if Cost_Tracker[len(Cost_Tracker)-30] - Cost_Tracker[len(Cost_Tracker)-1] < 0.5*10**-10:
             Restart = 1
#                 print("Restart Triggered by Not Made Much Improvement lately")
             break


     if len(Cost_Tracker) > 100: ## Not Made Much Improvement lately - on relative scale
         if np.abs((Cost_Tracker[len(Cost_Tracker)-30] - Cost_Tracker[len(Cost_Tracker)-1])/Cost_Tracker[0]) < 0.5*10**-10:
             Restart = 1
#                 print("Restart Triggered by Not Made Much Improvement lately - relative scale")
             break
             
     if len(Cost_Tracker) > 100: ## Stagnation
         if Cost_Tracker[len(Cost_Tracker)-30] == Cost_Tracker[len(Cost_Tracker)-1]:
             Restart = 1
#                 print("Restart Triggered by Stagnation")
             break        
     
 
     
     if print_Cost == 1:
         print("Run Summary:")
         print("---------------------------------------------------------------------")
         print("CurrentCost is:",CurrentCost,"After Running Iteration",current_iteration,"\\",Total_Number_Of_Iterations)
         print("the Random_Search_exponential_Growth_Factors are :",Suggestion_Tools_1.Random_Search_exponential_Growth_Factor,               Suggestion_Tools_2.Random_Search_exponential_Growth_Factor,Suggestion_Tools_3.Random_Search_exponential_Growth_Factor,               Suggestion_Tools_4.Random_Search_exponential_Growth_Factor)
         print("the Targeted_Search_exponential_Growth_Factor is:",Suggestion_Tools_1.Targeted_Search_exponential_Growth_Factor   ,              Suggestion_Tools_2.Targeted_Search_exponential_Growth_Factor   ,Suggestion_Tools_3.Targeted_Search_exponential_Growth_Factor   ,              Suggestion_Tools_4.Targeted_Search_exponential_Growth_Factor   )    
         print("---------------------------------------------------------------------")
         print("")            
         
       
 return Best_Individual,Second_Best_Individual,CurrentCost,PreviousCost,Best_Cost_Change_From_Individuals1,current_iteration,Number_Of_Evaluations ,ForceStop,Restart,Cost_Tracker,Best_Quarter_Of_Individuals


# In[8]:



def Force_Suggested_Hyperparameters_To_be_Within_The_Allowed_Limits(Suggested_Hyper_Parameter_Samples_From_Two_Mode,Hyperparameter_Allowed_Values):

## The 0 second index is the lower bound of the allowed values, and the 1 second index is the upper bound

    Suggested_Hyper_Parameter_Samples_From_Two_Mode[0,:] = np.minimum(np.maximum(Suggested_Hyper_Parameter_Samples_From_Two_Mode[0,:],Hyperparameter_Allowed_Values[0,0]),Hyperparameter_Allowed_Values[0,1])
    Suggested_Hyper_Parameter_Samples_From_Two_Mode[1,:] = np.minimum(np.maximum(Suggested_Hyper_Parameter_Samples_From_Two_Mode[1,:],Hyperparameter_Allowed_Values[1,0]),Hyperparameter_Allowed_Values[1,1])
    Suggested_Hyper_Parameter_Samples_From_Two_Mode[2,:] = np.minimum(np.maximum(Suggested_Hyper_Parameter_Samples_From_Two_Mode[2,:],Hyperparameter_Allowed_Values[2,0]),Hyperparameter_Allowed_Values[2,1])
    Suggested_Hyper_Parameter_Samples_From_Two_Mode[3,:] = np.minimum(np.maximum(Suggested_Hyper_Parameter_Samples_From_Two_Mode[3,:],Hyperparameter_Allowed_Values[3,0]),Hyperparameter_Allowed_Values[3,1])
    Suggested_Hyper_Parameter_Samples_From_Two_Mode[4,:] = np.minimum(np.maximum(Suggested_Hyper_Parameter_Samples_From_Two_Mode[4,:],Hyperparameter_Allowed_Values[4,0]),Hyperparameter_Allowed_Values[4,1])
    Suggested_Hyper_Parameter_Samples_From_Two_Mode[5,:] = np.minimum(np.maximum(Suggested_Hyper_Parameter_Samples_From_Two_Mode[5,:],Hyperparameter_Allowed_Values[5,0]),Hyperparameter_Allowed_Values[5,1])
    Suggested_Hyper_Parameter_Samples_From_Two_Mode[6,:] = np.minimum(np.maximum(Suggested_Hyper_Parameter_Samples_From_Two_Mode[6,:],Hyperparameter_Allowed_Values[6,0]),Hyperparameter_Allowed_Values[6,1])
    Suggested_Hyper_Parameter_Samples_From_Two_Mode[7,:] = np.minimum(np.maximum(Suggested_Hyper_Parameter_Samples_From_Two_Mode[7,:],Hyperparameter_Allowed_Values[7,0]),Hyperparameter_Allowed_Values[7,1])

    
    Hyperparameters_1 = Suggested_Hyper_Parameter_Samples_From_Two_Mode[:,0]
    Hyperparameters_2 = Suggested_Hyper_Parameter_Samples_From_Two_Mode[:,1]
    Hyperparameters_3 = Suggested_Hyper_Parameter_Samples_From_Two_Mode[:,2]
    Hyperparameters_4 = Suggested_Hyper_Parameter_Samples_From_Two_Mode[:,3]

    return Hyperparameters_1,Hyperparameters_2,Hyperparameters_3,Hyperparameters_4,Suggested_Hyper_Parameter_Samples_From_Two_Mode



def Run_Optimization_With_TWO_Mode_Online_Hyperparameter_Addaptation(InitialParameters,fun,Total_Number_Of_Iterations,Number_Of_Iterations_With_Hyperparameters_Set,Hyperparameters_1,Hyperparameters_2,Hyperparameters_3,Hyperparameters_4,Suggestion_Tools_1,Suggestion_Tools_2,Suggestion_Tools_3,Suggestion_Tools_4,Hyperparameters_For_Online_Hyperparameter_Search,print_Cost,Total_Resources,Hyperparameter_Allowed_Values,Number_Of_Evaluations,Allowed_Number_Of_Function_Evaluations,ForceStop,Restart  ):

    Cost_Tracker = []
    Number_Of_Iterations_With_These_Hyperparameters = Number_Of_Iterations_With_Hyperparameters_Set

    
    Random_Search_Growth_Rate  = Hyperparameters_For_Online_Hyperparameter_Search[0]
    Random_Search_Period = Hyperparameters_For_Online_Hyperparameter_Search[1]
    Maximal_Random_Search =Hyperparameters_For_Online_Hyperparameter_Search[2]
    Targeted_Search_Growth_Rate = Hyperparameters_For_Online_Hyperparameter_Search[3]
    Targeted_Search_Decay_Rate = Hyperparameters_For_Online_Hyperparameter_Search[4]
    Significant_Change_Value = Hyperparameters_For_Online_Hyperparameter_Search[5]
    alpha = Hyperparameters_For_Online_Hyperparameter_Search[6]
    beta = Hyperparameters_For_Online_Hyperparameter_Search[7]

    ## Initialize the suggestion tools for the hyperparameter optimization
    Random_Search_exponential_Growth_Factor = 0
    Targeted_Search_exponential_Growth_Factor = 0
    Hyper_Parameter_ChangeVector = np.zeros(len(Hyperparameters_1))
    Hyper_Parameter_Cost_Change = 0
    current_iteration = 0
    
    
    ## This is to set the hyperparameters for optimization
    Number_Of_Samples = Number_Of_Algorithms
    Best_Individual_Hyperparameters = Hyperparameters_1
    Second_Best_Individual_Hyperparameters = Hyperparameters_1
    
    ## This is to set the function parameters for optimization    
    Best_Individual_Parameters =     InitialParameters
    Second_Best_Individual_Parameters =     InitialParameters
    PreviousCost = fun(InitialParameters)
    Best_Quarter_Of_Individuals = np.zeros((len(InitialParameters),int(Total_Resources/4)))

    ## Update Directional Pools with random gradients
    Random_Parameter_Gradients = np.random.rand(len(InitialParameters),Total_Resources)
    Random_Parameter_Gradients[:,0] = 0
    Cost_Change_From_Probe_Samples = -np.ones(Total_Resources)

    ## This is find the Cost of the first set of hyperparameters
    Best_Individual_Parameters,Second_Best_Individual_Parameters,CurrentCost,PreviousCost,Best_Cost_Change_From_Individuals1,current_iteration,Number_Of_Evaluations,ForceStop,Restart,Cost_Tracker,Best_Quarter_Of_Individuals  =  Optimize_With_TWO_Mode_Section(Best_Individual_Parameters,Second_Best_Individual_Parameters,PreviousCost,fun,current_iteration,Total_Number_Of_Iterations,Number_Of_Iterations_With_These_Hyperparameters,Hyperparameters_1,Hyperparameters_2,Hyperparameters_3,Hyperparameters_4,Suggestion_Tools_1,Suggestion_Tools_2,Suggestion_Tools_3,Suggestion_Tools_4,Number_Of_Evaluations,Allowed_Number_Of_Function_Evaluations,ForceStop ,Cost_Tracker,Best_Quarter_Of_Individuals)

    Cost_Function_For_Hyperparameter_Optimization =( np.mean(Best_Cost_Change_From_Individuals1,1) + np.max(Best_Cost_Change_From_Individuals1,1))/2
    index_min = np.argmin(Cost_Function_For_Hyperparameter_Optimization)
    Current_Hyperparameter_Cost = Cost_Function_For_Hyperparameter_Optimization[index_min ]
    Previous_Hyperparameter_Cost = Cost_Function_For_Hyperparameter_Optimization[index_min]

    Adaptive_Amplitude = np.ones((len(Hyperparameters_1),Number_Of_Samples))
    Cost_Function_For_Hyperparameter_Optimization = np.zeros(Number_Of_Samples)
    S = 0
    
       
    Number_Of_Hyperparameter_Iterations = np.ceil(Total_Number_Of_Iterations/Number_Of_Iterations_With_These_Hyperparameters)

    for i in range(0,int(Number_Of_Hyperparameter_Iterations)):
        
        if Restart  ==1:
            break

        ## Generate suggested hyperparameters    
        Suggested_Hyper_Parameter_Samples_From_Two_Mode, Random_Search_exponential_Growth_Factor, Targeted_Search_exponential_Growth_Factor =     Two_Mode_Optimizer(Best_Individual_Hyperparameters,Second_Best_Individual_Hyperparameters,Number_Of_Samples,Hyper_Parameter_ChangeVector,Hyper_Parameter_Cost_Change, Random_Search_exponential_Growth_Factor,Targeted_Search_exponential_Growth_Factor,Targeted_Search_Decay_Rate,Random_Search_Growth_Rate,Targeted_Search_Growth_Rate,Random_Search_Period,Significant_Change_Value,Maximal_Random_Search,Adaptive_Amplitude)
        ## Force the hyperparameters to be within the allowed range
        Hyperparameters_1,Hyperparameters_2,Hyperparameters_3,Hyperparameters_4,Suggested_Hyper_Parameter_Samples_From_Two_Mode = Force_Suggested_Hyperparameters_To_be_Within_The_Allowed_Limits(Suggested_Hyper_Parameter_Samples_From_Two_Mode,Hyperparameter_Allowed_Values)
 

        ## Update Directional Pools with random gradients
        Random_Parameter_Gradients = np.random.rand(len(InitialParameters),Total_Resources)
        Random_Parameter_Gradients[:,0] = 0
        Cost_Change_From_Probe_Samples = -np.ones(Total_Resources)

    
        ## Try to optimize with these hyperparameters (Maybe output the suggestion tools as well?)
        Best_Individual_Parameters,Second_Best_Individual_Parameters,CurrentCost,PreviousCost,Best_Cost_Change_From_Individuals1,current_iteration,Number_Of_Evaluations,ForceStop,Restart,Cost_Tracker,Best_Quarter_Of_Individuals  =  Optimize_With_TWO_Mode_Section(Best_Individual_Parameters,Second_Best_Individual_Parameters,PreviousCost,fun,current_iteration,Total_Number_Of_Iterations,Number_Of_Iterations_With_These_Hyperparameters,Hyperparameters_1,Hyperparameters_2,Hyperparameters_3,Hyperparameters_4,Suggestion_Tools_1,Suggestion_Tools_2,Suggestion_Tools_3,Suggestion_Tools_4,Number_Of_Evaluations ,Allowed_Number_Of_Function_Evaluations,ForceStop,Cost_Tracker,Best_Quarter_Of_Individuals)
       
    
        ## Calculate the Cost_Function_For_Hyperparameter_Optimization
        Cost_Function_For_Hyperparameter_Optimization =( np.mean(Best_Cost_Change_From_Individuals1,1) + np.min(Best_Cost_Change_From_Individuals1,1))/2
        Cost_Function_For_Hyperparameter_Optimization1 = np.array([Cost_Function_For_Hyperparameter_Optimization[0],Cost_Function_For_Hyperparameter_Optimization[1],Cost_Function_For_Hyperparameter_Optimization[2],Cost_Function_For_Hyperparameter_Optimization[3]])

        
        ## Find the Best set of hyperparameters and the change vectors
        index_min = np.argmin(Cost_Function_For_Hyperparameter_Optimization)
        Current_Hyperparameter_Cost = Cost_Function_For_Hyperparameter_Optimization[index_min ]
        Hyper_Parameter_Cost_Change = Current_Hyperparameter_Cost - Previous_Hyperparameter_Cost
        Hyper_Parameter_ChangeVector = Suggested_Hyper_Parameter_Samples_From_Two_Mode[:,index_min] - Best_Individual_Hyperparameters


        
        ## This is the RMS Amplitude Addaptation - For the hyperparameter search
        Cost_Change_From_Probe_Samples  = Cost_Function_For_Hyperparameter_Optimization - Cost_Function_For_Hyperparameter_Optimization[0]

        Change_Vectors = Suggested_Hyper_Parameter_Samples_From_Two_Mode - np.array([Best_Individual_Hyperparameters]).T

        Gradient =   Cost_Change_From_Probe_Samples/((np.maximum(Change_Vectors,-0.0001) + np.minimum(Change_Vectors,-0.0001)) + 0.001)

        S = beta*S + (1-beta)*(np.mean(np.abs(Gradient),1))**2
        S = S*np.sign(np.mean(Gradient,1))
        Adaptive_Amplitude1 = alpha/np.sqrt(np.fmax(S + alpha**2,0.1*alpha**2))
        Adaptive_Amplitude[:,:] = np.array([Adaptive_Amplitude1]).T ## this makes an parameter wise adaptive amplitude Based on RMS prob, which we can do element wise multiplication with the random searches.
        # This should minimize the amplitude of ones that have large variance, and increase those which have high variance.

        
        ## Save the Best Sets of hyperparameters and Cost Function
        Best_Individual_Hyperparameters = Suggested_Hyper_Parameter_Samples_From_Two_Mode[:,index_min]
        Cost_Function_For_Hyperparameter_Optimization[index_min ] = (10**10) + Cost_Function_For_Hyperparameter_Optimization[index_min ]
        index_min = np.argmin(Cost_Function_For_Hyperparameter_Optimization)
        Second_Best_Individual_Hyperparameters = Suggested_Hyper_Parameter_Samples_From_Two_Mode[:,index_min]
        

        Previous_Hyperparameter_Cost = Current_Hyperparameter_Cost
        

        
        if print_Cost == 1:
            print("hyperparameter Summary:")
            print("---------------------------------------------------------------------")
            print("After Running hyperparameter Iteration",i,"\\",int(Number_Of_Hyperparameter_Iterations))
            print("Cost_Function_For_Hyperparameter_Optimization",Cost_Function_For_Hyperparameter_Optimization1)          
            print("")
            print("Best_Individual_Hyperparameters",Best_Individual_Hyperparameters)
            print("")
            print("the Random_Search_exponential_Growth_Factor for the hyperparameter search is:",Random_Search_exponential_Growth_Factor)
            print("the Targeted_Search_exponential_Growth_Factor for the hyperparameter search is:",Targeted_Search_exponential_Growth_Factor)        
            print("---------------------------------------------------------------------")
            print("")            
            

    return Best_Individual_Parameters,CurrentCost,Number_Of_Evaluations ,ForceStop,Restart


# In[9]:


def Set_Hyperparameter_Ranges(Random_Search_Growth_Rate_Range, Random_Search_Period_Range   ,Maximal_Random_Search_Range,Targeted_Search_Growth_Rate_Range ,Targeted_Search_Decay_Rate_Range ,Significant_Change_Value_Range,alpha_Range ,beta_Range):
    Hyperparameter_Allowed_Values = np.zeros((11,2))

    Hyperparameter_Allowed_Values[0,:] = Random_Search_Growth_Rate_Range
    Hyperparameter_Allowed_Values[1,:] = Random_Search_Period_Range    
    Hyperparameter_Allowed_Values[2,:] = Maximal_Random_Search_Range    
    Hyperparameter_Allowed_Values[3,:] = Targeted_Search_Growth_Rate_Range    
    Hyperparameter_Allowed_Values[4,:] = Targeted_Search_Decay_Rate_Range   
    Hyperparameter_Allowed_Values[5,:] = Significant_Change_Value_Range
    Hyperparameter_Allowed_Values[6,:] = alpha_Range
    Hyperparameter_Allowed_Values[7,:] = beta_Range

    return Hyperparameter_Allowed_Values


def Mean_ChangeValue_fit(Maximal_Random_Search,Number_Of_Parameters):
    fit_Parameters = [-0.007921662,0.003755446,5.92*(10**(-5)),-2.45*(10**(-7)),0.250014494] # These were found numerically from the file "Finding_The_Segnificant_Change_Value.ipynb"
    return    fit_Parameters[0]*Maximal_Random_Search + fit_Parameters[1]*Maximal_Random_Search**2 + fit_Parameters[2]*Number_Of_Parameters + fit_Parameters[3]*Number_Of_Parameters**2 + fit_Parameters[4]*Maximal_Random_Search*Number_Of_Parameters
    


# In[10]:


def Orchestra_Optimize_With_Two_Mode_Section(x0,fun,Allowed_Number_Of_Function_Evaluations,print_Cost):

    Number_Of_Iterations = int(np.maximum(np.ceil(Allowed_Number_Of_Function_Evaluations/60) -2,1))

    Number_Of_Parameters = len(x0)
    InitialParameters = x0
    Total_Number_Of_Iterations = Number_Of_Iterations
    
    Number_Of_Iterations_With_Hyperparameters_Set = 5
    
    

    ############################## Setting up The Conductor (Resource allocation) algorithm ###################################
    Total_Resources = 64
    Number_Of_Algorithms = 4
    Mass =10
    Self_Spring_Constants = 10
    Interaction_Spring_Constants = 40

    Hyperparameters_For_Symmetric_Conductor = Set_Hyperparameters_For_Symmetric_Conductor(Total_Resources,Number_Of_Algorithms,Mass,Self_Spring_Constants,Interaction_Spring_Constants)
    Initial_Resource_Allocation = Hyperparameters_For_Symmetric_Conductor[:,0]
    Current_Resource_Allocation = Hyperparameters_For_Symmetric_Conductor[:,0]

    Changevector = np.array([0, 0 ,0 ,0]) ## Negative is bad
    Recomended_Resource_Allocation = Symmetric_Section_Conductor(Current_Resource_Allocation,Changevector,Hyperparameters_For_Symmetric_Conductor) 
    Current_Resource_Allocation = Recomended_Resource_Allocation

    #############################################################################################################################




    
    ######################### Setting the hyperparameter ranges ######################################################################
    Random_Search_Growth_Rate_Range   = np.array([0.01,10])
    Random_Search_Period_Range  = np.array([0.1,10])
    Maximal_Random_Search_Range  = np.array([0.0000001,10])
    Targeted_Search_Growth_Rate_Range  = np.array([0.01,5])
    Targeted_Search_Decay_Rate_Range  = np.array([0.1,100]) 
    Significant_Change_Value_Range  = np.array([10**-30,0.01]) 
    alpha_Range = np.array([10**-5,10**15]) ## (Here Do 10**random)
    beta_Range = np.array([0.1,1]) 
    Directional_Pool_Ratio_Range = np.array([0.0001,1]) 

    Hyperparameter_Allowed_Values = Set_Hyperparameter_Ranges(Random_Search_Growth_Rate_Range, Random_Search_Period_Range   ,Maximal_Random_Search_Range,Targeted_Search_Growth_Rate_Range ,Targeted_Search_Decay_Rate_Range ,Significant_Change_Value_Range,alpha_Range ,beta_Range )
    #############################################################################################################################



    ############################### Setting the Hyperparameters_For_Online_Hyperparameter_improvement  Algorithm #########################
    Hyperparameters_For_Online_Hyperparameter_Search = np.zeros(8)
    Hyperparameters_For_Online_Hyperparameter_Search[0] =5.0   # Random_Search_Growth_Rate 
    Hyperparameters_For_Online_Hyperparameter_Search[1] =5.0    # Random_Search_Period 
    Hyperparameters_For_Online_Hyperparameter_Search[2] =20.02   # Maximal_Random_Search 
    Hyperparameters_For_Online_Hyperparameter_Search[3] =0.5    # Targeted_Search_Growth_Rate
    Hyperparameters_For_Online_Hyperparameter_Search[4] =2.0    # Targeted_Search_Decay_Rate 
    Hyperparameters_For_Online_Hyperparameter_Search[5] =0.00001*Mean_ChangeValue_fit(Hyperparameters_For_Online_Hyperparameter_Search[2] ,8)    # Significant_Change_Value
    Hyperparameters_For_Online_Hyperparameter_Search[6] =(10**4)/Hyperparameters_For_Online_Hyperparameter_Search[2]    # alpha 
    Hyperparameters_For_Online_Hyperparameter_Search[7] =0.6     # beta

    #############################################################################################################################


    
    
    
    Number_Of_Evaluations  = 0
    ForceStop = 0
    Restart = 0
    Best_Cost = 10**100
    Number_Of_Restarts = 0
    while ForceStop ==0:    
        ########################initialize the different Suggestion_Tools For the Two Mode and the hyperparameters #######################
        Suggestion_Tools_1 = Create_New_Two_Mode_Suggestion_Tools_Object(Number_Of_Parameters)
        Suggestion_Tools_2 = Create_New_Two_Mode_Suggestion_Tools_Object(Number_Of_Parameters)
        Suggestion_Tools_3 = Create_New_Two_Mode_Suggestion_Tools_Object(Number_Of_Parameters)
        Suggestion_Tools_4 = Create_New_Two_Mode_Suggestion_Tools_Object(Number_Of_Parameters)

        Suggestion_Tools_1.Number_Of_Samples =      Current_Resource_Allocation[0]                                                                                                                     
        Suggestion_Tools_2.Number_Of_Samples =      Current_Resource_Allocation[1] 
        Suggestion_Tools_3.Number_Of_Samples =      Current_Resource_Allocation[2]                                                                                                                                  
        Suggestion_Tools_4.Number_Of_Samples =      Current_Resource_Allocation[3]  



        ## This is to initiailize the Hyperparameters - Note that they should all be positive numbers (Maybe make 6 and 7 be 10**rand(2))
        Initial_Hyper_Parameter_Samples = 10*np.random.rand(11,Number_Of_Algorithms)
        Hyperparameters_1,Hyperparameters_2,Hyperparameters_3,Hyperparameters_4,Initial_Hyper_Parameter_Samples = Force_Suggested_Hyperparameters_To_be_Within_The_Allowed_Limits(Initial_Hyper_Parameter_Samples,Hyperparameter_Allowed_Values)
        #############################################################################################################################





        Best_Parameters1,Best_Cost1,Number_Of_Evaluations,ForceStop,Restart  = Run_Optimization_With_TWO_Mode_Online_Hyperparameter_Addaptation(InitialParameters,fun,Total_Number_Of_Iterations,Number_Of_Iterations_With_Hyperparameters_Set,Hyperparameters_1,Hyperparameters_2,Hyperparameters_3,Hyperparameters_4,Suggestion_Tools_1,Suggestion_Tools_2,Suggestion_Tools_3,Suggestion_Tools_4,Hyperparameters_For_Online_Hyperparameter_Search,print_Cost,Total_Resources,Hyperparameter_Allowed_Values,Number_Of_Evaluations,Allowed_Number_Of_Function_Evaluations,ForceStop,Restart )
        

        if Restart == 1:
            Number_Of_Restarts +=1
#             print(Number_Of_Restarts)
#             print(Best_Cost1)

        #if np.mod(Number_Of_Restarts,2) ==1 :
        InitialParameters = Best_Parameters1

      #  if np.mod(Number_Of_Restarts,2) ==0 :
      #      InitialParameters = 10*(np.random.rand(len(Best_Parameters1)) - 0.5)
            
        Restart = 0
        if Best_Cost1 <= Best_Cost:

            Best_Cost = Best_Cost1
            Best_Parameters = Best_Parameters1


    return Best_Parameters


# In[11]:

from multiprocessing import Process

def fun(x): ## Ellispod
#     return np.sum((x-2)**2) + np.sum(np.sin(x))

    a = 1**(np.arange(D)/D)
    return np.sum(a*((x-X_Optimal)**2) )


D = 320
Number_Of_Function_Evaluations = 3000*D
Number_Of_Tests = 11
print_Cost = 0
Print_Two_Mode_Analytics = 0
test_Optimization = 0
if test_Optimization == 1:

    for p in range (Number_Of_Tests):
        x0 = 10*(np.random.rand(D)-0.5)
        X_Optimal = 10*(np.random.rand(D)-0.5)
        print("Initial Cost",fun(x0))
        Best_Parameters_After_Optimization = Orchestra_Optimize_With_Two_Mode_Section(x0,fun,Number_Of_Function_Evaluations,print_Cost)
        print("Final Cost",fun(Best_Parameters_After_Optimization))
        print()


# In[6]:


Test_Optimization = 0
if Test_Optimization == 1:



    from multiprocessing import Process

    def Run_Sphere_Test(D,Budget,Function_Number,Process_Number):
        Number_Of_Function_Evaluations = Budget*D
        x0 = 10*(np.random.rand(D)-0.5)
        X_Optimal = 10*(np.random.rand(D)-0.5)

        if Function_Number == 1:
            def fun(x): ## Sphere
                return np.sum(((x-X_Optimal)**2) )

        if Function_Number == 2:
            def fun(x): ## Ellispod
                a = 1000000**(np.arange(D)/D)
                return np.sum(a*((x-X_Optimal)**2) )

        if Function_Number == 3:
            def fun(x): ## Sum_Of_Powers
                result = 0
                for i in range(0,len(x)):
                    result +=     (np.abs(x[i]-X_Optimal[i]))**(2+4*(i/len(x)))    
                return result


        print("Initial Cost",fun(x0),"   Process Number  ",Process_Number)
        Best_Parameters_After_Optimization = Orchestra_Optimize_With_Two_Mode_Section(x0,fun,Number_Of_Function_Evaluations,print_Cost)
        print("Final Cost",fun(Best_Parameters_After_Optimization),"   Process Number  ",Process_Number)
        print()


    D = 5
    Budget = 3000
    Function_Number = 1
    Number_Of_Cores = 1
    Number_Of_Tests = 4


    if Number_Of_Cores == 1:
        for test in range (Number_Of_Tests):
            Run_Sphere_Test(D,Budget,Function_Number,test)

    if Number_Of_Cores == 4:
        if __name__ == '__main__':
            p1 = Process(target = Run_Sphere_Test, args=(D,Budget,Function_Number,1))
            p2 = Process(target = Run_Sphere_Test, args=(D,Budget,Function_Number,2))
            p3 = Process(target = Run_Sphere_Test, args=(D,Budget,Function_Number,3))
            p4 = Process(target = Run_Sphere_Test, args=(D,Budget,Function_Number,4))
            p1.start()
            p2.start()
            p3.start()
            p4.start()
            p1.join()
            p2.join()
            p3.join()
            p4.join()


# In[8]:




## Defining Cost Function
class Adaptive_Two_Mode:
    
    Name = 'Adaptive_Two_Mode'

    def __init__(self, Name):
        self.__name__ = 'Adaptive_Two_Mode'

        
    def Optimize(fun, x0, remaining_evals ): 


        Best_Parameters_After_Optimization = Orchestra_Optimize_With_Two_Mode_Section(x0,fun,remaining_evals,print_Cost)
        
  

        return Best_Parameters_After_Optimization
    
    
    
    
    
SOLVER = Adaptive_Two_Mode

print("Finished Loading: ", SOLVER.__name__)
print()
#print(SOLVER.__name__ == 'Adaptive_Two_Mode')

#def fun(x):
#    return np.sum((x-2)**2)

#D  = 10
#x0= np.random.rand(D)
# remaining_evals = 1000*D
#Number_Of_Function_Evaluations = 60*D

#print(fun(x0))
#x_min1  = SOLVER.Optimize(fun, x0,Number_Of_Function_Evaluations)
#print(fun(x_min1))


# In[13]:


def optimize(fun, x0, remaining_evals ): 


    Best_Parameters_After_Optimization = Orchestra_Optimize_With_Two_Mode_Section(x0,fun,remaining_evals,print_Cost)
    


    return Best_Parameters_After_Optimization

