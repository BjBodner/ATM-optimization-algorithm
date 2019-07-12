# ATM-optimization-algorithm
An optimization algorithm for high dimensional problems, such as those found in machine learning and physical sciences.

the ATM algorithm is designed to minimize functions which input a vector or parameters and output a single scalar value. (such as cost functions in machine learning).
Examples of how to use call the optimizer are in the file titled : 
" Example_Of_Calling_Adaptive_To_Mode_Optimizer.py".
We invite users to try to use this optimizer on different optimization problems and find as many issues as you can. This will help us improve the optimizer in future versions.

---



## Scalability Experiment
This is a test to demonstrate the linear scaling capabilities of the ATM algorithms.
That is the time to run the algorithm (as well as the memory required and the number of operations required) to
run a fixed number of function evaluations - scales linearly with the size of the search space.


To run a scalability experiment on your PC You can use the Jupyter Notebook titled "Timing_Experiment_ATM.ipynb" .
It will call the file "Adaptive_Two_Mode_Optimizer.py" which contains the ATM algorithm V1.0, which was benchmarked on the BBOB Noiseless testbed in April 2019.

If you would prefer to run it on Google Colab, you can use "Timing_Experiment_For_Colab.ipynb".
It will ask you to upload the file "Adaptive_Two_Mode_Optimizer_For_Colab.py", which is the same version ( ATM algorithm V1.0) but does not import the COCOEX library. After the file is uploaded the experiment should run smoothly

Note that the full experiment can take up to 15 minute to run.

Please let me know if there are any issues with this test .

---

## Goals for future versions:
1. More customization: ability to set the number of iterations, set a threshold for convergence, number of parallelized function evaluations and more..
2. Examples of optimizing nerual networks using numpy
2. Examples of optimizing nerual networks using PyTorch
3. Examples of optimizing nerual networks using Tensorflow


Enjoy!

For questions and comments feel free to contact Benjy at:
benjamin_bodner@brown.edu
