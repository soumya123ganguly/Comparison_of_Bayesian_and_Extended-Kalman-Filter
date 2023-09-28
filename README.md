# Comparison_of_Bayesian_and_Extended-Kalman-Filter
We write MATLAB codes to run an approximate Bayesian filter and compare its conditional mean and variance with an Extended Kalman filter. 

# Tasks performed : 
For a detailed report please consult the file 'project report.pdf'. For reference to equations/variables used to describe the tasks, please consult 'Problem Formulation.pdf'.
1. Task 1: Use Bayes’ Rule to derive the Bayesian filter.
   
2. Task 2: Write a MATLAB program to run an approximate Bayesian filter for given values of (d; h) (to see what these things mean, please consult 'Problem Formulation.pdf'. This program will
need to be cognizant of:
(i) The range of feasible values for x_t as a function of d and t.
(ii) The range of feasible values for y_t as a function of d; t and h.
(iii) The number of sample points, or density of samples, that you want for each pdf. This will dictate computational complexity. But it is a simple problem.

3. Task 3: Use system equations (1)-(2) to produce a set of sample state and output values. You will test your estimator against these output measurements yt and against these truth samples of xt. Note that MATLAB’s rand function returns a white U[0; 1] set of random variables.

4. Task 4: Run your estimator to determine how well you can estimate x_t for specific values of d, h and t. Note the Bayesian filter returns a conditional density. You will need to develop your own figure of merit for the performance of the filter. Compare the predicted and filtered conditional densities.

5. Task 5: Playing with d and h to explore the behavior of the Bayesian filter in different noise environments. 

6. Task 6: Compare the Bayesian filter to the Extended Kalman Filter in terms of conditional mean and conditional variance.


# How to run the code to obtain the results 
1. You may run 'baysianfilter.m' to see by a plot, what a Bayesian filter predicts and compare with the true state value.
2. Run 'code.m'
3. Input a number as the range of the uniform distribution of disturbance,
4. Input a number as the range of the uniform distribution of noise,
5. Obtain the plots (needed in the tasks) and the predicted state for the Bayesian filter. 
