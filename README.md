# high-dimension-sircovid-MRes
This is the repositoray for the MRes project. We also have a test repo at https://github.com/WilliamLiu666/high-dimension-sircovid-test-repo

Project name: Improving real-time inference for high dimensional models\
Institution: Imperial College London\
Supervisor: Dr. Marc Baguelin and Dr. Edward Knock\
Student: Weizhi Liu


Project Aims\
   • Explore and quantify the effect of high dimensionality on the geometry of probability distributions appearing while performing inference in complex models\
   • Quantify the benefits of using fitting algorithms using gradient based samplers in particular Hamiltonian Monte Carlo\
   • Implement and tune a Hamiltonian Monte Carlo algorithm for a complex model\
   • Understand how the performance of random walk MCMC scale with dimensionality and the advantage of gradient based methods

Research Plan\
   • Initially the student will study our existing model and code, learning how to run the fitting of the model to the data.\
   • They will then adapt the code themselves to generate visualisations and summary statistics of the geometric structure of the posterior distributions of the model parameters\
   • Then the student will Implement the Hamiltonian Monte Carlo algorithm for a simple version of the model for an early time point in the pandemic and compare the performance of HMC with the existing adaptive MCMC algorithm using commonly used comparator for performance such as Effective Sample Size (ESS) per computational time.\
   • Finally the student will adapt the Hamiltonian Monte Carlo algorithm for the latest version of the model and compare the performance of HMC with the existing adaptive MCMC algorithm using commonly used comparator for performance such as Effective Sample Size (ESS) per computational time to see how both methods scale with dimension.

*Read and understand about RW-MCMC and Hamiltonian MCMC\
*Run our pipeline (from 1 region, using published repo?)\
*Get familiar with the RW-MCMC algorithm we are using (develop start mode, package and object)\
*Run simple HMC code on simple distribution\
*Runs gradient code for the pipeline\
*Way of parallelising the calculation of gradient?\
*Tuning the HMC\
*Compare algorithms performance


Things to test:\
-precision of numerical gradient; v-shape\
-ESS per computation time for RWM vs HMC\
-Complexity of the algorithm, how they scale with dimension as a function of likelihood & gradient calculation
