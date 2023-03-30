# Network-Impact-of-a-single-sample
Python functions for the calculation of the network impact parameters.

The "NetworkImpactPlots.py" script contains two general function used for both setups:

1) "Network_Impact_parameters_Calculate(reference_cohort, test_set)" - used to calculate the network impact parameters given the reference cohort (m samples) and the test cohort, which should contain all samples from the reference cohort and one tested sample (m+1 samples).
2) "create_net_intersected(mat, sig_p_val=1e-3)" - used to create a network from a given cohort ("mat") with a certain significance threshold (p-value threshold = sig_p_val).

And two setup-specific functions, which takes a path of an .xlsx file  (on your computer) containing the relevant cohorts of samples. Each sheet in the file should 
contain a cohort of m (or (m+1)) samples and n species where the samples are stores in the rows, suth that the cohort is of size m*n (or (m+1)*n). The sheets in the file which contans the cohort should have the respective name (mentioned below): 

1) "semi_supervised_plots(file_path)" - plots the semi supervised Fig 3.(a-c).
  This function takes a path to an .xlsx file containing three sheets. Each sheet contains a cohort of samples (no headers) with the following names: 
        a) "Reference cohort" - contains all samples in the reference cohort.
        b) "GLV Steady State cohort" - contains all GLV steady states used as test samples.
        c) "Shuffled cohort" - contains all the shuffled profiles used as test samples.                                                    
2) supervised_plots(file_path) - plot the supervised Fig 5.(a-d).
  This function takes a path to an .xlsx file containing four sheets. Each sheet contains a cohort of samples (no headers) with the following names: 
        a) "Reference cohort A" - contains all the samples in reference cohort of model A (GLV or otherwise).
        b) "Test cohort A" - contains all samples used as test samples for reference cohort A.
        c) "Reference cohort B" - contains all the samples in reference cohort of model B (GLV or otherwise).
        d) "Test cohort B" - contains all samples used as test samples for reference cohort B.
  

You may use your own data or use the following files (found in this repository):

1) Semi supervised example cohorts - m=100, 100 test samples.xlsx
2) Supervised example cohorts - m=50, 100 test samples.xlsx

which contains cohorts of simulated samples, configured to suit the appropriate setup (semi supervised and supervised).
