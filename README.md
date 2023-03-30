# Network-Impact-of-a-single-sample
Python functions for the calculation of the network impact parameters.

The "NetworkImpactPlots.py" script contains two general function used for both setups:

1) Network_Impact_parameters_Calculate(reference_cohort, test_set) - used to calculate the network impact parameters given the reference cohort (m samples) and the test cohort, which should contain all samples from the reference cohort and one tested sample (m+1 samples).
2) create_net_intersected(mat, sig_p_val=1e-3) - used to create a network from a given cohort ("mat") with a certain significance threshold (p-value threshold = sig_p_val).

And two setup-specific functions:

1) semi_supervised_plots(file_path) - which plots the semi supervised Fig 3.(a-c),
  This function takes a .xlsx file containing three sheets where each one contains a cohort of samples: 
                                                                                                        a) Reference cohort
                                                                                                        b) GLV Steady State cohort
                                                                                                        c) Shuffled cohort
                                                                                                        
  
  the file path on your computer to the following function
3)

You may use the following files (found in this repository):
1) Semi supervised example cohorts - m=100, 100 test samples.xlsx
2) Supervised example cohorts - m=50, 100 test samples.xlsx
which contains cohorts of simulated samples, configured to suit the appropriate setup (semi supervised and supervised).
