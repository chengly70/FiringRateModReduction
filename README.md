# FiringRateMod_Reduction

This directory contains Matlab program files to implement our reduction method for coupled firing rate neural network models subject to correlated background noise.

Citation: Barreiro & Ly (2017). A Practical Approximation Method for Firing Rate Models of Coupled Neural Networks with Correlated Inputs. Physical Review E xx, xx
DOI: 

All implementations of the method uses a standard sigmoidal transfer function (described in the paper), users can augment this to other transfer functions by manually changing our code.

Sub-directories:
LargerNetwork/ — contains an alternative method for larger sized (> ~100) networks where only a specified subset of covariance entries are saved. Contains 2 files: iter_method_large.m (our reduction method) and mc_WC_large.m (analogous Monte Carlo simulation)
LowOrder_FP/ — implements the ‘Lowest Order Approximation’ method based on a moment closure reduction, described in the Appendix of our paper. Only has the reduction method function iterLowO_Ref.m ; to compare to Monte Carlo, see mc_WC.m below

The main method described in the paper above is in the function file:
iter_method.m — see comments for Inputs/Outputs.  

mc_WC.m — implements Monte Carlo simulations to compare with the results of iter_method.m

We include several driver files that run the iter_method.m function and saves the results in a corresponding .mat file.
To reproduce the plots in the paper (with different instances of random/quenched parameters), first run scptAn_[].m THEN run scptMC_[].m

1) Network with 2 neurons, run scptAn_2nrn.m; this relies on d2FF_b.mat file (parameter file) saves a mat file dAn2_A.mat after finished.
To run/save Monte Carlo of the same network, run scptMC_2nrn.m
These results are shown in Fig 1 of our paper above

2) Network II, Fig 2.
run scptAn_netII.m first, saves dAn_netTw_aut.mat. To compare with Monte Carlo, run scptMC_netII.m AFTER running scptAn_netII.m (saves results in dmc_netTw_aut.mat)

3) Network III, Figs 3.
run scptAn_netTr.m first, saves dAn_netEI.mat. To compare with Monte Carlo, run scptMC_netTr.m AFTER running scptAn_netTr.m, insuring method converges, etc. (save results in dmc_netEI.mat)

4) Figure 5 in paper with disparate tau values.
scptAn_testTau.m and scptMC_testTau.m (see files for corresponding .mat files saved)

5) Figure 7 in paper. Relies on method in sub-directory ~/LowOrder_FP/ iterLowO_Ref.m. Run this function as you would run the function iter_method.m, following the structure of the files scptAn_[].m — to compare against a particular instance of a network, use the file scptMC_testTau.m as a guide (load the particular network parameters first).

Figures 4 & 6 are self-explanatory once you understand the program files.

VISUALIZATION:
Including one script file plots_netII.m to generate figures (analytic method & Monte Carlo).  To run this file without new instances, download dAn_netTw_aut.mat and dmc_netTw_aut.mat.  Augment this file to get other plots.


