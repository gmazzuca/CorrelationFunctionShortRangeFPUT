# Correlation_function_Short_range_FPUT
Python script to numerically compute and plot the correlation functions for the periodic FPUT chain

In short_range_FPUT_correlation.py you will find the main program that numerically compute the correlation functions for the FPUT chain, i.e. an Hamiltonian system with short range interaction with Hamiltonian of the form:

H = \sum_{j=0}^{N-1} (p_j^2/2 + \sum_{s=1}^m k_s((q_{j+s} - q_j)^2/2 + chi*(q_{j+s} - q_j)^3/3 + gamma*(q_{j+s} - q_j)^4/4) ))

The script is thought to be used on a cluster in order to produce a lot of data that one has to process and then plot using the python script Scaling_Non_Linear_chain.py.
This last file produce a series of plot of the correlation functions.

Both files contains a brief description of what they do and the parameters that they need.

The script LongRangeOscillatorChain.py instead produce a series of plot of the correlations functions for some short range harmonic chain, i.e. hamiltonian system of the form

H = \sum_{j=0}^{N-1} (p_j^2/2 + \sum_{s=1}^m k_s((q_{j+s} - q_j)^2/2) 

The plots are of three different kind:

- simple correlation functions;
- scaling of the extream peak of the correlation;
- scaling of the central peak of the correlation;

