
				README
				======
Date: 
27/08/2015

This is the readme for models associated with the paper:

Alex Pavlides, S. John Hogan and Rafal Bogacz. (2015) Computational
models describing possible mechanisms for generation of excessive beta
oscillations in Parkinson's disease. PLOS Computational Biology

Description: 
This set of Matlab files produce a figure showing Results of simulations of 
the resonance model (panels A and B) and feedback model (panels C and D). A, C) 
Each of the six panels shows the activity of the STN, GPe and cortical populations 
as a function of time. The labels to the left indicate if a row shows simulations 
of an intact model, or a model with particular connections blocked. 

In simulations in panel A the following parameters were used: 

wSG =4.87, wGS=1.33, wCS=9.98, wSC=8.93, wGG=0.53, wCC=6.17, C=172.18, Str=8.46, 
TCC=4.65, tau_E=11.59, tau_I=13.02, B_E=17.85, B_I=9.87, M_E=75.77 and M_I=205.72. 

In simulations in panel C the following values were used: 

wSG=2.56, wGS=3.22, wCS=6.60, wSC= 0.00, wGG=0.90, wCC=3.08, C=277.94, Str=40.51, 
TCC=7.74, tau_E=11.69, tau_I=10.45, B_E=3.62, B_I=7.18, M_E=71.77 and M_I=276.39. 

B, D) The comparison between experimental and simulated statistics of the oscillations.

File list:
main.m
minfunction.m
model_eqs.m
generate_fig.m
frenquency.m
peakdet.m
rotateXLabels.m

Data:
weights.mat


Operating instructions:
Ensure files and data are availabe in Matlab. To run type 'main' at the command line. 
The figure should be generated. 


Contact: If there are any questions please contact Dr Rafal Bogacz by email: rafal.bogacz@ndcn.ox.ac.uk
