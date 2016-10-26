% Script for the creation of an SBmodel of the Novak-Tyson cell cycle 
% model, described in J. theor. Biol. (1998) 195, 69-85
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sbm = SBmodel('NovakTysonModel.txt')

disp('The variable containing the Novak-Tyson model is: sbm');