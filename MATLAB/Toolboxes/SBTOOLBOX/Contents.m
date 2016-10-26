% Systems Biology Toolbox
% Version 1.7 (>=R14SP1) 29-January-2007
%
% DOWNLOAD & DOCUMENTATION 
% ========================
% The toolbox can be downloaded from http://www.sbtoolbox.org
% Extensive documentation and examples are available at the same address.
% In case of questions, suggestions, etc. please contact me at
% henning@fcc.chalmers.se.
%
% Installation of Toolbox
% =======================
%   installSB           - Installation script for the SBtoolbox. Edit the 
%                         script to match your system and run it
%   setvarSB            - Definition of variables used across different
%                         functions of the toolbox
%                       
% Model creation and manipulation
% ===============================
%   SBmodel             - Creating a new SBmodel
%   SBstruct            - Returns the internal data structure of an SBmodel
%   SBedit              - Graphical user interface for editing SBmodels in 
%                         an ODE representation
%   SBeditBC            - Graphical user interface for editing SBmodels in
%                         a more biochemically oriented representation
%
% Data creation and manipulation
% ==============================
%   SBdata                  - Creating a new SBdata data object
%   SBstruct                - Returns the internal data structure of an SBdata
%                             object
%   SBgetdata               - Allows to extract information about the
%                             measurement data stored in an SBdata object
%   SBexportCSVdata         - Exporting an SBdata object to a CSV file
%   SBexportXLSdata         - Exporting an SBdata object to an Excel file 
%   SBvisualizedata         - Visualizing data in an SBdata object graphically
%   SBgetsamplingtimedata   - Determining the sampling times used in a given
%                             SBdata object
%
% Model information
% =================
%   SBstates            - Returns information about states in an SBmodel
%                         (statenames, a cell-array with names of states in model
%   SBinitialconditions - Sets or returns initial conditions of the states
%                         in the model
%   SBparameters        - Returns parameter names and values in an SBmodel
%                         or ODE file model. Also used to change parameter
%                         values
%   SBvariables         - Returns information about variables in an
%                         SBmodel (variable names and formulas, but 
%                         also the variable values for given state)
%   SBreactions         - Returns information about reactions in an SBmodel
%                         (reaction names and formulas of kinetic laws, but 
%                         also the reaction rates for given state)
%   SBfunctions         - Returns information about functions in an SBmodel
%                         (functions names, arguments, and formulas)
%   SBevents            - Returns information about events in an SBmodel
%                         (names,triggers,assignment variables, and assignment 
%                         formulas)
%   SBfunctionsMATLAB   - Returns information about MATLAB functions in an
%                         SBmodel
%
% Export of SBmodel
% =================
%   SBcreateODEfile     - Converting an SBmodel to an ODE file
%   SBcreateTempODEfile - Same as SBcreateODEfile but ODE file is created
%                         in the systems temporary directory
%   deleteTempODEfileSB - Deletes the temporary ODE file 
%   SBcreateXPPAUTfile  - Converting an SBmodel to an XPPAUT ODE file
%   SBcreateTEXTfile    - Converting an SBmodel to a ODE text file description
%   SBcreateTEXTBCfile  - Converting an SBmodel to a biochemical oriented text 
%                         file description
%   SBexportSBML        - Exporting SBmodel to SBML Level 2 Version 1
%   SBconvert2MA        - Converting an SBmodel only containing reactions
%                         with mass action kinetics to a structure
%                         containing information about stoichiometry,
%                         kinetic parameters, and initial conditions
%
% Simulation Functions
% ====================
%   SBsimulate              - Deterministic simulation of an SBmodel or an 
%                             ODE file
%   SBexperiment            - Allows to do in silico experiments
%   SBstochsim              - Stochastic simulation of SBmodels, only
%                             containing reactions with mass action
%                             kinetics (currently only working on windows
%                             systems)
% 
% Plotting Functions
% ==================
%   SBplot              - (GUI) Plots simulation data (e.g., returned 
%                         by SBsimulate)
%   SBplot2             - (GUI) Plots different kind of data where a bar
%                         diagram representation is useful. So far mainly
%                         used for displaying results from parameter
%                         sensitivity analysis 
%
% Simple Analysis Functions
% =========================
%   SBsteadystate           - Determines the steady-state of an SBmodel or an 
%                             ODE file model, dealing also with singular
%                             systems
%   SBjacobian              - Determines the Jacobian of an SBmodel or an ODE 
%                             file
%   SBmoietyconservations   - Determines the moitey conservations and/or other 
%                             conservations that are present in a model
%   SBstoichiometry         - Determines the stoichiometric matrix for the 
%                             given model
%   SBmakeirreversible      - Converting all reversible reactions in an
%                             SBmodel to irreversible ones
% 
% Identification Functions
% ========================
%   SBparameterestimation   - Allows to estimate parameter values for given
%                             measurements and model structure
%   SBnetworkident          - Allows to identify a local network connectivity
%                             matrix from measurements of involved components 
%
% Model Reduction Functions
% =========================
%   SBreducemodel       - Reduces a singular model to a non-singular by 
%                         deleting algebraic realtions
%
% Bifurcation analysis
% ====================
%   SBxppaut            - Starts XPPAUT with the given XPPAUT ODE file
%   SBplotxppaut        - Plots bifurcation data file saved from AUTO/XPPAUT
%
% Parameter Sensitivity Analysis
% ==============================
%   SBsensdataosc       - Generating data for the parameter sensitivity 
%                         analysis of oscillating systems
%   SBsensdataoscevents - Generating data for the parameter sensitivity 
%                         analysis of oscillating systems in the case that
%                         events are present in the model
%   SBsensamplitude     - Parameter sensitivity analysis of the oscillation
%                         amplitude. Uses data generated by SBsensdataosc
%   SBsensperiod        - Parameter sensitivity analysis of the oscillation
%                         period. Uses data generated by SBsensdataosc
%   SBsensdatastat      - Generating data for the parameter sensitivity 
%                         analysis of the steady-state of systems
%   SBsensstat          - Parameter sensitivity analysis of the
%                         steady-state values of states, variables, and 
%                         reaction rates (can be seen as a generalized MCA)
%   SBmca               - Metabolic Control Analysis (MCA). Function
%                         calculating Flux Control Coefficients,
%                         Concentration Control Coefficients, and
%                         Elasticity Coefficients 
%
% Localization of mechanisms leading to complex behaviors
% =======================================================
%   SBlocbehavcomp      - Determines the importance of components in the 
%                         given biochemical system in the creation of an 
%                         observed complex behavior such as multiple 
%                         steady-states and sustained oscillations.
%   SBlocbehavinteract  - Determines the importance of direct interactions
%                         between components in the given biochemical system 
%                         in the creation of an observed complex behavior 
%                         such as multiple steady-states and sustained 
%                         oscillations.
%   SBlocbehavinteract2 - In principle the same as SBlocbehavinteract, but 
%                         possible to use for open-loop unstable systems.
%                         See help text for more information.
%
% Optimization and nonlinear solver functions
% ===========================================
%   simplexSB           - Local minimization function using downhill
%                         simplex method (Nelder-Mead)
%   simannealingSB      - Global minimization function based on simulated
%                         annealing
%   fsolveSB            - Solver for nonlinear equations
%
% String handling functions
% =========================
%   explodePCSB         - auxiliary function allowing to decompose a
%                         string expression into separated elements. The separation 
%                         character can be specified. Commas within parentheses 
%                         expressions are not considered
%   extractPSB          - This function looks for the top level parentheses in the
%                         given text string and returns the substring that
%                         is located between these parentheses
%
% Special functions that can be used in models
% ============================================
%   unitstepSB          - produces a unitstep 
%   heavisideSB         - same as unitstepSB
%   unitpulseSB         - produces a unitpulse for given on and off times
%   andSB               - logical and 
%   orSB                - logical or
%   piecewiseSB         - implementation of the SBML/MathML piecewise operator
%
% Other functions
% ===============
%   isparameterSB       - checks if a given "name" is a parameter in given model
%   isstateSB           - checks if a given "name" is a state in given model
%   isvariableSB        - checks if a given "name" is a variable in given model
%   isreactionSB        - checks if a given "name" is a reaction in given model
%
%   stateindexSB        - returns the number of given state in given model
%   variableindexSB     - returns the number of given variable in given model
%   reactionindexSB     - returns the number of given reaction in given model

% Information:
% ============
% Systems Biology Toolbox for MATLAB
% Copyright (C) 2005-2007 Fraunhofer-Chalmers Research Centre
% Author of the toolbox: Henning Schmidt, henning@hschmidt.de
% 
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details. 
% 
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
