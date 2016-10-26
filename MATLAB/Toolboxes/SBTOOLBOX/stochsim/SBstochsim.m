function [varargout] = SBstochsim(varargin)
% SBstochsim: Stochastic simulation of SBmodels that only contain mass
% action kinetics. The function can either simulate the time course of the
% species trajectories, or determine the probability density function (PDF)
% for chosen species at a given time point. Three different simulation
% algorithms are available. The SBmodel needs to have a certain format,
% explained below.
%
% THIS FUNCTION WORKS CURRENTLY ONLY ON WINDOWS SYSTEMS!
%
% USAGE:
% ======
% [output] = SBstochsim(model,time)         
% [output] = SBstochsim(model,time,method)         
% [output] = SBstochsim(model,time,method,speciesPDF)         
% [output] = SBstochsim(model,time,method,speciesPDF,numberSimulations)         
% [output] = SBstochsim(model,time,method,speciesPDF,numberSimulations, other)         
%
% model: SBmodel or MAstructure.
%       There are certain limitations on an SBmodel used for stochastic
%       simulation. Please read further below.
%       An MAstructure is a MATLAB data structure defined by:
%                        MAstructure.N: stoichiometric matrix
%        MAstructure.kineticParameters: vector of kinetic parameters
%        MAstructure.initialConditions: vector of initial conditions
%                  MAstructure.species: cell-array containing species names
% time: End time for simulation OR time at which to determine the PDF for
%       the species.
% method: Numerical value from 1 to 3.
%       1: SSA
%       2: Binomial tau-leap
%       3: Poisson tau-leap
% speciesPDF: Cell-array specifying the names of the species for which to
%       determine the PDFs.
% numberSimulations: Number of simulations to be performed for the
%       computation of the PDFs.
% other: Vector with additional options.
%       other = [numberIterations, outputFrequency, coarseGrainingFactor]
%       numberIterations: Total number of iterations
%       outputFrequency: Frequency of output
%       coarseGrainingFactor (r): A tradeoff factor between CPU and accuracy. 
%           Increase in r leads to decrease in CPU and accuracy. The max r
%           should be chosen not only in terms of accuracy but also in order 
%           not to violate stability characteristics. Proposed values of r
%           are given in the references below.
%
% DEFAULT VALUES:
% ===============
% method: SSA
% numberSimulations: 1 for trajectory simulation, 100 for PDF generation
% numberIterations: 10000000
% outputFrequency: 100
% coarseGrainingFactor: 0.2
%
% Output Arguments:
% =================
% If no output arguments are given, the result of the simulation is
% plotted. Either using SBplot for the trajectory simulation, or as bar
% plot for the PDFs.
% 
% If an output argument is specified, the simulation results are returned
% as a structure.
% For trajectory simulation the structure is given as:
%                     time: vector of time instants
%                  species: cell-array of species names
%            speciesValues: matrix with simulation results
%         simulationMethod: name of the simulation method
%        numberSimulations: number of simulations (always 1)
%         numberIterations: number of iterations
%          outputFrequency: output frequency used 
%     coarseGrainingFactor: coarse graining factor used
%
% For PDF generation the structure is given as:
%                     time: time at which the PDF is generated
%                  PDFdata: structure with the fields species,
%                           populationSize, and count
%         simulationMethod: name of the simulation method
%        numberSimulations: number of simulations (always larger than 1)
%         numberIterations: number of iterations
%          outputFrequency: output frequency used 
%     coarseGrainingFactor: coarse graining factor used
%
% The stochastic simulator is implemented as a windows exectuable, developed 
% by Abhijit Chatterjee from the group of Dr Dion Vlachos at the University
% of Delaware (http://dion.che.udel.edu/). This software can perform the
% stochastic simulation algorithm (SSA), the binomial t-leap and the
% Poisson t-leap simulations for a well-mixed reacting system. For further
% details on the stochastic simulator please refer to the following
% publications:
% 
% 1. A. Chatterjee, K. Mayawala, J. S. Edwards, and D. G. Vlachos, 
% "Time accelerated Monte Carlo simulations of biological networks 
% using the binomial t-leap method", Bioinformatics, 21(9), 2136-2137 (2005). 
%
% 2. A. Chatterjee, D. G. Vlachos, and M. A. Katsoulakis, "Binomial 
% distribution basedt-leap accelerated stochastic simulation", J. Chem.
% Phys., 122, 02411-1-6 (2005).
%
% FORMAT OF THE SBmodel:
% ======================
% SBmodels that can be used for stochastic simulation need to follow some
% rules:
% 1) All reaction kinetics need to be of mass action type and be defined in
%    the following syntax:     'ReactionName' = 'kineticParameter' * ...
%                      or:     'ReactionName' = 'numeric value' * ...
% 2) All reactions have to be irreversible. Eventually you can use the 
%    function SBmakeirreversible to convert your model
% 3) The reactions can at maximum have 2 substrates and 2 products.
% 4) The right hand side of the ODEs needs only to consist of reaction rate
%    terms and eventually stoichiometric coefficients. This is necessary in
%    order to be able to determine the stoichiometric matrix.
%    More information about the required syntax can be found in the help
%    text of the function SBstoichiometry
% 5) No variables, functions, events, functionsMATLAB are allowed to be
%    present.

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
% Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307,
% USA.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHECK OPERATING SYSTEM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isunix,
    error('Sorry, the stochastic simulation is currently only available on windows systems.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHECK IF STOCHASTIC SIMULATOR IS PRESENT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if exist('stochastic_sim.exe') ~= 2,
    error(sprintf('Please downlad the stochastic simulator from\nhttp://dion.che.udel.edu/multiscale/stochastic_sim.exe and\nplace the file ''stochastic_sim.exe'' in the\nSBTOOLBOX/stochsim folder.'));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HANDLING THE VARIABLE INPUT ARGUMENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% some default values
methodChoice = 1;
numberSimulations = 1;
numberIterations = 10000000;
outputFrequency = 100;
coarseGrainingFactor = 0.2;
other = [0 0 0];
speciesGeneratePDF = {};
if nargin == 2,
    model = varargin{1};
    maxSimulationTime = varargin{2};
elseif nargin == 3,
    model = varargin{1};
    maxSimulationTime = varargin{2};
    methodChoice = varargin{3};
elseif nargin == 4,
    model = varargin{1};
    maxSimulationTime = varargin{2};
    methodChoice = varargin{3}; 
    speciesGeneratePDF = varargin{4};
elseif nargin == 5,
    model = varargin{1};
    maxSimulationTime = varargin{2};
    methodChoice = varargin{3}; 
    speciesGeneratePDF = varargin{4};
    numberSimulations = varargin{5};
elseif nargin == 6,
    model = varargin{1};
    maxSimulationTime = varargin{2};
    methodChoice = varargin{3}; 
    speciesGeneratePDF = varargin{4};
    numberSimulations = varargin{5};
    other = varargin{6};
else
    error('Incorrect number of input arguments.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ERROR CHECKING AND PROCESSING INPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(class(model),'SBmodel'),
    type = 'SBmodel';
elseif isfield(model,'N') && isfield(model,'kineticParameters') && isfield(model,'initialConditions') && isfield(model,'species'),
    type = 'MAstructure';
else
    error('Unknown model type.');
end
      
if methodChoice < 1 || methodChoice > 3,
    error('Incorrect choice of method. (1...3)');
end
if ~isempty(speciesGeneratePDF) && numberSimulations == 1,
    numberSimulations = 100;  % set default value for PDF calculation
end
if numberSimulations > 1 && isempty(speciesGeneratePDF),
    error('For PDF calculations species names have to be provided.');
end
if length(other) == 3,
    if other(1) ~= 0,
        numberIterations = other(1);
    end
    if other(2) ~= 0,
        outputFrequency = other(2);
    end
    if other(3) ~= 0,
        coarseGrainingFactor = other(3);
    end
else
    error('If specified, the ''other'' input argument needs to have three elements.');
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GETTING NEEDED DATA FOR THE SIMULATOR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if numberSimulations > 1,
    simulationType = 2;
    generatePDF = 1;
else 
    simulationType = 1;
    generatePDF = 0;
end
if methodChoice == 1,
    method = 'SSA';
elseif methodChoice == 2,
    method = 'Binomial tau-leap';
elseif methodChoice == 3,
    method = 'Poisson tau-leap';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SAVE CURRENT PATH
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
oldPath = pwd;
stochsimPath = fileparts(which('stochastic_sim.exe'));

if strcmp(type,'SBmodel'),
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CONVERT MODEL TO MA KINETICS (ONLY IF SBMODEL)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    modelMA = SBconvert2MA(model);
    N = modelMA.N;
    kineticParameters = modelMA.kineticParameters;
    initialConditions = modelMA.initialConditions;
    speciesNames = modelMA.species;
else
    N = model.N;
    kineticParameters = model.kineticParameters;
    initialConditions = model.initialConditions;
    speciesNames = model.species;   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHECK IF A MAXIMUM OF 2 SUBSTRATES AND/OR 2 PRODUCTS IN A SINGLE REACTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(find(sum(N>0)>2)) || ~isempty(find(sum(N<0)>2)),
    error('At least one reaction has more than 2 substrates and/or products.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DELETE ALL INPUT/OUTPUT FILES IN STOCHSIM FOLDER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd(stochsimPath);
delete *.txt

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CREATE THE NEEDED MODEL INPUT FILES FOR STOCHASTIC_SIM.EXE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reactants.txt
fid = fopen('reactants.txt','w');
for k = 1:length(speciesNames),
    fprintf(fid,'%i  %s\n',k,speciesNames{k});
end
fclose(fid);

% stoichiometry.txt (transposed stoichiometric matrix)
Ntranspose = N';
dlmwrite('stoichiometry.txt',Ntranspose,' ');

% input_connectivity.txt & input_stoic.txt
input_connective_data = zeros(size(Ntranspose));
input_stoic_data = zeros(size(Ntranspose));
for k = 1:size(Ntranspose,1),
    indicesNeg = find(Ntranspose(k,:)<0);
    indicesPos = find(Ntranspose(k,:)>0);
    indices = [indicesNeg indicesPos];
    input_connective_data(k,1:length(indices)) = indices;
    input_stoic_data(k,1:length(indices)) = Ntranspose(k,indices);
end
% add reaction number s in first columns
input_connective_data = [[1:size(Ntranspose,1)]' input_connective_data];
input_stoic_data = [[1:size(Ntranspose,1)]' input_stoic_data];
% write to file
dlmwrite('input_connectivity.txt',input_connective_data,' ');
dlmwrite('input_stoic.txt',input_stoic_data,' ');

% input_dyna.txt
input_dyna_data = [[1:length(kineticParameters)]' kineticParameters(:)];
dlmwrite('input_dyna.txt',input_dyna_data,' ');

% input_species.txt - generate only if not single trajectory and species
% names given for which to generate the pdf
if ~isempty(speciesGeneratePDF) && numberSimulations > 1,
    fid = fopen('input_species.txt','w');
    for k = 1:length(speciesGeneratePDF),
        index = strmatch(speciesGeneratePDF{k},speciesNames,'exact');
        if isempty(index),
            delete *.txt;
            cd(oldPath);
            error(sprintf('Species ''%s'', for which a PDF should be calculated, is not present in the model.',speciesGeneratePDF{k}));
        end
        fprintf(fid,'%i\n',index);
    end
    fclose(fid);    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CREATE THE NEEDED SIMULATOR INPUT FILE FOR STOCHASTIC_SIM.EXE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
input_data = '';
fid = fopen('input.txt','w');
fprintf(fid,'%d \t\tNumber of species\n',length(speciesNames));
fprintf(fid,'%d \t\tNumber of reactions\n',size(Ntranspose,1));
fprintf(fid,'%d \t\tSimulation algorithm: (1) SSA (2) binomial tau-leap (3) poisson tau-leap\n',methodChoice);
fprintf(fid,'%d \t\tSimulation type : (1) single trajectory (2) PDF generation\n',simulationType);
fprintf(fid,'%d \t\tNumber of trajectories\n',numberSimulations);
fprintf(fid,'%g \t\tMaximum simulation time\n',maxSimulationTime);
fprintf(fid,'%d \t\tTotal number of iterations\n',numberIterations);
fprintf(fid,'%d \t\tFrequency of output\n',outputFrequency);
fprintf(fid,'%1.1f \t\tCoarse-graining factor\n',coarseGrainingFactor);
fprintf(fid,'%d \t\tRandom number generator seed (enter a positive integer)\n',floor(rand*1000));
fprintf(fid,'%d \t\tNumber of species for generating probability distribution fn\n',length(speciesGeneratePDF));
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CREATE THE INITIAL POPULATION FILE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
input_data_population = [[1:length(initialConditions)]' round(initialConditions(:))];
dlmwrite('input_population.txt',input_data_population,' ');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% START THE SIMULATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!stochastic_sim
rehash;  % needed to be aware of the presence of the generated files

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHECK IF OUTPUT FILES HAVE BEEN GENERATED
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~generatePDF,
    % check for out_traj.txt
    fid = fopen('out_traj.txt','r');
    if fid == -1,
        delete *.txt;
        cd(oldPath);
        error('No trajectory information generated.');
    end
    lineNumberStart = 1;
    while isempty(strfind(fgetl(fid),'Time, Species')),
        lineNumberStart = lineNumberStart + 1;
    end
    fclose(fid);
    % read the CSV data from out_traj.txt
    M = csvread('out_traj.txt',lineNumberStart,0);
    time = M(:,1);
    speciesValues = M(:,2:end);
else
    % check for out_pdf.txt
    PDFdata = [];
    if length(speciesGeneratePDF) > 0,
        fid = fopen('out_pdf.txt','r');
        if fid == -1,
            delete *.txt;
            cd(oldPath);            
            error('No pdf information generated. You might need to choose a shorter simulation time.');
        end
        % parse data from PDF file for the different species
        for k = 1:length(speciesGeneratePDF),
            while isempty(strfind(fgetl(fid),'#')), end
            fgetl(fid); % get to next line
            line = fgetl(fid);
            populationSize = [];
            count = [];
            while isempty(strfind(line,'--')),
                % parse the line
                values = sscanf(line,'%d');
                populationSize = [populationSize; values(1)];
                count = [count; values(2)];
                line = fgetl(fid);
            end
            PDFdata(k).species = speciesGeneratePDF{k};
            PDFdata(k).populationSize = populationSize;
            PDFdata(k).count = count;
        end
        fclose(fid);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ERROR CHECKING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~generatePDF,
    if size(speciesValues,2) ~= length(speciesNames),
        delete *.txt;
        cd(oldPath);
        error('Simulation was not successful');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HANDLE VARIABLE OUTPUT ARGUMENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargout == 0,
    if ~generatePDF,
        % Plot trajectories
        datanames = speciesNames;
        datavalues = speciesValues;
        SBplot(time,datavalues,datanames);
    else
        % Plot PDFs if available
        pdfdata = [];
        % determine number of rows and columns for subplot window
        cols = ceil(sqrt(length(PDFdata)));
        rows = ceil(length(PDFdata)./cols);
        for k = 1:length(speciesGeneratePDF),
            subplot(rows,cols,k);
            bar(PDFdata(k).populationSize,PDFdata(k).count);
            title(sprintf('PDF for %s',PDFdata(k).species));
            ylabel('Counts');
            xlabel('Population Size');
        end
    end
elseif nargout == 1,
    % construct output structure
    if ~generatePDF,
        output.time = time;
        output.species = speciesNames;
        output.speciesValues = speciesValues;
    else
        output.time = maxSimulationTime;
        output.PDFdata = PDFdata;
    end
    output.simulationMethod = method;
    output.numberSimulations = numberSimulations;
    output.numberIterations = numberIterations;
    output.outputFrequency = outputFrequency;
    output.coarseGrainingFactor = coarseGrainingFactor;
    varargout{1} = output;
else
    delete *.txt;
    cd(oldPath);
    error('Incorrect number of output arguments!');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DELETE ALL FILES AND RETURN TO START DIRECTORY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
delete *.txt;
cd(oldPath);
return

