function [N, componentNames, deltaT] = SBnetworkident(varargin)
% SBnetworkident: this function allows to identify a local network
% connectivity matrix from measurements of involved components. The
% function warns the user in case that the sampling time has been chosen
% too large to catch the fastest dynamics of the system or if moitey
% conservations are present and gives suggestions what to do in these cases.
% So far only experiment data can be handled in which constant
% perturbations have been applied to the system at a certain time point. 
% All used measurements are required to be collected using the same
% equidistant sampling time.
%
% USAGE:
% ======
% [N, componentNames, deltaT] = SBnetworkident(data1,data2,...,datan)         
% [N, componentNames, deltaT] = SBnetworkident(data1,data2,...,datan,startTimeSteps)         
% [N, componentNames, deltaT] = SBnetworkident(data1,data2,...,datan,startTimeSteps,numberTimeSteps)         
% [N, componentNames, deltaT] = SBnetworkident(data1,data2,...,datan,startTimeSteps,numberTimeSteps,excludeMeasurements)         
% [N, componentNames, deltaT] = SBnetworkident(data1,data2,...,datan,startTimeSteps,numberTimeSteps,excludeMeasurements,tol)         
%
% data1-n: SBdata objects with measurement data from perturbation experiments
% startTimeSteps: Vector containing as many elements as experiments are used
%   for the estimation. Each element describes the number of the sampling
%   instant from which on the corresponding data is to be used for
%   estimation. Typically this sampling instant should be the one after the
%   instant at which the perturbation has been applied. The first time step
%   has number 1.
% numberTimeStep: Vector containing as many elements as experiments are used
%   for the estimation. Each element describes the number of sampling
%   instants used for estimation. If set to "[]" the default values are
%   used.
% excludeMeasurements: Cell-array with names of measured components that are 
%   to be excluded in the identification. Necessary in cases of
%   moietyconservations or in cases where the used sampling time is to
%   large to capture the fastest dynamics of the system.
% tol: tolerance for the decision if equally spaced time steps in time
%   vector (used for the determination of the sampling time). Value given 
%   in percent of the mean of the time steps for the different
%   components.
%
% DEFAULT VALUES:
% ===============
% startTimeSteps: 2 (It is assumed that the measurements have started at the
%   time instant at which the perturbations have been applied.) 
% numberTimeSteps: all available time steps are used for the
%   identification (starting at startTimeSteps).
% excludeMeasurements: {} (no exclusions)
% tol: 1%
%
% Output Arguments:
% =================
% N:  network connectivity matrix
% componentNames: names of the components that have been used for the 
%   network identification. The ordering corresponds the ordering of the
%   connectivity matrix N.
% deltaT: Sampling time, underlying the estimation
%
% THEORY:
% =======
% Under the following assumptions the network connectivity matrix N
% corresponds to an estimation of the Jacobian of the system at the
% considered steady-state: 
%   - sampling time chosen small enough
%   - all components in the network are measured
%   - perturbations chosen small enough such that the response of the
%     system can be seen as linear
%
% More information about the method, implemented here, can be found in:
% H. Schmidt, K.-H. Cho, and E.W. Jacobsen (2005) "Identification of small
% scale biochemical networks based on general type system perturbations",
% FEBS Journal, 272, 2141-2151

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

% tolerance factor for the determination of large sampling times / moiety
% conservations
eigValueTolerance = 1;
warning off;
% info
disp('Network identification');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROCESS VARIABLE INPUT ARGUMENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin == 0,
    error('Wrong number of input arguments.');
end
% determine the number of SBdata input arguments (need to come first)
numberExperiments = 0;
for k = 1:nargin,
    if strcmp('SBdata',class(varargin{k})),
        numberExperiments = numberExperiments + 1;
    else
        break;
    end
end
if numberExperiments == 0,
    error('No experiment data given.');
end
if nargin-numberExperiments == 0,
    startTimeSteps = [];    
    numberTimeSteps = [];
    excludeMeasurements = {};
    tol = 1;
elseif nargin-numberExperiments == 1,
    startTimeSteps = varargin{numberExperiments+1};
    numberTimeSteps = [];
    excludeMeasurements = {};
    tol = 1;
elseif nargin-numberExperiments == 2,
    startTimeSteps = varargin{numberExperiments+1};
    numberTimeSteps = varargin{numberExperiments+2};
    excludeMeasurements = {};
    tol = 1;
elseif nargin-numberExperiments == 3,
    startTimeSteps = varargin{numberExperiments+1};
    numberTimeSteps = varargin{numberExperiments+2};
    excludeMeasurements = varargin{numberExperiments+3};
    tol = 1;
elseif nargin-numberExperiments == 3,
    startTimeSteps = varargin{numberExperiments+1};
    numberTimeSteps = varargin{numberExperiments+2};
    excludeMeasurements = varargin{numberExperiments+3};
    tol = varargin{numberExperiments+4};
end
% check length of pertStartTimeSteps
if length(startTimeSteps) ~= numberExperiments,
    error('Number of elements in ''startTimeSteps'' does not match number of experiments.');
end
% set startTimeSteps
if isempty(startTimeSteps),
    startTimeSteps = 2*ones(1,numberExperiments);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHECK EXPERIMENTAL DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Need to have same sampling time
% Need to have components specified in the same order
% Also determine deltaT and the measurement matrices
% Adjust data to startTimeSteps
% Check that no NaN values present in measurement matrices
measurements = {};
[dummy1,dummy2,measurements{1},dummy3,dummy4] = SBgetdata(varargin{1});
[compareDeltaT,compareNames, dummy1, dummy2, dummy3] = SBgetsamplingtimedata(varargin{1},tol);
allDeltaT = compareDeltaT;
for k = 1:numberExperiments,
    [dummy1,dummy2,measurementsk,dummy3,dummy4] = SBgetdata(varargin{k});
    % adjust measurementsk using the data in startTimeSteps
    measurementsk = measurementsk(startTimeSteps(k):end,:);
    if ~isempty(numberTimeSteps),
        % check that this number of time steps is available
        if size(measurementsk,1) >= numberTimeSteps,
            measurementsk = measurementsk(1:numberTimeSteps(k),:);
        else
            error('The chosen number of time steps for the %d-th experiment data\nis larger than the measured time steps. Max: %d',k,size(measurementsk,1));
        end
    end
    measurements{k} = measurementsk;
    % check if NaN present in measurements
    if ~isempty(find(isnan(measurements{k}))),
        error('NaN elements present in data of the %d-th experiment data.',k);
    end
    [deltaTk,componentNamesk, dummy1, dummy2, dummy3] = SBgetsamplingtimedata(varargin{k},tol);
    allDeltaT = [allDeltaT, deltaTk];
    if length(compareNames) ~= length(componentNamesk),
        error('Number of measured components in first experiment does not match the number in %d-th experiment.',k);
    end
    for k2 = 1:length(compareNames),
        if ~strcmp(compareNames{k2},componentNamesk{k2}),
            error('Measured components in first experiment do not match components in %d-th experiment.\nThe order of the components needs to be exactly the same.',k);
        end
    end
end
% check that only one single sampling time is used for all data
% check first if zero element present in allDeltaT
if ~isempty(find(allDeltaT == 0)),
    error('Measurement data at least partially contains no equidistant sampling.');
end
% no non-equidistant sampling present 
% now check if different sampling times in different measurements
deltaT = checksamplingtimevectorSB(allDeltaT,tol);
if deltaT == 0,
    error('The different data objects contain measurements collected with different sampling times.');
end
% deltaT is now determined and will be used in the following
% measurements is a cell-array containing all measurement data in matrices

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HANDLE MEASURED COMPONENTS THAT ARE NOT TO BE 
% TAKEN INTO ACCOUNT FOR THE IDENTIFICATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
excludeMeasurementsIndices = [];
for k1 = 1:length(excludeMeasurements),
    foundName = 0;
    for k2 = 1:length(compareNames),
        if strcmp(excludeMeasurements{k1},compareNames{k2}),
            excludeMeasurementsIndices(k1) = k2;
            foundName = 1;
            break;
        end
    end
    if foundName == 0,
        error('Excluded measurement ''%s'' is not available as measured component.',excludeMeasurements{k1});
    end
end
useComponentIndices = setdiff(1:length(compareNames),excludeMeasurementsIndices);
for k = 1:length(measurements),
    measurementsk = measurements{k};
    measurements{k} = measurementsk(:,useComponentIndices);
end
% take away the names from the component cell arry if not used
componentNames = compareNames(useComponentIndices);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DETERMINE THINGS THAT ARE NEEDED FOR THE ALGORITHM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% number of measured components
n = length(componentNames);
% number of measurements for different experiments
mi = [];
for k = 1:length(measurements),
    mi(k) = size(measurements{k},1);
end
% number of experiments
r = length(measurements);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONSTRUCT THE Ri AND Mi MATRICES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ri = [xi_mi, xi_mi-1, ..., xi_2]
% Mi = [xi_mi-1, xi_mi-2, ..., xi_1]
% xi_j: measurement vector in i-th experiment at j-th timestep
Ri = {};
Mi = {};
for k = 1:r,
    measurementsk = measurements{k};
    [nrtimesteps,nrcomponents] = size(measurementsk);
    X = rot90(measurementsk,3);
    Ri{k} = X(:,1:nrtimesteps-1);
    Mi{k} = X(:,2:nrtimesteps);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONSTRUCT THE R MATRIX
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% R = [R1,R2, ..., Rr]
R = [];
for k = 1:r,
    R = [R Ri{k}];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONSTRUCT THE M MATRIX (assuming single constant perturbation present
% already in first time step for each experiment)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Miall = []; 
% construct measurement and result matrices
for k = 1:r,
    Miall = [Miall Mi{k}];
end
M = [Miall; zeros(r,size(Miall,2))];
% fill in the ones where they belong
indexStart = 1;
indexEnd = 0;
for k = 1:r,
    indexEnd = indexEnd + size(Mi{k},2);
    M(size(Miall,1)+k,indexStart:indexEnd) = ones(1,indexEnd-indexStart+1);
    indexStart = indexEnd+1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHECK THAT THE MINIMUM NUMBER OF MEASUREMENTS IS AVAILABLE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The minimum number of measurements required to perform this 
% estimation is given by the fact that the M matrix needs to have 
% at least as many columns as rows.
if size(M,1) > size(M,2),
    error('To few measured time steps to perform the identification.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DO THE CALCULATION OF Ad and A
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
AdUdIdent = R*pinv(M);
AdIdent = AdUdIdent(:,1:n);
UdIdent = AdUdIdent(:,n+1:size(AdUdIdent,2));
AIdentZOH = logm(AdIdent)/deltaT;
AIdentEuler = (AdIdent-eye(n))/deltaT;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHECK PROBLEM WITH TO SLOW SAMPLING TIME / MOIETY CONSERVATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If the magnitude of an eigenvalue of the estimated discrete jacobian is
% significantly smaller than exp(-1), this indicates that a problem
% occurred and that the measurement matrix is probably very close to
% singular.
numberSingularities = length(find(abs(eig(AdIdent)) < eigValueTolerance*exp(-1)));
if numberSingularities > 0,
    disp('The sampling time might be too large to capture the fastest dynamics');
    disp('of the system. Or there might be moiety conservations present.');
    disp('This can lead to an arbitrarily wrong identification result.');
    disp('It is however possible to estimate a reduced network connectivity');
    disp(sprintf('matrix. To do this you need to exclude %d measured components',numberSingularities));
    disp('from the identification.');
    % determine the singular output vectors
    [U,S,V] = svd(M);
    Us = U(:,size(U,2)-numberSingularities+1:end)';
    % keep only the part belonging to the measured components
    Us = Us(:,1:n);
    for k = 1:size(Us,1),
        Us(k,:) = abs(Us(k,:)/max(abs(Us(k,:))));
    end
    Us(find(abs(Us)<1e-10)) = 0;
    maxElements = [];
    for k = 1:size(Us,1),
        [maxValue,maxValueIndex] = max(Us(k,:));
        maxElements = [maxElements, maxValueIndex];
    end
    maxElementsUnique = unique(maxElements);
    disp('In order to determine a reduced network connectivity matrix, you should');
    disp('run the algorithm again and specify the following measured components as');
    disp('excluded from the identification.');
    disp(' ');
    text = '';
    for k = 1:length(maxElementsUnique),
        text = sprintf('%s%s\n',text,componentNames{maxElementsUnique(k)});
    end
    disp(text);
    % might be that the same element is the largest in more than one
    % vector. Then the procedure is to be performed again.
    if length(maxElementsUnique) ~= length(maxElements),
        disp('After exclusion of the above components the algorithm might tell you to');
        disp('exclude additional components. These then have to be added to the other');
        disp('excluded ones.');
        disp(' ');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHECK WHICH ESTIMATION TO RETURN (ZOH OR EULER)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Euler estimation is returned in case that the ZOH estimation leads to a 
% complex result
if max(max(abs(imag(AIdentZOH)))) ~= 0,
    N = AIdentEuler;
    disp('Due to a complex estimation result using the zero-order-old discretization');
    disp('the Euler discretization has been used in this case. Complex results are due');
    disp('to at least one negative real eigenvalue of the discrete network connectivity');
    disp('matrix, caused by a too long sampling time or excessive measurement noise.');
else
    N = AIdentZOH;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PREPARE THE RETURN VALUES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% N
% deltaT
% componentNames
return