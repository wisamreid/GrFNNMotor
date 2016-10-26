function [X, FVAL, EXITFLAG] = genalgSB(FUN,X0,varargin)
% genalgSB: Evolutionary optimization algorithm
% 
% USAGE:
% ======
% [X,FVAL,EXITFLAG] = genalgSBPD(FUN,X)
% [X,FVAL,EXITFLAG] = genalgSBPD(FUN,X,OPTIONS)
%
% FUN: Function to optimize
% X: Starting Guess
% OPTIONS: structure containing options for the algorithm:
%       OPTIONS.maxgenerations: Number of generations to be evolved
%       OPTIONS.maxtime: Maximum time (in minutes) for optimization
%       OPTIONS.initpop_sigma: Parameter/initial condition std. dev. for initial population
%       OPTIONS.pop_size: Population size (EVEN numer)
%       OPTIONS.p_cross: Crossover probability
%       OPTIONS.p_mut: Mutation probability
%       OPTIONS.mut_sigma: Std. dev. for mutation
%       OPTIONS.highbounds: (NOT USED YET) vector containing upper bounds for parameters
%           Instead of a vector highbounds can also be a scalar > 1. In the
%           latter case the highbounds are determined by multiplying this
%           scalar to the initial parameter guess X. This latter case works
%           only for positive X(i) ... if negative present then the user
%           needs to specify a whole lowbounds vector. 
%       OPTIONS.lowbounds: (NOT USED YET) vector containing lower bounds for parameters.
%           Instead of a vector lowbounds can also be a scalar < 1. In the
%           latter case the lowbounds are determined by multiplying this
%           scalar to the initial parameter guess X. This latter case works
%           only for positive X(i) ... if negative present then the user
%           needs to specify a whole highbounds vector. 
%       OPTIONS.outputFunction: (NOT USED YET) string with output function name. If
%               not given or if empty then no output function will be used.
%               This output function can be used to display data, or to 
%               save optimization information to a file. The function needs
%               to be in the MATLAB path. Its calling syntax is:
%                   'outputFunction'(bestparameters,bestfunctionvalue,currentsimplex)
%       OPTIONS.silent: =0: output of info, =1: no output
%
% DEFAULT VALUES:
% ===============
% OPTIONS.maxgenerations:   500
% OPTIONS.maxtime:          120
% OPTIONS.initpop_sigma:    0.5
% OPTIONS.pop_size:         100
% OPTIONS.p_cross:          0.8
% OPTIONS.p_mut:            0.3
% OPTIONS.mut_sigma:        1.5
% OPTIONS.lowbounds:        0.001  => lowbounds = 0.001*X 
% OPTIONS.highbounds:       1000  => highbounds = 1000*X 
% OPTIONS.outputFunction:   no output function ('')
% OPTIONS.silent:           0 (no output of info)
%
% Output Arguments:
% =================
% X: Found solution
% FVAL: Value of the function FUN at X
% EXITFLAG: 1=success, 0=not found

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

global nfunk
nfunk = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin == 2,
    OPTIONS = [];
elseif nargin == 3,
    OPTIONS = varargin{1};
else
    error('Incorrect number input of arguments.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set default values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EXITFLAG = 1;
maxgenerations = 500;
maxtime = 120;
initpop_sigma = 0.5;
pop_size = 100;
p_cross = 0.8;
p_mut =  0.3;
mut_sigma = 1.5;
lowbounds = 0.001; 
highbounds = 1000;  
outputFunction = '';
silent = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Handle OPTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% silent
if isfield(OPTIONS,'silent'),
    if ~isempty(OPTIONS.silent),
        silent = OPTIONS.silent;
    end
end
% maxgenerations 
if isfield(OPTIONS,'maxgenerations '),
    if ~isempty(OPTIONS.maxgenerations ),
        maxgenerations  = OPTIONS.maxgenerations ;
    end
end
% initpop_sigma 
if isfield(OPTIONS,'initpop_sigma'),
    if ~isempty(OPTIONS.initpop_sigma),
        initpop_sigma = OPTIONS.initpop_sigma;
    end
end
% maxtime
if isfield(OPTIONS,'maxtime'),
    if ~isempty(OPTIONS.maxtime),
        maxtime = OPTIONS.maxtime;
    end
end
% pop_size
if isfield(OPTIONS,'pop_size'),
    if ~isempty(OPTIONS.pop_size),
        pop_size = OPTIONS.pop_size;
    end
end
% p_cross
if isfield(OPTIONS,'p_cross'),
    if ~isempty(OPTIONS.p_cross),
        p_cross = OPTIONS.p_cross;
    end
end
% p_mut
if isfield(OPTIONS,'p_mut'),
    if ~isempty(OPTIONS.p_mut),
        p_mut = OPTIONS.p_mut;
    end
end
% mut_sigma
if isfield(OPTIONS,'mut_sigma'),
    if ~isempty(OPTIONS.mut_sigma),
        mut_sigma = OPTIONS.mut_sigma;
    end
end
% lowbounds
if isfield(OPTIONS,'lowbounds'),
    if ~isempty(OPTIONS.lowbounds),
        if length(OPTIONS.lowbounds) == ndim,
            lowbounds = OPTIONS.lowbounds;
        else
            if length(OPTIONS.lowbounds) == 1,
                if ~isempty(find(X0<0)),
                    error('Please specify a vector for the lower bounds using OPTIONS.lowbounds.');
                end
                lowbounds = X0*OPTIONS.lowbounds;
            else
                error('The OPTIONS.lowbounds setting is not correct.');
            end
        end
    end
else
    if ~isempty(find(X0<0)),
        error('Please specify a vector for the lower bounds using OPTIONS.lowbounds.');
    end
end
% highbounds
if isfield(OPTIONS,'highbounds'),
    if ~isempty(OPTIONS.highbounds),
        if length(OPTIONS.highbounds) == ndim,
            highbounds = OPTIONS.highbounds;
        else
            if length(OPTIONS.highbounds) == 1,
               if ~isempty(find(X0<0)),
                    error('Please specify a vector for the higher bounds using OPTIONS.highbounds.');
               end
                highbounds = X0*OPTIONS.highbounds;
            else
                error('The OPTIONS.highbounds setting is not correct.');
            end
        end
    end
else
    if ~isempty(find(X0<0)),
        error('Please specify a vector for the higher bounds using OPTIONS.highbounds.');
    end
end
lowbounds = lowbounds(:)';
highbounds = highbounds(:)';
% outputFunction
outputFunction = '';
if isfield(OPTIONS,'outputFunction'),
    if ~isempty(OPTIONS.outputFunction),
        outputFunction = OPTIONS.outputFunction;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHECK BOUND VIOLATION FOR INITIAL GUESS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
indexXhi = find(X0 > highbounds);
if ~isempty(indexXhi),
    error('Initial guess does violate high parameter bounds.');
end
indexXlo = find(X0 < lowbounds);
if ~isempty(indexXlo),
    error('Initial guess does violate low parameter bounds.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIAL POINT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cost = feval(FUN,X0);
nfunk = nfunk + 1;
if ~silent,
    disp(' Nr Gen    Nr Fun Eval             Fitness');
    disp(sprintf('  %5.0f          %5.0f        %12.6g            %s', 0, nfunk, cost));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize population and plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n_param = length(X0);
init_individual(1:n_param) = X0;
population = (ones(pop_size,1) * init_individual) .* (10.^(initpop_sigma .* randn(pop_size, n_param)));
fitness = zeros(1,pop_size);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% evolve population
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic; % start timer
for gen = 1:maxgenerations
    if toc/60 > maxtime,
        EXITFLAG = 0;
        if ~silent,
            disp('Exceeded maximum time.');
        end
        break;
    end    
    % rate individuals
    for i = 1:pop_size
        individual = population(i,:);        
        cost =  eval(sprintf([FUN,'([',num2str(individual),'])']));
        fitness(i) = 1/cost;
        nfunk = nfunk + 1;
    end
    % create new generation
    temp_population = population;
    [best_value(gen) best_index(gen)] = max(fitness);
    temp_population(1,:) = population(best_index(gen),:); % elitism
    temp_population(2,:) = population(best_index(gen),:); % elitism
    for i = 3:2:pop_size
        individual1 = selectIndividual(fitness,population);
        individual2 = selectIndividual(fitness,population);
        if p_cross > rand
            [cross_individual1, cross_individual2] = crossover(individual1,individual2);
            temp_population(i,:) = cross_individual1;
            temp_population(i+1,:) = cross_individual2;
        else
            temp_population(i,:) = individual1;
            temp_population(i+1,:) = individual2;
        end
    end
    % mutate
    for i = 3:pop_size
        temp_population(i,:) = mutate(temp_population(i,:),p_mut,mut_sigma);
    end
    population = temp_population;
    % disp best value of generation
    if ~silent,
        disp(sprintf('  %5.0f          %5.0f        %12.6g            %s', 0, nfunk, cost));
    end
end
% get best ind. and return it
% rate individuals
for i = 1:pop_size
    individual = population(i,:);
    cost = cost2(individual);
    fitness(i) = 1/cost;
end
[best_value(gen) best_index(gen)] = max(fitness);
X = population(best_index(gen),:);
FVAL = eval(sprintf([FUN,'([',num2str(X),'])']));                                         
nfunk = nfunk + 1;
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% select individual
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function individual = selectIndividual(fitness,population)
r = rand*sum(fitness);
inc = fitness(1);
i = 1;
while 0 == 0
    if inc > r
        individual = population(i,:);
        break
    end
    i = i+1;
    inc = inc+fitness(i);
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% crossover
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [new_individual1, new_individual2] = crossover(individual1,individual2)
l = length(individual1);
random_order = randperm(l);
break_point = floor(rand*(l-1)) + 2;
new_individual1(random_order(1:break_point-1)) = individual1(random_order(1:break_point-1));
new_individual1(random_order(break_point:l)) = individual2(random_order(break_point:l));
new_individual2(random_order(1:break_point-1)) = individual2(random_order(1:break_point-1));
new_individual2(random_order(break_point:l)) = individual1(random_order(break_point:l));
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mutate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mutated_individual = mutate(individual,p_mut,mut_sigma)
mutated_individual = individual .* 10.^( (p_mut > rand(1,length(individual))) .* (mut_sigma .* randn(1,length(individual))) );
return