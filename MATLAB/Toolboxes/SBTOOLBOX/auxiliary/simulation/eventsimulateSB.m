function [t,x,te,xe,ie] = eventsimulateSB(model,method,tspan,ic,options,eventFunction,eventAssignmentFunction);
% eventsimulateSB: Function realizing the simulation of systems with events
% that cause discrete state changes. Implementation as auxiliary function
% in order to be able to reuse the code when simulating systems with events
% during parameter sensitivity analysis.

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
% initialize simulation results
t = [];
x = [];
te = [];
xe = [];
ie = [];
% add the eventFunction to the integrator options
options = odeset(options,'Events',eventFunction);


% simulate the system until end time is reached
% in case the user provided a defined time instant vector when results are
% to be returned, it is necessary to delete 2 event instants from the
% vectors ...
% in case that the event time is not a measured time instant then delete
% the last element from the previous piece and the first element from the
% following piece
% in case that the event time is a member of tspan, the last two elements
% of the previous simulation data should be deleted. 
tspansim = tspan;
deleteFirstElementNext = 0;
while 1,
    [tpiece,xpiece,tepiece,xepiece,iepiece] = eval(sprintf('feval(@%s,@%s,tspansim,ic,options);',method,model));
    % save the state at the event
    pieceEndState = xpiece(end,:);
    pieceEndTime = tpiece(end);
    if ~isempty(tepiece),
        eventTime = tepiece(end);
        eventIndex = iepiece(end);
        if length(tspan) > 2,
            % ONLY DO CORRECTIONS IF tspan > 0 (USER DEFINED TIME STEPS)
            if isempty(find(tspan==eventTime)),
                % event time is not member of tspan
                % take away last entry in results
                if tpiece(end) < tspan(end),
                    % take away last entry only in the case where last time < tspan(end)!!!
                    tpiece = tpiece(deleteFirstElementNext+1:end-1);
                    xpiece = xpiece(deleteFirstElementNext+1:end-1,:);
                else
                    tpiece = tpiece(deleteFirstElementNext+1:end);
                    xpiece = xpiece(deleteFirstElementNext+1:end,:);
                end
                % take away the first entry in the
                deleteFirstElementNext = 1;
            else
                % event time is member of tspan
                % take away last two entries in results
                if tpiece(end) < tspan(end),
                    % take away last entries only in the case where last time < tspan(end)!!!
                    tpiece = tpiece(deleteFirstElementNext+1:end-2);
                    xpiece = xpiece(deleteFirstElementNext+1:end-2,:);
                else
                    tpiece = tpiece(deleteFirstElementNext+1:end);
                    xpiece = xpiece(deleteFirstElementNext+1:end,:);
                end
                deleteFirstElementNext = 0;
            end

        end
    end
    % collect simulation data
    t = [t; tpiece];
    x = [x; xpiece];
    te = [te; tepiece(:)];
    ie = [ie; iepiece(:)];
    xe = [xe; xepiece];
    % check if integration is finished - then break the loop and continue
    if pieceEndTime >= max(tspansim),
        break;
    end
    % set new initial conditions
    ic = feval(eventAssignmentFunction,eventIndex,eventTime,pieceEndState);
    % set new tspansim
    if length(tspansim) == 2,
        % if tspansim has 2 elements:
        tspansim = [eventTime tspansim(2)];
    else
        % if tspansim given as vector of time instants
        tspansimhelp = tspansim(find(tspansim > eventTime));
        tspansim = [eventTime tspansimhelp(:)'];
    end
end









return