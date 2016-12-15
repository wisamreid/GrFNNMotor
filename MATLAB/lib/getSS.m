function [r_star, psi_star] = getSS(f_osc, f_input, alpha, beta1, beta2, epsilon, F, display_flag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 
%  getSS: Calculate steady state amplitude and phase 
%         for a given oscillator under periodic forcing
% 
% Author: Wisam Reid
%  Email: wisam@ccrma.stanford.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 
% Function Definition:
% 
%   getSS(f_osc, f_input, alpha, beta1, beta2, epsilon, F, display_flag)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% 
% Arguments: 
% 
%           f_osc:  (float) Oscillator Frequency
%         f_input:  (float) Forcing Input Frequency
%           alpha:  (float) Hopf Oscillator Parameter
%           beta1:  (float) Hopf Oscillator Parameter
%           beta2:  (float) Hopf Oscillator Parameter
%         epsilon:  (float) Hopf Oscillator Parameter
%               F:  (float) Forcing Amplitude [F = 0 is autonomous] 
%    display_flag:  [optional] (boolean) print argument details to console
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% 
% Returns:  [r_star, psi_star]
% 
%          r_star: A float, steady state gain
%        psi_star: A float, steady state phase
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%% Handle Arguments

switch nargin
    case nargin == 7
        display_flag = 1;
    case nargin ~= 8 ||  nargin ~= 7       
        disp('Error running getSS')
        disp('Please provide the right number of arguments')
        disp('Function Call: ')
        disp('getSS(f_osc, f_input, alpha, beta1, beta2, epsilon, F, print_flag)')
        fprintf('\n')
        disp('type "help getSS" for more details')
        fprintf('\n')
        return
end 

% Check input arguments
if F < 0
  error('F > 0')
end


if beta2 == 0 && epsilon ~= 0
    % set epsilon = 0 when b2 = 0 to avoid problem in root finding
    epsilon = 0; 
end

if epsilon == 0 && beta2 ~= 0
    error('epsilon cannot be zero if beta2 is not zero')
end

%% Displays

if display_flag
    % Display
    disp(['Running getFP:'])
    fprintf('\n')
    disp(['The Oscillator Frequency is: ',num2str(f_osc), ' Hz'])
    disp(['The Input Frequency is: ',num2str(f_input), ' Hz'])
    fprintf('\n')
    disp(['alpha is: ',num2str(alpha)])
    disp(['beta1 is: ',num2str(beta1)])
    disp(['beta2 is: ',num2str(beta2)])
    disp(['epsilon is: ',num2str(epsilon)])
    disp(['F is: ',num2str(F)])
    fprintf('\n')
end

%% Utilities

% boolean for solution warning messages
print_warnings = 0; 

%% Find ss amplitudes numerically for forced or autonomous oscillators

if F % oscillator is driven by a complex sinusoid
    
        % calculate Omega
    Omega = 2*pi*(f_osc - f_input);

    r = sqrt(roots([ ...
              (beta1-beta2)^2*epsilon^2,...
              -2*(beta1-beta2)*(beta1-alpha*epsilon)*epsilon,...
              beta1^2-4*alpha*beta1*epsilon+(epsilon*alpha^2+2*alpha*beta2+epsilon*Omega^2)*epsilon,...
              -2*alpha^2*epsilon+2*alpha*beta1-(epsilon*F^2+2*Omega^2)*epsilon,...
              alpha^2+2*epsilon*F^2+Omega^2,...
              -F^2]));

    % take only real roots
    r = r(find(abs(imag(r)) < eps('single'))); 
    r = real(r);

    % remove multiple roots
    r = sort(unique(r),'descend'); 

else % autonomous oscillator
    
    % this is r_dot = 0 with some algebraic manipulation
    r = roots([epsilon*(beta2-beta1), 0, beta1-epsilon*alpha, 0, alpha, 0]);

    % only unique real values
    r = real(unique(r(find(abs(imag(r)) < eps('single')))));

    r = r(find(r >=0)); % no negative amplitude

end

% take r's below the asymptote
if beta2
  r = r(find(r < 1/sqrt(epsilon))); 
end

%% Find corresponding relative phrases

% For now we are only going to find the steady-state phase for driven 
% oscillators
if F 
    % get the sign of psi
    signPsi = (Omega >= 0)*2 - 1;

    psi = signPsi*acos(-(alpha*r + beta1*r.^3 + beta2*epsilon*r.^5./(1-epsilon*r.^2))/F);

    if print_warnings

      maxImag = max(abs(imag(psi)));

      if  maxImag > eps('single')*100
        disp(['Warning: significant nonzero imaginary part in psi ('...
          num2str(maxImag) ') for Omega = ' num2str(Omega)])
      end

    end

    psi = real(psi);
else
    psi = [];
end

%% output

r_star = r;
psi_star = psi;

% % ========================================================
% function drdotdr = slope(r, a, b1, b2, e)
% drdotdr = a + 3*b1*r.^2 + (5*e*b2*r.^4-3*e^2*b2*r.^6)./((1-e*r.^2).^2);


% %% Get Jacobian matrix
% J11 = alpha + 3*beta1*r.^2 + epsilon*beta2*r.^4.*(5-3*epsilon*r.^2)./(1-epsilon*r.^2).^2;
% J12 = -F*sin(psi);
% J21 = F*sin(psi)./r.^2;
% J22 = -F*cos(psi)./r;
% delta = J11.*J22 - J12.*J21; % determinant of Jacobian
% tau = J11 + J22; % trace of Jacobian
% chDet = tau.^2 - 4*delta; % determinant of characteristic eq
% 
% %% Calculate stability type
% % 4 = stable node, 3 = stable spiral, 2 = unstable node,
% % 1 = unstable spiral, 0 = saddle point
% stabType = 2*sign(delta) -1*sign(tau) + .5*sign(chDet) + .5;
% stabType(find(stabType < 0)) = 0;
% stability = (stabType >= 3); % 1 = stable, 0 = unstable
% 
% %% Prepare output
% if All % both stable and unstable fixed points
%   rStar = r;
%   psiStar = psi;
% else % only stable fixed points
%   indStab = find(stability);
%   rStar = r(indStab);
%   psiStar = psi(indStab);
%   stability = stability(indStab);
%   stabType = stabType(indStab);
% end
% 
% % ========================================================
% function drdotdr = slope(r, a, b1, b2, e)
% drdotdr = a + 3*b1*r.^2 + (5*e*b2*r.^4-3*e^2*b2*r.^6)./((1-e*r.^2).^2);

% ========================================================
% ========================================================
% ========================================================
% ========================================================

%% Determine parameter regime

% % First, we need to know which bifurcation parameter regime we 
% % are operating in
% regime = getOscParamRegime(alpha, beta1, beta2, epsilon);

%%
% Next, we can find the boundary between stable nodes and stable spirals by solving:
%     
%     discriminant = 0
%     r_dot = 0 
%     psi_dot = 0 
%     
% simultaneously where T and Delta  are the trace and determinant of the 
% Jacobian matrix evaluated at a fixed point (see Section 2). 
% We find that the boundary is at

% %% Prepare output
% 
% 
% indStab = find(stab
% ility);
% r_star = r(indStab);
% psi_star = psi(indStab);
% 
% 
% % stability = stability(indStab);
% % stabType = stabType(indStab);
% 
% % % Return all SS Solutions for 
% % % both stable and unstable fixed points
% % r_star = r;
% % psi_star = psi;


%% Now we handle each quadrant correctly

% % % % % first quadrant
% % % % % do nothing
% % % % 
% % % % % second quadrant
% % % % if psi_star > pi/2 && psi_star < pi
% % % %     % do nothing
% % % % end
% % % % 
% % % % % third quadrant 
% % % % if psi_star > pi && psi_star < 3*pi/2
% % % %     psi_star = psi_star - pi;
% % % %     r_star = -1*r_star;
% % % % end
% % % % 
% % % % % fourth quadrant 
% % % % if psi_star > 3*pi/2
% % % %     psi_star = mod(psi_star + pi, 2*pi);
% % % %     r_star = -1*r_star;
% % % % end
% % % % 
% % % % % edge cases
% % % % % if psi_star == 0 || psi_star == pi
% % % % %    r_star = -1*r_star; 
% % % % % end
% % % % 
% % % % % % print
% % % % % r_star
% % % % % psi_star

