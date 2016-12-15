function [ ] = plotAmplitudeVectorFeild(alpha, beta1, beta2, epsilon)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 
% plotAmplitudeVectorFeild: 
%         Plots the amplitude vector feild for a given
%         oscillator
% 
% Author: Wisam Reid
%  Email: wisam@ccrma.stanford.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 
% Function Definition:
% 
%   plotAmplitudeVectorFeild(alpha, beta1, beta2, epsilon)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% 
% Arguments: 
% 
%           alpha:  (float) Hopf Oscillator Parameter
%           beta1:  (float) Hopf Oscillator Parameter
%           beta2:  (float) Hopf Oscillator Parameter
%         epsilon:  (float) Hopf Oscillator Parameter
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% 
% Returns:  N/A
%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Handle Arguments

switch nargin
    case nargin ~= 4       
        disp('Error running plotAmplitudeVectorFeild')
        disp('Please provide the right number of arguments')
        disp('Function Call: ')
        disp('plotAmplitudeVectorFeild(alpha, beta1, beta2, epsilon)')
        fprintf('\n')
        disp('type "help plotAmplitudeVectorFeild" for more details')
        fprintf('\n')
        return
end

%%

for z0 = [0:0.2:1];
f = 1;
    
% parameters for time
fs = 10000;
dur = 10; % in seconds
T = 1/fs;
time = 0:T:dur;

% ntime = length(time);

% % initialize empty arrays
% t_out = zeros(size(time));
% z_out_1 = zeros(size(time));
% 
% % initialization values
% t_out(1) = time(1);
% z_out_1(1) = z0;

dzdt = @(t,z,alpha,beta1,beta2,epsilon,f)  ...
    z*(alpha + 1i*2*pi*f + beta1*abs(z)^2 + ...
    (epsilon*beta2*abs(z)^4)/(1-epsilon*abs(z)^2)); 

% integrate the ode to obtain the amplitude
[t_out,z_out] = ode45(@(t,z) dzdt(t,z,alpha,beta1,beta2,epsilon,f),time,z0);


figure(1)
plot((abs(z_out(1:end-1))),diff(abs(z_out)));
axis square
hold on
title('Oscillator 1')
grid on
% subplot(2,2,2)
% %     plot((abs(z_out_1(1:end-1))),diff(abs(z_out_1)));
% axis square
% hold on
% plot((abs(z_out(1:end-1))),diff(abs(z_out)));
% title('Oscillator 1')
% grid on
% subplot(2,2,3:4)
% %     plot(time,abs(z_out_1))
% grid on
% hold ons
% plot(time,abs(z_out))
title('Magnitude')
% 
% figure(2)
% plot3(time,real(z_out),imag(z_out))
% hold on
% grid on
end

%% Define r_dot

% % r values
% r = 0:0.01:1;

% r_dot = @(t) alpha*r + beta1*r.^3 + ((epsilon*beta2*r.^5)/(1 - epsilon*r.^2));

% % parameters for time
% fs = 1000;
% end_t = 10; % in seconds
% T = 1/fs;
% time = 0:T:end_t;
% 
% % Initial Condition
% r0 = 0.001;
% 
% % solve r_dot
% r_dot = @(t,r,alpha,beta1,beta2,epsilon)  ...
%     [alpha*r(1) + beta1*r(1)^3 + ((epsilon*beta2*r(1)^5)/(1-epsilon*r(1)^2))];
% 
% % integrate the ode to obtain the amplitude
% [time,r] = ode45(@(t,r) r_dot(t,r,alpha,beta1,beta2,epsilon),time,[r0]);

%% Plot

regime = getOscParamRegime(alpha, beta1, beta2, epsilon); % Check for a known regime


% figure(1);
% % subplot(1,nfreqs,i)
% plot(time,r)
% title('Yo bro')
% ylim([0 1.3])
% grid on
% xlabel('Time (s)')
% ylabel('Magnitude') 

% figure(2)    
% polar(pi,1.3) 
% hold on
% polar(z(:,2),z(:,1))
% hold on
% grid on
% % pause

% figure(1)
% if regime == 1
%     plot(r,r_dot)
%     hold on 
%     plot([0 max(r)],[0 0],'k-')
%     ylim([min(r_dot),0.2])
% elseif regime == 2
%     plot(r,r_dot)
%     hold on 
%     plot([0 max(r)],[0 0],'k-')
%     ylim([min(r_dot),max(0.2,max(r_dot)+ 0.1)])
% elseif regime == 3
%     plot(r,r_dot)
%     hold on 
%     plot([0 max(r)],[0 0],'k-')
%     ylim([min(-0.1,min(r_dot)),max(0.2,max(r_dot)+ 0.1)])
% elseif regime == 4
%     plot(r,r_dot)
%     hold on 
%     plot([0 max(r)],[0 0],'k-')
%     ylim([min(-0.1,min(r_dot)),max(0.2,max(r_dot)+ 0.1)])
% else
%     display('Warning: The oscillator parameter regime could not be classified')
% end

end

