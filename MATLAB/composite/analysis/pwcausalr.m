function [pp,cohe,Fx2y,Fy2x,Fxy]=pwcausalr(x,Nr,Nl,porder,fs,freq)   
% Using Geweke's method to compute the causality between any two channels 
% 
%   x is a two dimentional matrix whose each row is one variable's time series 
%   Nr is the number of realizations, 
%   Nl is the length of every realization 
%      If the time series have one ralization and are stationary long, just let Nr=1, Nl=length(x) 
%   porder is the order of AR model 
%   fs is sampling frequency 
%   freq is a vector of frequencies of interest, usually freq=0:fs/2 % CK: WRONG!! freq must be a scalar, else the for loop doesn't work.
% 
%   Fx2y is the causality measure from x to y 
%   Fy2x is causality from y to x 
%   Fxy is instantaneous causality between x and y 
%        the order of Fx2y/Fy2x is 1 to 2:L, 2 to 3:L,....,L-1 to L.  That is, 
%        1st column: 1&2; 2nd: 1&3; ...; (L-1)th: 1&L; ...; (L(L-1))th: (L-1)&L. 

% revised Jan. 2006 by Yonghong Chen 
% Note: remove the ensemble mean before using this code 

[L,N] = size(x); %L is the number of channels, N is the total points in every channel 

%phase=0;
index = 0; 
for i = 1:L-1 
    for j = i+1:L 
        index = index + 1; 
        y(1,:) = x(i,:); 
        y(2,:) = x(j,:);   
        [A2,Z2] = armorf(y,Nr,Nl,porder); %fitting a model on every possible pair 
        eyx = Z2(2,2) - Z2(1,2)^2/Z2(1,1); %corrected covariance 
        exy = Z2(1,1) - Z2(2,1)^2/Z2(2,2); 
        f_ind = 0; 
        for f = 0:1:freq 
            f_ind = f_ind + 1; 
            [S2,H2] = spectrum_AR(A2,Z2,porder,f,fs); 
            pp(i,f_ind) = abs(S2(1,1)*2);      % revised 
            if (i==L-1)&(j==L) 
                pp(j,f_ind) = abs(S2(2,2)*2);  % revised 
            end 
            cohe(index,f_ind) = real( abs(S2(1,2))^2 / S2(1,1)/S2(2,2) );   
            Fy2x(index,f_ind) = log(abs(S2(1,1))/abs(S2(1,1)-(H2(1,2)*eyx*conj(H2(1,2)))/fs)); %Geweke's original measure 
            Fx2y(index,f_ind) = log(abs(S2(2,2))/abs(S2(2,2)-(H2(2,1)*exy*conj(H2(2,1)))/fs)); 
            Fxy(index,f_ind) = log(abs(S2(1,1)-(H2(1,2)*eyx*conj(H2(1,2)))/fs)*abs(S2(2,2)-(H2(2,1)*exy*conj(H2(2,1)))/fs)/abs(det(S2))); 
            
        end

    end 
end 
