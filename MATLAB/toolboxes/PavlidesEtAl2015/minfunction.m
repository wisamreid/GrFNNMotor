function [k, Features_opt, time_info, WEIGHT_OF_FREQ] = minfunction(var)

% [k, Features_opt] = minfunction(var)
% This function contains the cost function to be minimised. 
%
% INPUTS: var
% var - variables of the model, which include: 
% wgs=var(1)
% wsg=var(2)
% wgg=var(3)
% C=var(4)
% str=var(5)
% wcs=var(6)
% wsc=var(7)
% wc=var(8)
% Tc=var(9)
% taue = var(10)
% taui = var(11)
% Be = var(12)
% Bi = var(13)
% Me = var(14)
% Mi = var(15)

% OUTPUTS: [k, Features_opt, time_info]
% k            - cost function value.
% Features_opt - An optional call that returns a call array of features 
%                from the fitting procedure.
% time_info    - An optional call that contains information about time
%                delays, history, and tspan

if var(:) >= 0
    
% time_info   
totaltime = 1;
tspan = [0,totaltime];
history = [0.1, 0.1, 0.1, 0.1];

Tcc = var(9)/1000; %rescale variable

   %lag = [(Tsg,Tgs),Tgg   , Tcs      , Tsc       , Tcc]
    lag = [6*10^-3, 4*10^-3, 5.5*10^-3, 21.5*10^-3, Tcc];

WEIGHT_OF_FREQ = 20;    %This value was used in most model fittings. Only while creating file fullmodel_weakerlongloop_longdelays the value of 20 was used.
NumCond = 6;   % number of model conditions
flagSTN = 0;   % when flagSTN = 1 the adjusted STN input is included in the model.
flagC   = 0;   % when flagC = 1 the adjusted Cortical input is included in the model.
AdjSTN  = 0;   % AdjSTN not generated on first call of model_eqs so pass 0 as augument instead.
AdjC    = 0;   % AdjC not generated on first call of model_eqs so pass 0 as augument instead.

for i = 1:NumCond
   
    %store var;
    CompleteVar = var; 
    
    if i == 1;
        % Full model
    end
    
    if i == 2;
        % Model with wgs = 0
        var(1) = 0;
    end
    
    if i == 3;
        % Model with wsg = 0
        var(2) = 0;
    end
    
    if i == 4;
        % Model with wcs = 0
        var(6) = 0;
        flagSTN = 1;
    end
    
    if i == 5;
        % Model with wsc = 0
        var(7) = 0;
        flagC = 1;
    end
    
    if i == 6;
        % Model with str = 0
        var(5) = 0;
    end
    
    sol = dde23('model_eqs',lag,history,tspan,[],var,flagSTN,flagC,AdjSTN,AdjC);
     
    if i == 1;
        
        x1 = sol.x;
        y1 = sol.y;
        
        minSTN  = min(y1(1,(round(length(y1)/2)):end));
        maxSTN  = max(y1(1,(round(length(y1)/2)):end));
        meanSTN = mean(y1(1,(round(length(y1)/2)):end));
        
        minGP  = min(y1(2,(round(length(y1)/2)):end));
        maxGP  = max(y1(2,(round(length(y1)/2)):end));
        meanGP = mean(y1(2,(round(length(y1)/2)):end));
        
        meanE = mean(y1(3,(round(length(y1)/2)):end));     
   
            FMminSTN  = 5   - minSTN;
            FMmeanSTN = 65  - meanSTN;
            FMmaxSTN  = 125 - maxSTN;
            FMminGP   = 45  - minGP;
            FMmeanGP  = 100 - meanGP;
            FMmaxGP   = 155 - maxGP;
           
        freq = frequency(x1, y1, meanSTN, totaltime);
        FMfreq = 14 - freq;
        
        AdjSTN = meanE * var(6);  % AdjSTN = meanFRe*wcs
        AdjC = meanSTN * var(7);  % AdjC = meanFRs * wsc
    end
    
    if i == 2;
        % Model with wgs = 0
        x2 = sol.x;
        y2 = sol.y;
        
        zeroWGS_minSTN = min(y2(1,(round(length(y2)/2)):end));
        zeroWGS_maxSTN = max(y2(1,(round(length(y2)/2)):end));
        zeroWGS_minGP  = min(y2(2,(round(length(y2)/2)):end));
        zeroWGS_maxGP  = max(y2(2,(round(length(y2)/2)):end));
        
        %for STN
        zeroWGS_STN = zeroWGS_maxSTN - zeroWGS_minSTN;
        %for GP
        zeroWGS_GP  = zeroWGS_maxGP - zeroWGS_minGP;
        
    end
    
    if i == 3;
        % Model with wsg = 0
        x3 = sol.x;
        y3 = sol.y;
        
        zeroWSG_minSTN = min(y3(1,(round(length(y3)/2)):end));
        zeroWSG_maxSTN = max(y3(1,(round(length(y3)/2)):end));
        zeroWSG_minGP  = min(y3(2,(round(length(y3)/2)):end));
        zeroWGS_maxGP  = max(y3(2,(round(length(y3)/2)):end));
        
        %for STN
        zeroWSG_STN = zeroWSG_maxSTN - zeroWSG_minSTN;
        %for GP
        zeroWSG_GP = zeroWGS_maxGP - zeroWSG_minGP;
               
    end
    
    if i == 4;
        % Model with wcs = 0
        x4 = sol.x;
        y4 = sol.y;
        
        zeroCtxSTN_minSTN = min(y4(1,(round(length(y4)/2)):end));
        zeroCtxSTN_maxSTN = max(y4(1,(round(length(y4)/2)):end));
        zeroCtxSTN_minGP  = min(y4(2,(round(length(y4)/2)):end));
        zeroCtxSTN_maxGP  = max(y4(2,(round(length(y4)/2)):end));
        
        %for STN
        zeroCtxSTN_STN  = zeroCtxSTN_maxSTN - zeroCtxSTN_minSTN;
        %for GP
        zeroCtxSTN_GP  = zeroCtxSTN_maxGP - zeroCtxSTN_minGP;
               
    end
    
    if i == 5;
        % Model with wsc = 0
        x5 = sol.x;
        y5 = sol.y;
        
    end
    
    if i == 6;
        % Model with str = 0
        x6 = sol.x;
        y6 = sol.y;
        
        zeroSTR_minSTN = min(y6(1,(round(length(y6)/2)):end));
        zeroSTR_maxSTN = max(y6(1,(round(length(y6)/2)):end));
        zeroSTR_minGP  = min(y6(2,(round(length(y6)/2)):end));
        zeroSTR_maxGP  = max(y6(2,(round(length(y6)/2)):end));
        
        % After blocking striatum firing rates increase. We impose a penalty
        % if firing rates after blocking are lower than before blocking.
                      
            if 120 - abs(zeroSTR_maxSTN-zeroSTR_minSTN) <= 0;
              zeroSTR_STN =  0;
            else
              zeroSTR_STN = 120 - abs(zeroSTR_maxSTN-zeroSTR_minSTN);
            end
            
            if 110 - abs(zeroSTR_maxGP-zeroSTR_minGP) <= 0;
               zeroSTR_GP =  0;
            else
               zeroSTR_GP = 110 - abs(zeroSTR_maxGP-zeroSTR_minGP);
            end 
        
    end
    
    %restore var to full model
    var = CompleteVar; 
    %turns off any additional parameter
    flagSTN = 0;
    flagC = 0; 
end


%% Cost function
    
    condition1 = FMminSTN^2 + FMmeanSTN^2 + FMmaxSTN^2 +...
        FMminGP^2 + FMmeanGP^2 + FMmaxGP^2 + (WEIGHT_OF_FREQ*FMfreq)^2;
    condition2 = zeroWGS_STN^2 + zeroWGS_GP^2;
    condition3 = zeroWSG_STN^2 + zeroWSG_GP^2;  
    condition4 = zeroCtxSTN_STN^2 + zeroCtxSTN_GP^2;  
    condition5 = zeroSTR_STN^2 + zeroSTR_GP^2;
    
    k = condition1 + condition2 + condition3 + condition4 + condition5;
    
    % return info
    Features_opt = {minSTN, meanSTN, maxSTN, minGP, meanGP, maxGP, freq...
        , x1, y1, x2, y2, x3, y3, x4, y4, x5, y5, x6, y6};
    
    time_info = {lag,tspan,history};
else 
    k = 10*10^100;
    Features_opt = {};
end

end




