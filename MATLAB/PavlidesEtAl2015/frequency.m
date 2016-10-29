function freq = frequency(x, y, threshold, totaltime)

% function freq = frequency(x, y, threshold,totaltime)
% Counts period of oscillation but does not discriminate between different
% frequency bands, unlike a fast fourier transform (FFT). However, when
% simulation times are short this algorithm is more acccurate than FFT.
% When totaltime is long it performs no better than FFT. It works by
% counting the number of peaks, including partial oscillation periods at
% the beginning and end of the time series

% INPUTS 
% x         - time vector
% y         - amplitude vector
% threshold - minimum level at which oscillations are detected
% totaltime - maximum simulation time

% OUTPUTS
% freq      - frequency


[up, down] = peakdet(y(1,:), threshold, x);

if length(up) > 2 || length(down) > 2
    
    if length(up) == length(down)
        
        % count the number of peaks
        integer_freq = length(up);
        
        period = mean(diff(up(:,1)));
        
        % add any partial oscillations at the front or end of oscillation.
        front_adj = up(1,1)/period;
        end_adj = (totaltime - up(end,1))/period;
        
        % Sum and divide by totaltime for frequency.
        freq = (integer_freq + front_adj + end_adj-1)/totaltime;
        
    end
    
    if length(up) < length(down)
        
        integer_freq = length(up);
        
        period = mean(diff(up(:,1)));
        
        front_adj = up(1,1)/period;
        end_adj = (totaltime - up(end,1))/period;
        
        freq = (integer_freq + front_adj + end_adj-1)/totaltime;
        
    end
    
    if length(up) > length(down)
        
        integer_freq = length(down);
        
        period = mean(diff(down(:,1)));
        
        front_adj = down(1,1)/period;
        end_adj = (totaltime - down(end,1))/period;
        
        freq = (integer_freq + front_adj + end_adj-1)/totaltime;
        
    end
    
else
    freq = 0;
end

if isempty(up) || isempty(down)
    freq = 0;
end
end