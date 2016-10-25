%% SPECTRUM
% 
% This code, based on spectrum.py, calculates various
% spectrum-related quantities.
%
% This function actually contains several different functions, which are
% selected by the first argument, which must be a string. The choices are:
%     'tf': Calculate the frequency vector from the time vector
%     'pl': Calculate the "plain" FFT with optional smoothing
%     'mt': Calculate the spectrum via the multitaper method (pmtm)
%     'ar': Calculate the spectrum via BSMART's autoregressive method
%     'gr': Calculate the spectrogram
%     'gp': Plot the spectrogram
%
% Each of these options is described in detail below.
%
%
% TF
%
% This function converts a time vector into a frequency vector.
%
% Usage:
%   F=spectrum('tf',x,[fs])
% where
%   F is the vector of frequencies
%   x is the time vector
%   fs is the optional sampling rate (default calculated from x)
%
%
% PL
%
% This script returns the plain, smoothed FFT power spectrum for an input vector x.
%
% Usage:
%     P,F = spectrum('pl',x,[rate,smooth,padding])
% where:
%     P = power spectrum
%     F = corresponding frequencies
%     x = input time series
%     rate = data sampling rate (default 200 Hz)
%     smooth = amount of smoothing (default 100)
%     pad = amount of padding to avoid edge effects (default 0)
%
%
% MT
% 
% This script returns the multitaper power spectrum for an input vector x.
% 
% Usage:
%     P,F = spectrum('mt',x,[rate,bandwidth])
% where:
%     P = power spectrum
%     F = corresponding frequencies
%     x = input time series
%     rate = data sampling rate (default 200 Hz)
%     bandwidth = time*bandwidth product of the multitapers (default 4)
% 
%
% AR
% 
% This script returns the autoregressive power spectrum for an input vector
% x. It requires BSMART (specifically, armorf and spectrum_AR) to function.
% 
% Usage:
%     P,F = ar(x,[rate,order,maxfreq,ntrials])
% where:
%     P = power spectrum
%     F = corresponding frequencies
%     x = input time series
%     rate = data sampling rate (default 200 Hz)
%     order = polynomial order for autoregressive fitting (default 12)
%     maxfreq = maximum frequency to return (default 0 = calculate automatically)
%     ntrials = number of trials to separate the signal into (default 1)
%
%
% GR
% 
% This script plots a spectrogram of a given input vector or file.
% 
% Usage:
%     [T,F,arr]=spectrum('gr',x,[rate,window,maxfreq,tsmooth,fsmooth])
% where:
%     T       is a 1xM vector of times (y-axis)
%     F       is a 1xN vector of frequencies (x-axis)
%     arr     is an NxM array of spectral power at each time, frequency point
%     x       is the input time series
%     rate    is the sampling rate in Hz (default 200)
%     window  is the window length in s (default 2)
%     maxfreq is the maximum frequency to plot in Hz (default 50)
%     tsmooth is the amount of smoothing in the time domain (default 2)
%     fsmooth is the amount of smoothing in the frequency domain (default 2)
%
%
% GP
%
% Same input arguments as GR, except plots the spectrogram instead of
% returning it.
%
% Usage:
%     spectrum('gp',x,[rate,window,maxfreq,tsmooth,fsmooth])
%
% Version: 2011dec09 by Cliff Kerr (cliffk@neurosim.downstate.edu)

function varargout=spectrum(varargin)

noutputargs=nargout; % Store the number of output arguments

whichtype=varargin{1};
switch whichtype
    case 'tf'
        varargout{1}=timefreq(varargin{2:end});
    case 'pl'
        [varargout{1},varargout{2}]=pl(varargin{2:end});
    case 'ar'
        [varargout{1},varargout{2}]=ar(varargin{2:end});
    case 'mt'
        [varargout{1},varargout{2}]=mt(varargin{2:end});
    case 'gr'
        [varargout{1},varargout{2},varargout{3}]=makespectrogram(varargin{2:end});
    case 'gp'
        varargout={}; plotspectrogram(varargin{2:end});
    otherwise
        error('Requested function does not and has never existed!')
end



    function F=timefreq(varargin)
        x=varargin{1};
        if nargin>=2 && ~isempty(varargin{2}), fs=varargin{2}; else fs=1/(x(2)-x(1)); end
        maxfreq=fs/2.0; % Maximum frequency
        minfreq=fs/length(x); % Minimum and delta frequency -- simply the inverse of the length of the recording in seconds
        F=0:minfreq:maxfreq; % Create frequencies evenly spaced from 0:minfreq:maxfreq
    end




    function [P,F]=pl(varargin)
        
        x=varargin{1};
        if nargin>=2 && ~isempty(varargin{2}), rate=varargin{2};   else rate=200; end
        if nargin>=3 && ~isempty(varargin{3}), smooth=varargin{3}; else smooth=100; end
        if nargin>=4 && ~isempty(varargin{4}), pad=varargin{4};    else pad=0; end
        
        if noutputargs>1, F=timefreq(x,rate); else F=[]; end
        P=abs(fft(x)).^2/length(x);  % Raw, normalized FFT
        kernel=[0.25,0.5,0.25]; % Convolution kernel
        if pad>0 % Pad ends
            padded=zeros(length(P)+2*pad,1);
            padded(1:pad)=P(1);
            padded(end-pad:end)=P(end);
            padded(pad+1:end-pad)=P;
            if smooth>0, for q=1:smooth, padded=conv(padded,kernel,'same'); end, end % Convolve FFT with kernel
            P=padded(pad+1:end-pad);
        else
            if smooth>0, for q=1:smooth, P=conv(P,kernel,'same'); end, end % Convolve FFT with kernel
        end
        P=P(1:floor(length(x)/2)+1); % Trim so they're always the same size
    end
        



    function [P,F]=mt(varargin)
        x=varargin{1};
        if nargin>=2 && ~isempty(varargin{2}), rate=varargin{2};      else rate=200;    end
        if nargin>=3 && ~isempty(varargin{3}), bandwidth=varargin{3}; else bandwidth=4; end
        [P,F]=pmtm(x,bandwidth,[],round(rate));
    end






    function [P,F]=ar(varargin)
        
        x=varargin{1};
        if nargin>=2 && ~isempty(varargin{2}), rate=varargin{2};    else rate=200;  end
        if nargin>=3 && ~isempty(varargin{3}), order=varargin{3};   else order=10;  end
        if nargin>=4 && ~isempty(varargin{4}), maxfreq=varargin{4}; else maxfreq=0; end
        if nargin>=5 && ~isempty(varargin{5}), ntrials=varargin{5}; else ntrials=1; end
        
        % Set up frequency
        if maxfreq==0, F=timefreq(x,rate); % Define the frequency points
        else F=0:maxfreq; % ...or just pick them
        end
        
        % Figure out the length of each trial
        ptspertrial=length(x)/ntrials;
        if mod(length(x),ntrials), warning('spectrum:notdivisible','Data not divisible by proposed number of trials!'), keyboard, end
        
        x=reshape(x,1,[]); % For some reason armorf requires a row vector
        P=zeros(size(F)); % Initialize output -- same size as frequency
        [A,Z]=armorf(x,ntrials,ptspertrial,order); % Calculate autoregressive fit
        for i=1:length(F) % Loop over frequencies
            S=spectrum_AR(A,Z,order,F(i),rate); % Calculate spectrum
            P(i)=abs(S(1).^2); % Calculate and store power
        end
        
    end









    function [T,F,spectarr]=makespectrogram(varargin)
        
        ts=varargin{1};
        if nargin>=2 && ~isempty(varargin{2}), rate=varargin{2};    else rate=200;   end
        if nargin>=3 && ~isempty(varargin{3}), window=varargin{3};  else window=2;   end
        if nargin>=4 && ~isempty(varargin{4}), maxfreq=varargin{4}; else maxfreq=50; end
        if nargin>=5 && ~isempty(varargin{5}), tsmooth=varargin{5}; else tsmooth=2;  end
        if nargin>=6 && ~isempty(varargin{6}), fsmooth=varargin{6}; else fsmooth=2;  end
        
        % Handle input arguments
        if maxfreq==0 || maxfreq>rate/2, maxfreq=rate/2; end % Set maximum frequency if none entered or if too big
        npts=length(ts); % How long the data are
        maxtime=npts/rate; % Maximum time
        ts=ts-mean(ts); % Remove mean
        
        % Calculating spectra
        nwind=floor(maxtime/window); % Number of windows
        lwind=floor(window*rate); % Length of each window
        spectarr=zeros(lwind/2,nwind);
        for i=1:nwind
            tstart=lwind*(i-1)+1; % Initial time point
            tfinish=lwind*(i); % Final timepoint
            thists=ts(tstart:tfinish); % Pull out the part of the time series to make this spectrum
            tmparr=abs(fft(thists));
            spectarr(:,i)=tmparr(1:lwind/2);
        end
        
        % Perform smoothing
        tblur=[0.25,0.5,0.25];
        fblur=tblur;
        for i=1:nwind % Smooth in frequency
            for j=1:fsmooth
                spectarr(:,i)=conv(spectarr(:,i),fblur,'same');
            end
        end
        for i=1:lwind/2 % Smooth in time
            for j=1:tsmooth
                spectarr(i,:)=conv(spectarr(i,:),tblur,'same');
            end
        end
        
        % Calculate time and frequency limits
        finalfreq=floor(window*maxfreq);
        F=(0:finalfreq-1)/window;
        T=(0:nwind-1)*window;
        spectarr=spectarr(1:finalfreq,:);
        
    end
        





    function plotspectrogram(varargin)
        
        [F,T,spectarr]=makespectrogram(varargin{:});
        
        disp('Plotting...')
        surf(T,F,spectarr)
        xlabel('Time (s)')
        ylabel('Frequency (Hz)')
        
        disp('Done.')
        
    end

        


end

