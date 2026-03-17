function [S,t,f,Serr]=mtspecgramc_win(data,idx_win,params)
%MTSPECGRAMC Multi-taper spectrogram for continuous data.
%
%   [S,t,f,Serr] = mtspecgramc(data,movingwin,params)
%
%   Computes the time-frequency spectrum of a continuous signal by applying
%   Chronux multi-taper spectral estimation within sliding windows.
%
%   INPUTS
%       data
%           Continuous data arranged as:
%               samples x channels/trials
%
%       movingwin
%           Two-element vector specifying the sliding-window configuration:
%               [window_length  step_size]
%           The units must be consistent with params.Fs. For example, if
%           movingwin is given in seconds, then Fs must be in Hz.
%
%       params
%           Structure of analysis parameters with fields:
%
%           tapers
%               Either precomputed DPSS tapers, or one of the following:
%
%               (1) [TW K]
%                   TW : time-bandwidth product
%                   K  : number of tapers to use, with K <= 2*TW - 1
%
%               (2) [W T p]
%                   W : bandwidth
%                   T : duration of the data
%                   p : integer such that 2*TW - p tapers are used
%
%                   In this form, T must equal movingwin(1). Units of W and
%                   T must be consistent with each other and with params.Fs.
%
%               Default: form (1) with TW = 3 and K = 5
%
%           pad
%               Padding factor for the FFT. Allowed values are -1,0,1,2,...
%               -1 : no padding
%                0 : pad to next power of 2
%                1 : pad to twice the next power of 2, etc.
%               Default: 0
%
%           Fs
%               Sampling frequency. Default: 1
%
%           fpass
%               Frequency band of interest in the form [fmin fmax].
%               Default: [0 Fs/2]
%
%           err
%               Error bar option:
%                   [1 p] : theoretical error bars
%                   [2 p] : jackknife error bars
%                   [0 p] or 0 : no error bars
%               Default: 0
%
%           trialave
%               If 1, average across channels/trials.
%               If 0, return separate spectra for each channel/trial.
%               Default: 0
%
%   OUTPUTS
%       S
%           Spectrogram.
%           If trialave = 0:
%               time x frequency x channels/trials
%           If trialave = 1:
%               time x frequency
%
%       t
%           Time vector corresponding to the center of each analysis window.
%
%       f
%           Frequency vector.
%
%       Serr
%           Error bars, returned only when err(1) >= 1.
%
%   NOTES
%       - This function uses sliding windows and calls mtspectrumc on each
%         windowed segment.
%       - Window centers are reported in the same time units implied by Fs.

if nargin < 2; error('Need data and window parameters'); end;
if nargin < 3; params=[]; end;

% Parse analysis parameters and fill in defaults where needed.
[tapers,pad,Fs,fpass,err,trialave,params]=getparams(params);

[nw,Nwin] = size(idx_win);

T_win = Nwin/Fs;

% Validate that the taper duration matches the moving-window duration when
% tapers are supplied in the [W T p] form.
if length(params.tapers)==3 & T_win~=params.tapers(2);
    error('Duration of data in params.tapers is inconsistent with movingwin(1), modify params.tapers(2) to proceed')
end

% Error bars cannot be requested if err(1) indicates no error computation.
if nargout > 3 && err(1)==0; 
%   Cannot compute error bars with err(1)=0. change params and run again.
    error('When Serr is desired, err(1) has to be non-zero.');
end;

% Ensure the data are arranged as samples x channels/trials.
data=change_row_to_column(data);

% Extract data dimensions.
[N,Ch]=size(data);

% Compute FFT length, allowing optional zero-padding.
nfft=max(2^(nextpow2(Nwin)+pad),Nwin);

% Construct the frequency grid restricted to the requested passband.
f=getfgrid(Fs,nfft,fpass); 
Nf=length(f);

% Generate/check DPSS tapers for the requested window length.
params.tapers=dpsschk(tapers,Nwin,Fs); % check tapers


% Preallocate output arrays based on whether spectra are averaged across
% channels/trials.
if trialave
    S = zeros(nw,Nf);
    if nargout==4; Serr=zeros(2,nw,Nf); end
else
    S = zeros(nw,Nf,Ch);
    if nargout==4; Serr=zeros(2,nw,Nf,Ch); end
end

% Loop over windows and compute the multi-taper spectrum for each segment.
for n=1:nw
   indx=idx_win(n,:);%winstart(n):winstart(n)+Nwin-1;
   datawin=data(indx,:);
   
   if nargout==4
     % Compute spectrum and corresponding error bars for this window.
     [s,f,serr]=mtspectrumc(datawin,params);
     Serr(1,n,:,:)=squeeze(serr(1,:,:));
     Serr(2,n,:,:)=squeeze(serr(2,:,:));
   else
     % Compute spectrum only.
     [s,f]=mtspectrumc(datawin,params);
   end
   
   % Store the spectrum for the current window.
   S(n,:,:)=s;
end

% Remove singleton dimensions introduced during preallocation.
S=squeeze(S); 
if nargout==4;Serr=squeeze(Serr);end

% Compute the center time of each window in seconds (or consistent Fs units).
winmid=(idx_win(:,1) + idx_win(:,end))/2;%winstart+round(Nwin/2);
t=winmid/Fs;