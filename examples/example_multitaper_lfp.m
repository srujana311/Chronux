clear all
close all
clc

thisFile = matlab.desktop.editor.getActiveFilename;
thisFolder = fileparts(thisFile);

repo_folder = fileparts(thisFolder); 

addpath(genpath(repo_folder));

project_folder = '/Users/lilychamakura/Documents/Lab_related/GitHub/XAI_project';
addpath([project_folder,'/Subfunctions'])

%% Load data

data_folder = [project_folder,'/DBSTRD010/Preprocessing/data/BlankThr_15'];
filename = 'Preprocessed_rBipolar_Blk538.mat';

load([data_folder,'/',filename],'X_preproc','cfg_preproc');

fs = cfg_preproc.info.seeg.sampling_freq_Hz;              % Sampling frequency (1 kHz)

[nPts,n_seeg] = size(X_preproc);
cnt_chan = randi(n_seeg);
test_signal = zscore(X_preproc(:,cnt_chan));

%% Welch PSD

% [Pxx_welch,freq] = pwelch(test_signal,1*fs,round(0.5*fs),2^nextpow2(length(test_signal)),fs);
% 
% figure;
% plot(freq,Pxx_welch);
% xlim([0,100])

%% Multi-taper spectrum

% tic
% params = [];
% params.pad = 0;
% params.tapers = [4,7];
% params.Fs = fs;
% [Pxx_mt,freq] = mtspectrumc(test_signal,params);
% toc

% hold on;
% plot(freq,Pxx_mt);
% xlim([0,100])

%% Spectrogram

params = [];
params.pad = 0;
params.Fs = fs;
params.tapers = [3,5];
params.fpass = [0,100];

T_win = 1;
T_step = T_win/4;

tic
[S,t,f] = mtspecgramc(X_preproc,[T_win,T_step],params);
toc

figure; 
imagesc(t,f,log(S(:,:,cnt_chan)'))
% %%
band_limits = [1 4;
               4 7;
               8 12;
               13 30;
               31 55];

S_mn = [];
for cnt_band = 1:size(band_limits,1)

    f_low = band_limits(cnt_band,1);
    f_hi = band_limits(cnt_band,2);

    flag = f >= f_low & f <= f_hi;

    S_mn(:,cnt_band) = mean(S(:,flag,cnt_chan),2);

end

figure; plot(log(S_mn(:,3)))

%%

filename = 'HilbertEnv_Alpha_rBipolar_Blk538.mat';

load([data_folder,'/',filename],'X_envelope');

N_win=round(fs*T_win); % number of samples in window
N_step=round(T_step*fs); % number of samples to step through (50% overlap)
% nfft=max(2^(nextpow2(N_win)+pad),Nwin);
% f=getfgrid(Fs,nfft,fpass); Nf=length(f);

N = size(X_envelope,1);
winstart=1:N_step:N-N_win+1;
nw=length(winstart); 

S_hilbert = [];
for n=1:nw
   indx=winstart(n):winstart(n)+N_win-1;
   datawin=X_envelope(indx,:);
    S_hilbert(n,:) = mean(datawin.^2,1);
end

%%

% t_start = randi(550);
% time_window = [t_start t_start + 30]
% figure;
subplot(211)
plot_windowed_timeseries(S_hilbert(:,cnt_chan),1/T_step,time_window); hold on;

subplot(212)
plot_windowed_timeseries(S_mn(:,3),1/T_step,time_window); hold on;

% figure;
% 
% plot_windowed_timeseries(S_hilbert(:,cnt_chan),1/T_step); hold on;
% plot_windowed_timeseries(10*S_mn(:,3),1/T_step);

%%

figure;
scatter(S_hilbert(:,cnt_chan),S_mn(:,3))