clear all
close all
clc

thisFile = matlab.desktop.editor.getActiveFilename;
thisFolder = fileparts(thisFile);

repo_folder = fileparts(thisFolder); 

addpath(genpath(repo_folder));

%% Create test signal

fs = 1000;              % Sampling frequency (1 kHz)
t = 0:1/fs:10-1/fs;     % 10 seconds of data

% 1. A stable 20 Hz Beta oscillation
beta_tone = 1.5 * sin(2*pi*20*t);

% 2. A frequency chirp (moving from 5 Hz to 45 Hz)
% This tests the temporal tracking of your spectrogram
chirp_sig = 2.0 * chirp(t, 5, 10, 45); 

% 3. Additive White Gaussian Noise (SNR test)
noise = 1.2 * randn(size(t));

% Combine them
test_signal = beta_tone + chirp_sig + noise;

% Visualization of the raw trace
plot(t(1:500), test_signal(1:500));
title('Synthetic LFP Test Signal (First 500ms)');
xlabel('Time (s)'); ylabel('Amplitude');

%% Welch PSD

[Pxx_welch,freq] = pwelch(test_signal,1*fs,round(0.5*fs),2^nextpow2(length(test_signal)),fs);

figure;
plot(freq,Pxx_welch);
xlim([0,100])

%% Multi-taper spectrum

tic
params = [];
params.pad = 0;
params.tapers = [4,7];
params.Fs = fs;
[Pxx_mt,freq] = mtspectrumc(test_signal,params);
toc

hold on;
plot(freq,Pxx_mt);
xlim([0,100])

%% Spectrogram
params.tapers = [3,5];
params.fpass = [0,50];
tic
[S,t,f] = mtspecgramc(test_signal,[1,0.5],params);
toc

figure; 
imagesc(t,f,S')
