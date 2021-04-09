%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script test methods for frequency tagging on fake data
% Ana Flo January 2021
% Applies the same method but using data of sleeping 5-month-old infants,
% and the signal is simulated using the scripts provide by:
% https://data.mrc.ox.ac.uk/data-set/simulated-eeg-data-generator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all;

%% Declare variables
MethodToTest    = [1,1,1,1];
NbDat           = 40;          % Number of random data repetition
fs              = 500;           % sampling frequency
WordDuration    = 0.9;   % Word duration
duration        = 800*WordDuration;    % 800 repetitions of the word
epoch           = 12*WordDuration;          % epochs of 10s % --> I replaced 6 by 12
overlapp        = 11/12*epoch;     % overlapp in seconds of 5/6 % --> I modified it
binsSNR         = [(-3:-1) (1:3)]; % number of frequency bins on both sides to compute SNR
noise2signal    = 50;  % ratio between the noise and the signal (1/SNR)
folderout       = 'datatagging';
fmin            = 0.2;
fmax            = 6;
ftarget         = 1.11111*(1:6);
el              = 11;  % electrode to use

%% Load the data
load(fullfile('datanoise','dataNoise_fh2e-01.mat'))
dataNoise = dataNoise(1:NbDat,:);

%% Define the simulated data

% Syllabic response simulation
signalSyllfrq = 2.5;  % frequency of the cosine (width of the peak)
signalSyllpeak = (WordDuration*fs)/3*[1 2 3]-75;  % position of the peaks within a word
ERPsyll = zeros(1,WordDuration*fs);
for i=1:length(signalSyllpeak)
    ERPsyll = ERPsyll + createsignal(WordDuration*fs, fs, signalSyllfrq, signalSyllpeak(i));
end
ERPsyll = repmat(ERPsyll,[1 duration/WordDuration]);
dataSyll = dataNoise + 1/noise2signal * repmat(ERPsyll,[NbDat,1]);

% Syllabic+Word response simulation
signalSyllfrq = 2.5;  % frequency of the cosine (width of the peak)
signalSyllpeak = (WordDuration*fs)/3*[1 2 3]-75;  % position of the peaks within a word
signalWordfrq = 0.8;  % frequency of the cosine (width of the peak)
signalWordpeak = (WordDuration*fs)/3+10;  % position of the peaks within a word  (WordDuration*fs)/2-0.03*fs
ERPsyllword = zeros(1,WordDuration*fs);
for i=1:length(signalSyllpeak)
    ERPsyllword = ERPsyllword + createsignal(WordDuration*fs, fs, signalSyllfrq, signalSyllpeak(i));
end
for i=1:length(signalWordpeak)
    ERPsyllword = ERPsyllword + 2*createsignal(WordDuration*fs, fs, signalWordfrq, signalWordpeak(i));
end
ERPsyllword = repmat(ERPsyllword,[1 duration/WordDuration]);
dataSyllWord = dataNoise + 1/noise2signal * repmat(ERPsyllword,[NbDat,1]);


% Plot the simulated signal
n = 4;

figure('Position',[100 100 650 500],'Color',[1 1 1]),
subplot(2,1,1)
d = ERPsyll(1,1:n*WordDuration*fs);
time = (1:n*WordDuration*fs)/500;
plot(time,d,'k'), 
set(gca,'FontSize',12), 
xlabel('time (s)'), 
xlim([0 time(end)])
ylim([0 1.20])
title('Signal at the syllabic rate')

subplot(2,1,2)
d = ERPsyllword(1,1:n*WordDuration*fs);
time = (1:n*WordDuration*fs)/500;
plot(time,d,'k'), 
set(gca,'FontSize',12), 
xlabel('time (s)'), 
xlim([0 time(end)])
ylim([0 3.0])
title('Signal at the word rate')


%% Define the epochs

nsbj = size(dataNoise,1);
nsmpltot = size(dataNoise,2);
nsmpl = round(epoch*fs);
nsmplshift = round((epoch-overlapp)*fs);
nepI = floor(size(dataNoise,2)/round(epoch*fs));
iepO = (1:nsmplshift:(nsmpltot-nsmpl+1));
nepO = length(iepO);
   
% Noise - epochs with overlap
OEpNoise = nan(nsbj,nsmpl,nepO);
for i=1:nepO
    OEpNoise(:,:,i) = dataNoise(:,iepO(i):iepO(i)+nsmpl-1);
end

% Noise - epochs without overlap
IEpNoise = reshape(dataNoise(:,1:nepI*nsmpl),[nsbj nsmpl nepI]);

% Syll - epochs with overlap
OEpSyll = nan(nsbj,nsmpl,nepO);
for i=1:nepO
    OEpSyll(:,:,i) = dataSyll(:,iepO(i):iepO(i)+nsmpl-1);
end

% Syll - epochs without overlap
IEpSyll = reshape(dataSyll(:,1:nepI*nsmpl),[nsbj nsmpl nepI]);

% SyllWord - epochs with overlap
OEpSyllWord = nan(nsbj,nsmpl,nepO);
for i=1:nepO
    OEpSyllWord(:,:,i) = dataSyllWord(:,iepO(i):iepO(i)+nsmpl-1);
end

% Syll - epochs without overlap
IEpSyllWord = reshape(dataSyllWord(:,1:nepI*nsmpl),[nsbj nsmpl nepI]);

%% Method 1 : FFT power on each epoch
if MethodToTest(1)
    % Noise - On epochs with overlapp
    [pwrONoise, freq] = dft_power_eachepoch(OEpNoise, fs, ones(1,size(OEpNoise,3)),'fmin',fmin,'fmax',fmax);
    % Noise - On independant epochs
    [pwrINoise, freq] = dft_power_eachepoch(IEpNoise, fs, ones(1,size(IEpNoise,3)),'fmin',fmin,'fmax',fmax);
    % Syll - On epochs with overlapp
    [pwrOSyll, freq] = dft_power_eachepoch(OEpSyll, fs, ones(1,size(OEpSyll,3)),'fmin',fmin,'fmax',fmax);
    % Syll - On independant epochs
    [pwrISyll, freq] = dft_power_eachepoch(IEpSyll, fs, ones(1,size(IEpSyll,3)),'fmin',fmin,'fmax',fmax);
    % SyllWord - On epochs with overlapp
    [pwrOSyllWord, freq] = dft_power_eachepoch(OEpSyllWord, fs, ones(1,size(OEpSyllWord,3)),'fmin',fmin,'fmax',fmax);
    % SyllWord - On independant epochs
    [pwrISyllWord, freq] = dft_power_eachepoch(IEpSyllWord, fs, ones(1,size(IEpSyllWord,3)),'fmin',fmin,'fmax',fmax);
    
    % SNR
    pwrONoise_SNR = frqa_powernorm(pwrONoise, freq, 'binsnorm',binsSNR,'typeSNR','linear','ftarget',ftarget);
    pwrINoise_SNR = frqa_powernorm(pwrINoise, freq, 'binsnorm',binsSNR,'typeSNR','linear','ftarget',ftarget);
    pwrOSyll_SNR = frqa_powernorm(pwrOSyll, freq, 'binsnorm',binsSNR,'typeSNR','linear','ftarget',ftarget);
    pwrISyll_SNR = frqa_powernorm(pwrISyll, freq, 'binsnorm',binsSNR,'typeSNR','linear','ftarget',ftarget);
    pwrOSyllWord_SNR = frqa_powernorm(pwrOSyllWord, freq, 'binsnorm',binsSNR,'typeSNR','linear','ftarget',ftarget);
    pwrISyllWord_SNR = frqa_powernorm(pwrISyllWord, freq, 'binsnorm',binsSNR,'typeSNR','linear','ftarget',ftarget);
    
    
end

%% Method 2 : FFT power on averaged epochs
if MethodToTest(2)
    % Noise - On epochs with overlapp
    [pwravgONoise, freq] = dft_power(OEpNoise, fs, ones(1,size(OEpNoise,3)),'fmin',fmin,'fmax',fmax);
    % Noise - On independant epochs
    [pwravgINoise, freq] = dft_power(IEpNoise, fs, ones(1,size(IEpNoise,3)),'fmin',fmin,'fmax',fmax);
    % Syll - On epochs with overlapp
    [pwravgOSyll, freq] = dft_power(OEpSyll, fs, ones(1,size(OEpSyll,3)),'fmin',fmin,'fmax',fmax);
    % Syll - On independant epochs
    [pwravgISyll, freq] = dft_power(IEpSyll, fs, ones(1,size(IEpSyll,3)),'fmin',fmin,'fmax',fmax);
    % SyllWord - On epochs with overlapp
    [pwravgOSyllWord, freq] = dft_power(OEpSyllWord, fs, ones(1,size(OEpSyllWord,3)),'fmin',fmin,'fmax',fmax);
    % SyllWord - On independant epochs
    [pwravgISyllWord, freq] = dft_power(IEpSyllWord, fs, ones(1,size(IEpSyllWord,3)),'fmin',fmin,'fmax',fmax);
   
    % SNR
    pwravgONoise_SNR = frqa_powernorm(pwravgONoise, freq, 'binsnorm',binsSNR,'typeSNR','linear','ftarget',ftarget);
    pwravgINoise_SNR = frqa_powernorm(pwravgINoise, freq, 'binsnorm',binsSNR,'typeSNR','linear','ftarget',ftarget);
    pwravgOSyll_SNR = frqa_powernorm(pwravgOSyll, freq, 'binsnorm',binsSNR,'typeSNR','linear','ftarget',ftarget);
    pwravgISyll_SNR = frqa_powernorm(pwravgISyll, freq, 'binsnorm',binsSNR,'typeSNR','linear','ftarget',ftarget);
    pwravgOSyllWord_SNR = frqa_powernorm(pwravgOSyllWord, freq, 'binsnorm',binsSNR,'typeSNR','linear','ftarget',ftarget);
    pwravgISyllWord_SNR = frqa_powernorm(pwravgISyllWord, freq, 'binsnorm',binsSNR,'typeSNR','linear','ftarget',ftarget);
        
end

%% Method 3 : Phase Locking Value
if MethodToTest(3)
    % Noise - On epochs with overlapp
    [plvONoise, freq] = dft_plv(OEpNoise, fs, ones(1,size(OEpNoise,3)),'fmin',fmin,'fmax',fmax);
    % Noise - On Independants epochs
    [plvINoise, freq] = dft_plv(IEpNoise, fs, ones(1,size(IEpNoise,3)),'fmin',fmin,'fmax',fmax);
    % Syll - On epochs with overlapp
    [plvOSyll, freq] = dft_plv(OEpSyll, fs, ones(1,size(OEpNoise,3)),'fmin',fmin,'fmax',fmax);
    % Syll - On Independants epochs
    [plvISyll, freq] = dft_plv(IEpSyll, fs, ones(1,size(IEpNoise,3)),'fmin',fmin,'fmax',fmax);
    % SyllWord - On epochs with overlapp
    [plvOSyllWord, freq] = dft_plv(OEpSyllWord, fs, ones(1,size(OEpNoise,3)),'fmin',fmin,'fmax',fmax);
    % SyllWord - On Independants epochs
    [plvISyllWord, freq] = dft_plv(IEpSyllWord, fs, ones(1,size(IEpNoise,3)),'fmin',fmin,'fmax',fmax);

    % SNR (Z-score relative to adjacent frequencies)
    plvONoise_SNR = frqa_norm(plvONoise, freq, 'binsnorm',binsSNR,'typeSNR','mean','ftarget',ftarget);
    plvINoise_SNR = frqa_norm(plvINoise, freq, 'binsnorm',binsSNR,'typeSNR','mean','ftarget',ftarget);
    plvOSyll_SNR = frqa_norm(plvOSyll, freq, 'binsnorm',binsSNR,'typeSNR','mean','ftarget',ftarget);
    plvISyll_SNR = frqa_norm(plvISyll, freq, 'binsnorm',binsSNR,'typeSNR','mean','ftarget',ftarget);
    plvOSyllWord_SNR = frqa_norm(plvOSyllWord, freq, 'binsnorm',binsSNR,'typeSNR','mean','ftarget',ftarget);
    plvISyllWord_SNR = frqa_norm(plvISyllWord, freq, 'binsnorm',binsSNR,'typeSNR','mean','ftarget',ftarget);
   
end

%% method 4 ITC (Batterink and Paller)
if MethodToTest(4)
    frames = epoch*fs;
    tlimits = 1000*[1/fs,epoch];
    cycles = [1 45];
    freqs = [0.2,20.2];
    nfreqs = 200;
    nsbj = size(IEpNoise,1);
    % Noise - On epochs with overlapp
    itcwONoise = nan(nfreqs,nsbj);
    for i=1:size(OEpNoise,1)
        dat = squeeze(OEpNoise(i,:,:));
        [ersp,itc,powbase,times,freq,erspboot,itcboot] = newtimef(dat, frames, tlimits, fs, cycles,'freqs', freqs,'nfreqs',nfreqs,'plotersp','off','plotitc','off');
        itcwONoise(:,i) = mean(abs(itc), 2);
    end
    % Noise - On independant epochs
    itcwINoise = nan(nfreqs,nsbj);
    for i=1:size(IEpNoise,1)
        dat = squeeze(IEpNoise(i,:,:));
        [ersp,itc,powbase,times,freq,erspboot,itcboot] = newtimef(dat, frames, tlimits, fs, cycles,'freqs', freqs,'nfreqs',nfreqs,'plotersp','off','plotitc','off');
        itcwINoise(:,i) = mean(abs(itc), 2);
    end
    % Syll - On epochs with overlapp
    itcwOSyll = nan(nfreqs,nsbj);
    for i=1:size(OEpSyllWord,1)
        dat = squeeze(OEpSyll(i,:,:));
        [ersp,itc,powbase,times,freq,erspboot,itcboot] = newtimef(dat, frames, tlimits, fs, cycles,'freqs', freqs,'nfreqs',nfreqs,'plotersp','off','plotitc','off');
        itcwOSyll(:,i) = mean(abs(itc), 2);
    end
    % Syll - On independant epochs
    itcwISyll = nan(nfreqs,nsbj);
    for i=1:size(IEpSyll,1)
        dat = squeeze(IEpSyll(i,:,:));
        [ersp,itc,powbase,times,freq,erspboot,itcboot] = newtimef(dat, frames, tlimits, fs, cycles,'freqs', freqs,'nfreqs',nfreqs,'plotersp','off','plotitc','off');
        itcwISyll(:,i) = mean(abs(itc), 2);
    end
    % SyllWord - On epochs with overlapp
    itcwOSyllWord = nan(nfreqs,nsbj);
    for i=1:size(OEpSyllWord,1)
        dat = squeeze(OEpSyllWord(i,:,:));
        [ersp,itc,powbase,times,freq,erspboot,itcboot] = newtimef(dat, frames, tlimits, fs, cycles,'freqs', freqs,'nfreqs',nfreqs,'plotersp','off','plotitc','off');
        itcwOSyllWord(:,i) = mean(abs(itc), 2);
    end
    % SyllWord - On independant epochs
    itcwISyllWord = nan(nfreqs,nsbj);
    for i=1:size(IEpSyllWord,1)
        dat = squeeze(IEpSyllWord(i,:,:));
        [ersp,itc,powbase,times,freq,erspboot,itcboot] = newtimef(dat, frames, tlimits, fs, cycles,'freqs', freqs,'nfreqs',nfreqs,'plotersp','off','plotitc','off');
        itcwISyllWord(:,i) = mean(abs(itc), 2);
    end
    
end
