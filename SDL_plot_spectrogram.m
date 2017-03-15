close all
clear all

% Read file
[in, Fs] = audioread('SDRSharp_20170301_171509Z_868787500Hz_IQ_125k.wav');

% Allocate in-phase and quadrature components
x = (in(:,2) + 1i*in(:,1)).';

% Crop signal in time
n = 2.166*Fs:3.485*Fs;
x = x(n);

% Spectrogram parameters
window_length = 1024;
Nfft = window_length;

% SPECTROGRAM OPTION 1
% spectrogram(x, blackman(window_length), window_length/2, Nfft, Fs);

% SPECTROGRAM OPTION 2
[s, f, t] = spectrogram(x, blackman(window_length), window_length/2, Nfft, Fs);
% f = f(299:352); % crop frequency
% s = s(299:352,:); % crop spectrogram
surf(f,t,10*log10(abs(s.')),'EdgeColor','none')
axis xy; axis tight; colormap(jet); view(0,90);
ylabel('Time');
xlabel('Frequency (Hz)');
% xlim([1.3e6 2.3e6])
% ylim([0.001 0.81])