close all
clear all

% Read file
[in, Fs] = audioread('SDRSharp_20170301_172427Z_868712500Hz_IQ_125k.wav');

% Allocate in-phase and quadrature components
x = (in(:,2) + 1i*in(:,1)).';

% Spectrogram parameters
Nx = length(x);
window_length = 1024;
Nfft = 1024;

% Crop signal in time
n = 2.306*Fs:4.13*Fs; % 2.3055
x = x(n);

% Bring signal to baseband
t = 0:1/Fs:length(x)/Fs-1/Fs;
x = x.*cos(2*pi*1.577e6*t); % 1.456e6

% LoRa parameters
BW = 125e3;
SF = 6;
% Nfft = 2^SF;
chirp_rate = BW/2^SF;
Fs_prime = 2*BW;

% Chirp generation
fo = BW; % reverse go and f1 for an up-chirp
f1 = 0;
delta_t = 2^(2*SF)/BW; % 32.8e-3
t = 0:1/Fs:delta_t;
c = chirp(t,fo,delta_t,f1); % generate a down-chirp segment and concatenate multiple ones
c = [c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c];
c = c(1:length(x)); % crop to signal length
% Curently there is a mismatch in chirp BW, so the resulting spectrum
% is split in two,  half of it (non-consecutive pieces) is shifted to the right

% Decimation
% x = decimate(x, 9);

% De-chirping
de_chirped = x.*conj(c);
% de_chirped = c;


% Plot spectrogram
[s, f, t] = spectrogram(de_chirped, blackman(window_length), window_length/2, Nfft, Fs);
surf(f,t,10*log10(abs(s')),'EdgeColor','none')
axis xy; axis tight; colormap(jet); view(0,90);
ylabel('Time');
xlabel('Frequency (Hz)');

xlim([0 2*BW])
% ylim([0.001 0.81])
% printfigure('Caputred LoRa signal')

% % Plot spectrum
% X = fft(x, Nfft)/Nfft;
% f = Fs*linspace(0,1,Nfft);
% plot(f,10*log10(abs(X)))

% % Test to crop signal in frequency
% X = X(.5781/(2*pi)*Nfft:.6875/(2*pi)*Nfft);
% x = ifft(X);
