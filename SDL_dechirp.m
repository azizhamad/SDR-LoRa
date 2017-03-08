close all
clear all

% Read file
[in, Fs] = audioread('SDRSharp_20170301_172427Z_868712500Hz_IQ_125k.wav');

% Allocate in-phase and quadrature components
x = (in(:,1) + 1i*in(:,2))';

% Spectrogram parameters
Nx = length(x);
window_length = 1024;
Nfft = 1024;
% Nfft = 2^SF;

% Crop signal in time
n = 2.306*Fs:4.13*Fs;
x = x(n);

% Test to bring signal to baseband
t = 0:1/Fs:length(x)/Fs-1/Fs;
x = x.*cos(2*pi*1.577e6*t);

% LoRa parameters
BW = 125e3;
SF = 6;
chip_rate = BW/2^SF;
Fs_prime = 2*BW;

% Chirp generation
fo = 0;
f1 = BW;
delta_t = 32.8e-3;
t = 0:1/Fs:length(x)/Fs;
c = chirp(t,fo,delta_t,f1);

% Decimation
% x = decimate(x, 9);

% De-chirping
% de_chirped = x.*conj(c);
de_chirped = x;


% Plot spectrogram
[s, f, t] = spectrogram(de_chirped, blackman(window_length), window_length/2, Nfft);
surf(f,t,10*log10(abs(s')),'EdgeColor','none')
axis xy; axis tight; colormap(jet); view(0,90);
ylabel('Time');
xlabel('Frequency (Hz)');
% xlim([0 2.3e6])
% ylim([0.001 0.81])
% printfigure('Caputred LoRa signal')

% % Plot spectrum
% X = fft(x, Nfft)/Nfft;
% f = Fs*linspace(0,1,Nfft);
% plot(f,10*log10(abs(X)))

% % Test to crop signal in frequency
% X = X(.5781/(2*pi)*Nfft:.6875/(2*pi)*Nfft);
% x = ifft(X);
