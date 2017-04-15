close all
clear

%% Recorded signal import
[in, Fs] = audioread('SDRSharp_20170301_172427Z_868712500Hz_IQ_125k.wav');
% Allocate in-phase and quadrature components
x = complex(in(:,2), in(:,1)).';
clear in

%% LoRa parameters
BW = 125e3;
SF = 12;
Fc = 1.6385e6; % 1.64e6
symbol_time = 2^SF/BW; % 32.8e-3
symbols_per_frame = 57;

%% Signal channelization (DDC)
% Comment: Use the Digital down converter or polyphase FFT bank
% Bring signal to baseband
t = 0:1/Fs:length(x)/Fs-1/Fs;
x = x.*cos(2*pi*(Fc-BW/2)*t);
% Filter the signal
% freqs = [0 2*BW/Fs 2*BW/Fs*9/8 1];
% damps = [1 1 0 0];
% order = 50;
% b = firpm(order,freqs,damps);
% x = filter(b,1,x);
% Decimation
x = resample(x, 2*BW, Fs); % Output sampling frequency is 2*BW
Fs = 2*BW;
clear t

%% Chirp generation
f0 = 0;
f1 = BW;
t = 0:1/Fs:symbol_time - 1/Fs;
upChirp = chirp(t,f0,symbol_time,f1); % generate a down-chirp segment and concatenate multiple ones
downChirp = chirp(t,f1,symbol_time,f0);
% downChirp = downChirp .* exp(-1i*pi/64*ones(1,length(downChirp))); % add a phase

upChirp = repmat(upChirp,1,10);

%% Signal conditioning and synchronization
% Find the start of the signal
[corr, lag] = xcorr(x, upChirp);
% corrThresh = 20;
corrThresh = max(abs(corr))/4;
cLag = find(abs(corr) > corrThresh, 1); % Find the first peak above threshold
signalStartIndex = abs(lag(cLag)) + 9*symbol_time*Fs;
signalEndIndex = round(signalStartIndex + symbols_per_frame*symbol_time*Fs);

% signalEndIndex = find(abs(corr) > corrThresh, length(corr));
% signalEndIndex = lag(signalEndIndex(end)) + 5*symbol_time*Fs;

% Synchronize SFD
symbol_offset = 0.25; % 12.25 to skip preamble and SFD
signalStartIndex = round(signalStartIndex + symbol_offset*symbol_time*Fs);
% Crop signal in time
% % signalStartIndex = 2.3055*Fs;
% % signalEndIndex = 4.11*Fs;
x = x(signalStartIndex:signalEndIndex);
clear lag corr

%% De-chirping
downChirp = repmat(downChirp,1,ceil(length(x)/length(downChirp)));
downChirp = downChirp(1:length(x));
% upChirp = repmat(upChirp,1,ceil(length(x)/length(upChirp)));
% upChirp = upChirp(1:length(x));
de_chirped = x.*downChirp;

%% Comment:
% The resulting spectrum is split in two,  half of it (non-consecutive
% pieces) is shifted to the right. This is why it's needed to
% artificially overlap the two halves.

%% Spectrogram computation
signal = de_chirped;
Nfft = 2^(SF+1); % +1 because the spectrum is doubled (Nyquist) % 1024
window_length = Nfft; % same as symbol_time*Fs;
[s, f, t] = spectrogram(signal, blackman(window_length), 0, Nfft, Fs);

%% Spectrogram conditioning
if isreal(signal)
    Nfft = Nfft/2+1;
end

s = circshift(s,Nfft*3/4,1);

% Overlapping option 1: wrong, should only consider first half of
% spectrogram
% s_first = s(1:round(BW/Fs*Nfft),:); %s(1:BW/Fs*Nfft,:);
% s_second = s(round(BW/Fs*Nfft)+1:round(BW/Fs*Nfft)*2,:); % s(BW/Fs*Nfft+1:BW/Fs*Nfft*2-1,:);
% padding = zeros(Nfft-round(BW/Fs*Nfft),size(s,2));
% s_second_padded = [s_second; padding];
% s = s + s_second_padded; % add padding

% Overlapping option 2: correct
% s_first = s(1:round(BW/Fs*Nfft),:);
% s_second = s(round(BW/Fs*Nfft)+1:round(BW/Fs*Nfft)*2,:);
% s = s_first + s_second;
% f = f(1:round(BW/Fs*Nfft));

%% Spectrogram plotting
surf(f,t,10*log10(abs(s.')),'EdgeColor','none')
axis xy; axis tight; colormap(jet); view(0,90);
ylabel('Time');
xlabel('Frequency (Hz)');
xlim([0 BW])

%% Bit extraction
s = s(1:Nfft/2,:); % just to make sure
s = s(:,1:symbols_per_frame-2);
[~, symbols] = max(abs(s));
symbols = mod(symbols - round(mean(symbols(1:8))), 2^SF);
bits =  dec2base(symbols, 2);


%% Other code
% printfigure('Caputred LoRa signal')

% % Plot spectrum
% X = fft(x, Nfft)/Nfft;
% f = Fs*linspace(0,1,Nfft);
% plot(f,10*log10(abs(X)))

% % Test to crop signal in frequency
% X = X(.5781/(2*pi)*Nfft:.6875/(2*pi)*Nfft);
% x = ifft(X);

%% Other chirp generation options
% Manual chirp generation
% k = (f1-f0)/symbol_time;
% upChirp = sin(-1*2*pi*(f0*t+k/2*t.^2));
% downChirp = sin(2*pi*(f0*t+k/2*t.^2));
% % downChirp = downChirp .* exp(-1i*pi/64*ones(1,length(downChirp))); % add a phase

% Manual complex chirp generation
% k = (f1-f0)/symbol_time;
% upChirp = -1*2*pi*(f0*t+k/2*t.^2);
% downChirp = 2*pi*(f0*t+k/2*t.^2);
% upChirp = cos(upChirp) + 1i * sin(upChirp);
% downChirp = cos(downChirp) + 1i * sin(downChirp);
% upChirp = (1+1i) * upChirp;
% downChirp = (1+1i) * downChirp;
% % c = c .* exp(-1i*pi/64*ones(1,length(c))); % add a phase

% Channelized chirp generation
% fo = BW; % reverse fo and f1 for an up-chirp
% f1 = 0;
% Fs_prime = 2*BW;
% symbol_time = 1/chirp_rate; % 32.8e-3
% t = 0:1/Fs_prime:symbol_time;
% c = chirp(t,fo,symbol_time,f1); % generate a down-chirp segment and concatenate multiple ones
% c = repmat(c,1,ceil(length(x)/length(c)));
% c = c(1:length(x));
% % c = c .* exp(-1i*pi*ones(1,length(x))); % add a phase
%%
