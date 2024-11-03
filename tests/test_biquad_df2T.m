clc; close all; clearvars;
% biquad_df2T_example - Example as for the biquad_df2T class. Stimulus is
% generated and processed in blocks by the biquad_df2T class that uses
% the specific biquad structure filter called "direct form 2 transposed".
% Before processing, the biquad filter type and parameters are defined.
% Output is then stored in a buffer. Last step is plotting modulus and
% phase response of the filter and input and output signals in time domain.
% --------------------------
% Author:  Oscar Butler
% Project: MBiquad
% Date:    11.10.2023
% --------------------------

%% General settings
fs = 48e3;
blocksize = 128; % Buffer size

%% Stimulus
N = 1024;
Nfft = 2^14;
fSine = 1000;
t = 0:1/fs:(N-1)/fs;
f = 0:fs/Nfft:fs-1;
x = sin(2*pi*fSine*t);

%% Biquad parameters
filter_param.type = 0;           % Lowpass filter type
filter_param.numStages = 2;      % Biquad order
filter_param.gaindB = 0;         % Gain in dB
filter_param.freqCut = 2500.0;   % Significant frequency
filter_param.Q = 0.5;            % Quality factor

%% Array init
inBuffer = zeros(1,blocksize);

%% Init of LS filter instance
filter = biquad_df2T(filter_param, fs, blocksize);
H_filter = freqz(filter.coeffs(1:3),[1 filter.coeffs(4:5)],f,fs);

%% Block processing
for i=0:N/blocksize-1
    inBuffer = x(i*blocksize+1:(i+1)*blocksize);
    filter.mono_df2T(inBuffer);
    y(i*blocksize+1:(i+1)*blocksize) = filter.outputBuffer;
end

%% Plot filter modulus and phase responses, stimulus and filtered stimulus
figure;
subplot(221)
semilogx(f,20*log10(abs(H_filter)));
grid on;
xlim([10 20e3])
xlabel('Frequency [Hz]')
ylabel('Magnitude [dB]')

subplot(222)
semilogx(f,180/pi*angle(H_filter));
grid on;
xlim([10 20e3])
xlabel('Frequency [Hz]')
ylabel('Phase [deg]')

subplot(2,2,[3 4])
plot(x)
hold on;
plot(y)
grid on;
xlabel('Samples')
ylabel('Amplitude')
xlim([0 N])

sgtitle('Filter frequency response and in/out signals in time domain')









