%% Load the audio file
% Read the original audio file to be processed
[original_audio, sample_rate] = audioread('handel_J.wav');

%% Add Noise to Original File

% Define noise frequencies in Hz and amplitudes
noise_freqs = [1500, 3060, 3510, 3569.9, 3920, 5180, 5679.9, 6040, 6129.9, 6369.9];
noise_amps = 0.1 * ones(size(noise_freqs));  % Set amplitude for each noise component

% Generate noise signal based on specified frequencies and amplitudes
time_vector = (0:length(original_audio)-1) / sample_rate; % Time vector for signal
noise_signal = zeros(size(original_audio));  % Initialize noise signal
for i = 1:length(noise_freqs)
    noise_signal = noise_signal + noise_amps(i) * sin(2 * pi * noise_freqs(i) * time_vector');
end

% Combine noise with the original audio signal
noisy_audio = original_audio + noise_signal;

% Save the noisy audio to a file
audiowrite('NoisyAudio.wav', noisy_audio, sample_rate);

% Calculate the Fourier Transform to identify spectral locations of noise
N = length(noisy_audio);
noisy_audio_fft = fft(noisy_audio);
frequency_vector = (0:N-1) * (sample_rate / N); % Frequency vector for FFT

% Find noise peaks in the positive half of the spectrum
[~, noise_peaks] = findpeaks(abs(noisy_audio_fft((N/2+1):end)), frequency_vector((N/2+1):end), 'MinPeakHeight', 0.01);

% Display identified noise frequencies
disp("Frequencies of Noise (Hz):")
disp(noise_peaks);

disp("Discrete Normalized Frequencies of Noise:")
disp(noise_peaks / sample_rate);

%% START Filtering NoisyAudio

% Reload the noisy audio file for filtering steps
[noisy_audio, sample_rate] = audioread('NoisyAudio.wav');

% Define radii for IIR filters
radius_1 = 0.4;
radius_2 = 0.95;
L = length(noisy_audio);
noisy_audio_fft = fftshift(fft(noisy_audio)) / L;
frequency_vector = -sample_rate/2:sample_rate/L:sample_rate/2 - sample_rate/L;

%% A: Plot magnitude and phase response of audio file
figure
subplot(211);
plot(frequency_vector, abs(noisy_audio_fft));
ylabel('Magnitude Response')
xlabel('Hz')

% Plot phase response
subplot(212);
plot(frequency_vector, unwrap(angle(noisy_audio_fft))); 
ylabel('Phase Response')
xlabel('Hz')

%% B: Determine Spectral location associated with noise
[~, noise_peaks] = findpeaks(abs(noisy_audio_fft((length(noisy_audio_fft)+1)/2:end)), frequency_vector((length(noisy_audio_fft)+1)/2:end), 'MinPeakHeight', 0.01);

disp("Frequencies of Noise (Hz):")
disp(noise_peaks);
disp("Discrete Normalized Frequencies of Noise:")
disp(noise_peaks / sample_rate);

%% C: Design an FIR filter to remove noise
syms z
fir_numerator = 1;  % FIR filter numerator (product of terms for each noise frequency)
for i = 1:length(noise_peaks)
   fir_numerator = fir_numerator * (z^2 - 2 * cos(2 * pi * noise_peaks(i) / sample_rate) * z + 1);
end
fir_coefficients = sym2poly(fir_numerator);

% Display FIR filter properties
plotMagnitude(fir_coefficients, 1, "FIR Filter");
plotPoleZero(fir_coefficients, 1, "FIR Filter");

% Apply FIR filter
filtered_fir_audio = filter(fir_coefficients, 1, noisy_audio);

% %% D: Design an IIR filter with radius_1
% iir_denom_1 = 1;
% for i = 1:length(noise_peaks)
%    iir_denom_1 = iir_denom_1 * (z^2 - 2 * radius_1 * cos(2 * pi * noise_peaks(i) / sample_rate) * z + radius_1^2);
% end
% iir_coefficients_1 = sym2poly(iir_denom_1);
% 
% % Display IIR filter properties (radius = 0.4)
% plotMagnitude(fir_coefficients, iir_coefficients_1, "IIR Filter with r = 0.4");
% plotPoleZero(fir_coefficients, iir_coefficients_1, "IIR Filter with r = 0.4");
% 
% % Apply IIR filter with radius_1
% filtered_iir_audio_1 = filter(fir_coefficients, iir_coefficients_1, noisy_audio);

%% E: Design an IIR filter with radius_2
iir_denom_2 = 1;
for i = 1:length(noise_peaks)
   iir_denom_2 = iir_denom_2 * (z^2 - 2 * radius_2 * cos(2 * pi * noise_peaks(i) / sample_rate) * z + radius_2^2);
end
iir_coefficients_2 = sym2poly(iir_denom_2);

% Display IIR filter properties (radius = 0.95)
plotMagnitude(fir_coefficients, iir_coefficients_2, "IIR Filter with r = 0.95");
plotPoleZero(fir_coefficients, iir_coefficients_2, "IIR Filter with r = 0.95");

% Apply IIR filter with radius_2
filtered_iir_audio_2 = filter(fir_coefficients, iir_coefficients_2, noisy_audio);

%% Save the new audio files
audiowrite("FIRaudio_filtered.wav", filtered_fir_audio / max(abs(filtered_fir_audio(:))) * 0.5, sample_rate)
% audiowrite("IIR1audio_filtered.wav", filtered_iir_audio_1 / max(abs(filtered_iir_audio_1(:))) * 0.5, sample_rate)
audiowrite("IIR2audio_filtered.wav", filtered_iir_audio_2 / max(abs(filtered_iir_audio_2(:))) * 0.5, sample_rate)



function plotMagnitude(num,den,str)
w = 0:0.001:2*pi;
numw = 0;
denw = 0;
for i = 1:length(num)
    numw = numw + num(i).*exp(1j*w*(i-1));
end
for i = 1:length(den)
    denw = denw + den(i).*exp(1j*w*(i-1));
end

figure
subplot(211)
Hw = numw./denw;
plot(w/pi,abs(Hw))
xlabel('Frequency *pi')
title("Magnitude Response of " + str)

subplot(212)
plot(w/pi,angle(Hw))
xlabel('Frequency *pi')
title("Phase Response of " + str)
end

function [poleG,zeroG] = plotPoleZero(num,den,str)
syms z
poleG = roots(den);
zeroG = zeros(length(den)-1,1); 
zeroG(1:length(num)-1,1) = roots(num);
disp("Zeros:")
disp(zeroG');
pause(1)
disp("Poles:");
disp(poleG');
pause(1)
denz = 1;
numz = 1;
for i = 1:length(zeroG)
    numz = numz * ( 1 - zeroG(i)*z^-1);
end
for j = 1:length(poleG)
    denz = denz * ( 1 - poleG(j)*z^-1);
end
disp("Factored G(z)")
pretty(numz);
if den ~= 1
    disp("----------------------------------------------------------------------------------------------------")
    pretty(denz);
end
s = tf('s');
Gs = tf(num,den);
figure
pzplot(Gs*s^(length(den)-length(num)));
title("Zero-Pole Plot of " + str)
end
