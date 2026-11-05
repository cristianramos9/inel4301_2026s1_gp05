% inel4301_s086_gp05_sn05_tk05_mca03
% Modified by: Cristian N. Ramos Ramos
% Original by: Prof. Domingo Antonio Rodriguez
% CUANTIZADOR_MIDRISE - N-level uniform mid-rise quantizer.
% MATLAB R2020b or later. No specialized toolboxes required.
% The last cell includes Xmax; reconstruction values are cell centers.
% DTMF WORD USED: dtmfau
clear; clc; close all;
fprintf('\nUNIFORM MID-RISE QUANTIZER\n');
File_type = input('Enter 0 for WAV or 1 for numeric TXT: ');
if ~isscalar(File_type) || ~ismember(File_type,[0 1])
    error('File type must be 0 or 1.');
end
if File_type == 0
    filename = input('Enter WAV filename (e.g., dtmfad.wav): ','s');
    [xm,Fs] = audioread(filename);
else
    filename = input('Enter numeric TXT filename: ','s');
    xm = load(filename);
    Fs = input('Enter sampling frequency Fs (Hz): ');
end
if ~isnumeric(xm) || isempty(xm) || any(~isfinite(xm(:)))
    error('Input signal must contain finite numeric samples.');
end
if ~isscalar(Fs) || ~isfinite(Fs) || Fs <= 0
    error('Sampling frequency must be positive and finite.');
end
if ~isvector(xm) && size(xm,2) > 1
    fprintf('Multichannel signal: averaging channels to mono.\n');
    xm = mean(xm,2);
end
xm = xm(:);
fprintf('Available samples: %d; Fs = %.12g Hz\n',numel(xm),Fs);
length_of_xm = input('Samples to process (Enter for all): ');
if isempty(length_of_xm), length_of_xm = numel(xm); end
if ~isscalar(length_of_xm) || ~isfinite(length_of_xm) || ...
        length_of_xm ~= floor(length_of_xm) || ...
        length_of_xm < 1 || length_of_xm > numel(xm)
    error('Length must be an integer from 1 to %d.',numel(xm));
end
xm_cons = xm(1:length_of_xm);		% Original signal: xm_cons
Xmin = min(xm_cons);
Xmax = max(xm_cons);
N = input('Enter N quantization levels (power of 2, N >= 2): ');
if ~isscalar(N) || ~isfinite(N) || N < 2 || ...
        N ~= floor(N) || 2^round(log2(N)) ~= N
    error('N must be a power of 2 and at least 2.');
end
if Xmax == Xmin
    % Degenerate signal: avoid dividing by zero.
    Delta = 0;
    index = zeros(size(xm_cons));
    X_Quantized = xm_cons;		% Quantized signal: xm_cons
else
    Delta = (Xmax-Xmin)/N;
    index = floor((xm_cons-Xmin)/Delta);
    index = min(max(index,0),N-1); % Includes Xmax in the last cell.
    X_Quantized = Xmin + (index+0.5)*Delta;		% Quantized signal: xm_cons
end
time = (0:length_of_xm-1)'/Fs;
figure('Name','Mid-rise uniform quantization');
subplot(2,1,1); plot(time,xm_cons); grid on;
xlabel('Time (s)'); ylabel('Amplitude'); title('Original signal');
subplot(2,1,2); plot(time,X_Quantized,'r'); grid on;
xlabel('Time (s)'); ylabel('Amplitude'); title('Mid-rise quantized signal');
MSE = mean((xm_cons-X_Quantized).^2);

% Cristian Ramos
% Fundamental signal metrics for original signal
x = xm_cons;
arithmetic_mean = mean(x);
total_sig_energy = sum(abs(x).^2);
avg_sig_power = mean(abs(x).^2);
rms_amplitude = sqrt(mean(abs(x).^2));
variance = var(x,1);
eucl_norm = norm(x,2);

% Fundamental signal metrics for quantized signal
y = X_Quantized;
q_delta = str2double(sprintf('%.7f', Delta));		%	USE round(x, 7, 'significant');
q_arithmetic_mean = str2double(sprintf('%.7f', mean(y)));
q_total_sig_energy = sum(abs(y).^2);
q_avg_sig_power = mean(abs(y).^2);
q_rms_amplitude = sqrt(mean(abs(y).^2));
q_variance = var(y,1);
q_eucl_norm = norm(y,2);

q_mse_n4 = 0.2;
q_mse_n8 = 0.4;
q_mse_n16 = 0.5;
q_mse_n32 = 0.6;
q_mse_n64 = 0.7;
q_mse_n128 = MSE;

% Mid-rise quantizer MSE graph
q_levels = [4 8 16 32 64 128]
q_mse_levels = [q_mse_n4 q_mse_n8 q_mse_n16 q_mse_n32 q_mse_n64 q_mse_n128];
q_mse_title = sprintf('Mid-rise quantizer MSE (%d samples)', length_of_xm);
figure('Name', 'Mid-rise quantized MSE');
plot(q_levels, q_mse_levels, '-o'); grid on;
xlabel('Number of quantization levels, N'); ylabel('Mean squared quantization error (MSE)'); title(q_mse_title);

fprintf('\nInput range: [%.12g, %.12g]\n',Xmin,Xmax);
fprintf('N = %d, bits/sample = %d, Delta = %.12g\n',N,round(log2(N)),Delta);
fprintf('Output range: [%.12g, %.12g]\n',min(X_Quantized),max(X_Quantized));
fprintf('MSE = %.12g\n',MSE);
if Delta > 0
    fprintf('Highest reconstruction level: %.12g\n',Xmax-Delta/2);
end

% CSV File
% length_of_xm	-	SampleCount
% Fs		-	SamplingFrequency_Hz
% N		-	QuantizationLevels
% round(log2(N)) -	BitsPerSample
% Delta		-	QuantizationStep
% arithmetic_mean -	ArithmeticMean
% total_sig_energy -	TotalEnergy
% avg_sig_power	-	AveragePower
% rms_amplitude	-	RMS
% variance	-	Variance
% eucl_norm	-	EuclideanNorm
% MSE		-	MSE
bits_per_sample = round(log2(N))
%r0_title = ["Signal" "SampleCount" "SamplingFrequency_Hz" "QuantizationLevels" "BitsPerSample" "QuantizationStep" "ArithmeticMean" "TotalEnergy" "AveragePower" "RMS" "Variance" "EuclideanNorm" "MSE"];
%r1_original = ['Original' length_of_xm Fs N bits_per_sample Delta arithmetic_mean total_sig_energy avg_sig_power rms_amplitude variance eucl_norm MSE];
%r2_modified = ['MidRiseQuantized' length_of_xm Fs N bits_per_sample Delta q_arithmetic_mean q_total_sig_energy q_avg_sig_power q_rms_amplitude q_variance q_eucl_norm MSE];

r1_original = [01 length_of_xm Fs N bits_per_sample Delta arithmetic_mean total_sig_energy avg_sig_power rms_amplitude variance eucl_norm MSE];
r2_modified = [02 length_of_xm Fs N bits_per_sample q_delta q_arithmetic_mean q_total_sig_energy q_avg_sig_power q_rms_amplitude q_variance q_eucl_norm MSE];
%csvwrite('inel4301_s086_gp05_sn05_tk05_N128_mca03.csv', modified);

%csv_output = [r0_title; r1_original; r2_modified];
csv_output = [r1_original; r2_modified];

%writematrix(csv_output, 'inel4301_s086_gp05_sn05_tk05_N128_mca03.csv');
%dlmwrite('inel4301_s086_gp05_sn05_tk05_N128_mca03.csv', csv_output);

Ask = input('Play quantized audio? 1=yes, other=no: ');
if isequal(Ask,1), sound(X_Quantized,Fs); end
