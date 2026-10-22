% inel4301_syyy_gpzz_sn00_tk00_mca03
% Modified by:
% Original by: Prof. Domingo Antonio Rodriguez
% CUANTIZADOR_MIDRISE - N-level uniform mid-rise quantizer.
% MATLAB R2020b or later. No specialized toolboxes required.
% The last cell includes Xmax; reconstruction values are cell centers.
% DTMF WORD USED:
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
xm_cons = xm(1:length_of_xm);
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
    X_Quantized = xm_cons;
else
    Delta = (Xmax-Xmin)/N;
    index = floor((xm_cons-Xmin)/Delta);
    index = min(max(index,0),N-1); % Includes Xmax in the last cell.
    X_Quantized = Xmin + (index+0.5)*Delta;
end
time = (0:length_of_xm-1)'/Fs;
figure('Name','Mid-rise uniform quantization');
subplot(2,1,1); plot(time,xm_cons); grid on;
xlabel('Time (s)'); ylabel('Amplitude'); title('Original signal');
subplot(2,1,2); plot(time,X_Quantized,'r'); grid on;
xlabel('Time (s)'); ylabel('Amplitude'); title('Mid-rise quantized signal');
MSE = mean((xm_cons-X_Quantized).^2);
fprintf('\nInput range: [%.12g, %.12g]\n',Xmin,Xmax);
fprintf('N = %d, bits/sample = %d, Delta = %.12g\n',N,round(log2(N)),Delta);
fprintf('Output range: [%.12g, %.12g]\n',min(X_Quantized),max(X_Quantized));
fprintf('MSE = %.12g\n',MSE);
if Delta > 0
    fprintf('Highest reconstruction level: %.12g\n',Xmax-Delta/2);
end
Ask = input('Play quantized audio? 1=yes, other=no: ');
if isequal(Ask,1), sound(X_Quantized,Fs); end
