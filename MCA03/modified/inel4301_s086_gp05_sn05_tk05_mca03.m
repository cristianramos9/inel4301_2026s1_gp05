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
    X_Quantized = xm_cons;		% Quantized signal: X_Quantized
else
    Delta = (Xmax-Xmin)/N;
    index = floor((xm_cons-Xmin)/Delta);
    index = min(max(index,0),N-1); % Includes Xmax in the last cell.
    X_Quantized = Xmin + (index+0.5)*Delta;		% Quantized signal: X_Quantized 
end
time = (0:length_of_xm-1)'/Fs;
figure('Name','Mid-rise uniform quantization');

subplot(2,1,1); plot(time,xm_cons); grid on;
xlabel('Time (s)'); ylabel('Amplitude'); title('Original signal');
subplot(2,1,2); plot(time,X_Quantized,'r'); grid on;
xlabel('Time (s)'); ylabel('Amplitude'); title('Mid-rise quantized signal');
MSE = mean((xm_cons-X_Quantized).^2);

%==========================================================================
% Cristian Ramos
%==========================================================================

%==========================================================================
% DECLARING VARIABLES
%==========================================================================
% Fundamental signal metrics for original signal
%--------------------------------------------------------------------------
x = xm_cons;
orig_delta = "NaN";
orig_N = "NaN";
orig_bits_per_sample = "NaN";
orig_mse = "NaN";
arithmetic_mean = str2double(sprintf('%9f', mean(x)));
total_sig_energy = str2double(sprintf('%7f', sum(abs(x).^2)));
avg_sig_power = str2double(sprintf('%.9f', mean(abs(x).^2)));
rms_amplitude = str2double(sprintf('%.9f', sqrt(mean(abs(x).^2))));
variance = str2double(sprintf('%.8f', var(x,1)));
eucl_norm = str2double(sprintf('%.9f', norm(x,2)));

% Fundamental signal metrics for quantized signal
%--------------------------------------------------------------------------
y = X_Quantized;
q_delta = round(str2double(sprintf('%.7g', Delta)));		%	USE round(x, 7, 'significant');
q_bits_per_sample = round(log2(N));
q_arithmetic_mean = str2double(sprintf('%.9f', mean(y)));
q_total_sig_energy = str2double(sprintf('%.7f', sum(abs(y).^2)));
q_avg_sig_power = str2double(sprintf('%.9f', mean(abs(y).^2)));
q_rms_amplitude = str2double(sprintf('%.9f', sqrt(mean(abs(y).^2))));
q_variance = str2double(sprintf('%.9f', var(y,1)));
q_eucl_norm = str2double(sprintf('%.8f', norm(y,2)));
MSE = str2double(sprintf('%.9f', MSE));

%==========================================================================
% CALCULATING MSE
%==========================================================================
% Set up variables for computing MSE for six quantization values:
q_levels = [4 8 16 32 64 128];	% All quatization levels (N)
q_mse_all = zeros(size(q_levels));	% Empty list of MSE values for all N
% q_z_all - for storing all six quantized signals
q_z_all = zeros(length_of_xm, length(q_levels)); % Empty list of quantized signals for all N
q_delta_all = zeros(size(q_levels));	% Empty list of delta values for all N

%else
%    Delta = (Xmax-Xmin)/N;
%    index = floor((xm_cons-Xmin)/Delta);
%    index = min(max(index,0),N-1); % Includes Xmax in the last cell.
%    X_Quantized = Xmin + (index+0.5)*Delta;		% Quantized signal: X_Quantized 
%end
% Compute MSE for six quantization values:
for pos = 1:length(q_levels)
	% Calculate quantized signal for current N
	current_n = q_levels(pos);
	current_delta = (Xmax-Xmin)/current_n;

	sprintf('N = %d, Delta = %.f', current_n, current_delta);
	current_index = floor((x-Xmin)/current_delta);
	current_index = min(max(current_index,0),current_n-1);

	% The quantized signal for current N being evaluated:
	current_xq = Xmin + (current_index+0.5)*current_delta;

	% Populate list of values for Quantized Signal, Delta and MSE
	q_mse_all(pos) = mean((xm_cons-current_xq).^2);
	q_z_all(:,pos) = current_xq;
	q_delta_all(pos) = current_delta;
end

%==========================================================================
% PLOT ALL QUANTIZATION LEVELS GRAPH
%==========================================================================
%pos_n = length(q_levels);
for cmp_figure_pos = 1:2
	fig_title = sprintf('DTMF Quantization Comparison %d', cmp_figure_pos);
	figure('Name', fig_title);
	% Plot the original signal first:
	subplot(4,1,1); plot(time,xm_cons); grid on;
	xlabel('Time (s)'); ylabel('Amplitude'); title('Original signal');

	% Plot 3 out of 6 quantized signals per figure
	for cmp_plot_pos = 1:3
		if mod(cmp_figure_pos,2) == 1
			n_position = cmp_plot_pos;
		else
			n_position = cmp_plot_pos+3;
		end
		plot_title = sprintf('Mid-rise, N = %d (%d bits), \x0394_N = %.6g', q_levels(n_position), round(log2(q_levels(n_position))), q_delta_all(n_position));
		subplot(4,1,cmp_plot_pos+1); plot(time,q_z_all(:,n_position)); grid on;
		xlabel('Time (s)'); ylabel('Amplitude'); title(plot_title);
	end


end

%==========================================================================
% PLOT MSE GRAPH
%==========================================================================
% Plot: Mid-rise quantizer MSE graph
q_mse_title = sprintf('Mid-rise quantizer MSE (%d samples)', length_of_xm);
figure('Name', 'Mid-rise quantized MSE');
plot(q_levels, q_mse_all, '-o'); grid on;
xlabel('Number of quantization levels, N'); ylabel('Mean squared quantization error (MSE)'); title(q_mse_title);


fprintf('\nInput range: [%.12g, %.12g]\n',Xmin,Xmax);
fprintf('N = %d, bits/sample = %d, Delta = %.12g\n',N,round(log2(N)),Delta);
fprintf('Output range: [%.12g, %.12g]\n',min(X_Quantized),max(X_Quantized));
fprintf('MSE = %.12g\n',MSE);
if Delta > 0
    fprintf('Highest reconstruction level: %.12g\n',Xmax-Delta/2);
end

%==========================================================================
% CSV FILE
%==========================================================================
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

% PREPARING MATRIX ROWS:

% First row
r0_title = ["Signal" "SampleCount" "SamplingFrequency_Hz" "QuantizationLevels" "BitsPerSample" "QuantizationStep" "ArithmeticMean" "TotalEnergy" "AveragePower" "RMS" "Variance" "EuclideanNorm" "MSE"];

% Second row
r1_original = ["Original" length_of_xm Fs orig_N orig_bits_per_sample orig_delta arithmetic_mean total_sig_energy avg_sig_power rms_amplitude variance eucl_norm orig_mse];

% Third row
r2_modified = ["MidRiseQuantized" length_of_xm Fs N q_bits_per_sample Delta q_arithmetic_mean q_total_sig_energy q_avg_sig_power q_rms_amplitude q_variance q_eucl_norm MSE];

% Put all rows in single matrix
csv_output = [r0_title; r1_original; r2_modified];
% Create CSV file
writematrix(csv_output, 'inel4301_s086_gp05_sn05_tk05_N128_mca03.csv');

Ask = input('Play quantized audio? 1=yes, other=no: ');
if isequal(Ask,1), sound(X_Quantized,Fs); end
