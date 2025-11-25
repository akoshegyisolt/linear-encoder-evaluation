function y = alpha_trim_filter(x, k, p)     %This function was used to test the Alpha-trim filter on the signal
    N = length(x);
    y = zeros(size(x));
    pad_size = floor(k/2);
    x_pad = [zeros(pad_size,1); x; zeros(pad_size,1)]; % Simple zero padding

    for n = 1:N
        window = x_pad(n:n+k-1);
        y(n) = trimmean(window, p);
    end
end

clear
clc
close all


%% Parameters
fs      = 100;          % Sampling frequency in Hz
lambda  = 532e-9;       % Laser wavelength in meters (532 nm)

% Stepper + screw data
pitch_mm        = 0.5;      % mm/rotation of micrometer screw
pitch           = pitch_mm * 1e-3;          % m/rotation
steps_per_rev   = 200;                      % full steps / rotation
microsteps      = 64;                       % microsteps / full step

%% File dialog to load CSV
[filename, pathname] = uigetfile( ...
    {'*.csv','CSV files (*.csv)'; '*.*','All Files (*.*)'}, ...
    'Select interferometer data CSV file');

if isequal(filename,0)
    error('No file selected. Analysis aborted.');
end

fullpath = fullfile(pathname, filename);
data = readmatrix(fullpath);

%% Extract columns
id   = data(:, 1);
sigA = data(:, 2);
sigB = data(:, 3);
adc  = data(:, 4);

sigA = sigA(2:end);
sigB = sigB(2:end);
adc  = adc(2:end);

N = length(adc);
t = (0:N-1).' / fs;  % time vector (s)

%% Interferometer physics
disp_per_period      = lambda / 2;           % m per fringe
disp_per_period_nm   = disp_per_period * 1e9;
disp_per_period_um   = disp_per_period * 1e6;

%% Stepper/micrometer resolution
microstep_size        = pitch / (steps_per_rev * microsteps);   % m/microstep
microstep_size_nm     = microstep_size * 1e9;
microsteps_per_fringe = disp_per_period / microstep_size;

%% --- Set your actual motor step frequency here ---
f_stepper = 100;   % [microsteps / s]


%% Expected fringe frequency from mechanics
f_ideal = f_stepper / microsteps_per_fringe;    % Hz

%% Prepare real (measured) signal and compute FFT
adc_detr = adc - mean(adc);
w        = hann(N);
adc_win  = adc_detr .* w;

Y = fft(adc_win);
f = (0:N-1).' * (fs / N);

halfN   = floor(N/2);
f_pos   = f(1:halfN);
Y_pos   = Y(1:halfN);
amp_pos = abs(Y_pos);

[~, idx_peak] = max(amp_pos(2:end));
idx_peak      = idx_peak + 1;
f0_meas       = f_pos(idx_peak);

%% Generate ideal sinusoidal interferometer signal
% Using ideal fringe frequency from mechanics:
phi_ideal  = 2*pi*f_ideal .* t;   % phase = 2*pi*f*t
ideal_adc  = sin(phi_ideal);

% scale roughly to measured std (optional)
ideal_adc  = ideal_adc * std(adc_detr);

%% FFT of ideal signal
ideal_adc_win = (ideal_adc - mean(ideal_adc)) .* w;
Y_ideal       = fft(ideal_adc_win);
Y_ideal_pos   = Y_ideal(1:halfN);
amp_ideal_pos = abs(Y_ideal_pos);

[~, idx_peak_ideal] = max(amp_ideal_pos(2:end));
idx_peak_ideal      = idx_peak_ideal + 1;
f0_ideal_fft        = f_pos(idx_peak_ideal);

%% === Displacement vs time from instantaneous phase (handles two close peaks) ===
% Band-pass around the fringe frequency to get a narrowband signal
bw  = 0.8;                         % half-width of band in Hz (adjust if needed)
f1  = max(0.01, f_ideal - bw);
f2  = min(fs/2 - 0.01, f_ideal + bw);

adc_bp = bandpass(adc_detr, [f1 f2], fs);      

% Analytic signal and phase
z   = hilbert(adc_bp);
phi = unwrap(angle(z));                        % phase in radians

% Displacement vs time: x(t) = (phi/2pi)*disp_per_period
disp_t      = (phi - phi(1)) / (2*pi) * disp_per_period;  % m, start at 0
disp_t_um   = disp_t * 1e6;                               % µm
total_disp_phase_um = disp_t_um(end) - disp_t_um(1);

%% Display results
fprintf('Laser wavelength:                 %.0f nm\n', lambda * 1e9);
fprintf('Displacement per fringe period:   %.3f nm (%.3f µm)\n', ...
        disp_per_period_nm, disp_per_period_um);

fprintf('\nStepper/micrometer mechanics:\n');
fprintf('  Micrometer pitch:               %.3f mm/rev\n', pitch_mm);
fprintf('  Microstep size:                 %.3f nm/microstep\n', microstep_size_nm);
fprintf('  Microsteps per fringe period:   %.3f\n', microsteps_per_fringe);
fprintf('  Stepper rate (microsteps/s):    %.3f\n', f_stepper);
fprintf('  Expected fringe frequency:      %.5f Hz\n', f_ideal);

fprintf('\nFFT comparison:\n');
fprintf('  Measured dominant frequency:    %.5f Hz\n', f0_meas);
fprintf('  Ideal FFT peak frequency:       %.5f Hz\n', f0_ideal_fft);
fprintf('  Frequency difference:           %.5f Hz\n', f0_meas - f_ideal);

fprintf('\nEstimated total displacement (avg from FFT):\n');
fprintf('  Total displacement:     %.3f µm\n', total_disp_phase_um);

%% Plots
figure;

subplot(3,1,1);
plot(t, adc);
xlabel('Time [s]');
ylabel('ADC (measured)');
title('Measured interferometer signal');

subplot(3,1,2);
plot(t, ideal_adc);
xlabel('Time [s]');
ylabel('ADC (ideal)');
title('Ideal sinusoidal interferometer signal (from mechanics)');

subplot(3,1,3);
plot(f_pos, amp_pos, 'DisplayName','Measured');
hold on;
plot(f_pos, amp_ideal_pos, '--', 'DisplayName','Ideal');
hold off;
xlim([0 fs/2]);
xlabel('Frequency [Hz]');
ylabel('Amplitude');
title('Amplitude spectra: measured vs ideal');
l=legend('show');
l.Location="southeast";
grid on;


%% === Clean signals calculation ===
% Rolling average window (e.g., 10 samples). 
% Adjust 'window_size' based on how much smoothing you want.
window_size = 100; 
sigA_clean = movmedian(sigA, window_size);
sigB_clean = movmedian(sigB, window_size);


%% Displacement vs time with A and B on the same time axis
figure;
yyaxis left
plot(t, disp_t_um, 'k', 'LineWidth', 1.2, 'DisplayName', 'Displacement');
ylabel('Displacement [\mum]');

yyaxis right
% 1. Raw Signals (Magenta and Cyan)
plot(t, sigA, '-', 'Color', 'm', 'LineWidth', 0.5, 'DisplayName', 'Sig A (Raw)'); hold on;
plot(t, sigB, '-', 'Color', 'c', 'LineWidth', 0.5, 'DisplayName', 'Sig B (Raw)');

% 2. Clean Rolling Avg Signals (Red and Blue)
plot(t, sigA_clean, '-', 'Color', 'r', 'LineWidth', 1.5, 'DisplayName', 'Sig A (Clean)');
plot(t, sigB_clean, '-', 'Color', 'b', 'LineWidth', 1.5, 'DisplayName', 'Sig B (Clean)');

ylim('auto');
hold off;

ylabel('Signal A, B [ADC]');
xlabel('Time [s]');
title('Displacement and signals A/B vs time');
legend('show', 'Location', 'southeast');
grid on;

%% === Fit ideal sinusoids to sigA_clean and sigB_clean (spatial period 0.75 mm) ===

% Displacement from interferometer in mm
disp_mm = disp_t_um / 1000;   % µm -> mm

L_mm = 0.75;                  % encoder wavelength (mm)
k_spatial = 2*pi / L_mm;      % spatial angular frequency [rad/mm]

% Design matrix: sin(k*d), cos(k*d), and offset, as function of displacement
Xsp = [sin(k_spatial * disp_mm), ...
       cos(k_spatial * disp_mm), ...
       ones(N,1)];

% Least-squares fit for sigA_clean
paramsA   = Xsp \ sigA_clean;           % [A_sin; A_cos; offset]
sigA_fit  = Xsp * paramsA;              % fitted ideal sinusoid for A

% Least-squares fit for sigB_clean
paramsB   = Xsp \ sigB_clean;
sigB_fit  = Xsp * paramsB;              % fitted ideal sinusoid for B

%% Scale fitted signals to [300, 3796]
new_min = 300;
new_max = 3796;

% For sigA_fit
old_min_A = min(sigA_fit);
old_max_A = max(sigA_fit);
sigA_fit_scaled = (sigA_fit - old_min_A) * (new_max - new_min) / (old_max_A - old_min_A) + new_min;

% For sigB_fit
old_min_B = min(sigB_fit);
old_max_B = max(sigB_fit);
sigB_fit_scaled = (sigB_fit - old_min_B) * (new_max - new_min) / (old_max_B - old_min_B) + new_min;

sigA_fit_int = round(sigA_fit_scaled);
sigB_fit_int = round(sigB_fit_scaled);


%% Plot: original vs fitted (vs time, same time base)
figure;

subplot(2,1,1);
plot(t, sigA_clean, 'b', 'LineWidth', 1, 'DisplayName', 'sigA\_clean'); hold on;
plot(t, sigA_fit_scaled,   'r', 'LineWidth', 1.5, 'DisplayName', 'sigA\_fit (ideal)');
hold off;
xlabel('Time [s]');
ylabel('Signal A [ADC]');
title('Signal A: measured vs ideal sinusoid (0.75 mm spatial period)');
legend('show', 'Location', 'best');
grid on;

subplot(2,1,2);
plot(t, sigB_clean, 'b', 'LineWidth', 1, 'DisplayName', 'sigB\_clean'); hold on;
plot(t, sigB_fit_scaled,   'r', 'LineWidth', 1.5, 'DisplayName', 'sigB\_fit (ideal)');
hold off;
xlabel('Time [s]');
ylabel('Signal B [ADC]');
title('Signal B: measured vs ideal sinusoid (0.75 mm spatial period)');
legend('show', 'Location', 'best');
grid on;

%Delete the following rows in order to perform calculations on the original data
sigA_clean = sigA_fit_int;
sigB_clean = sigB_fit_int;

%% Zhao et al.
Us = sigB_clean-mean(sigB_clean);
Us = Us / max(abs(Us));
Uc = sigA_clean-mean(sigA_clean);
Uc = Uc / max(abs(Uc));

%Linearization
term1 = abs(Us);
term2 = abs(Uc);
term3 = (1/sqrt(2)) * abs(Us + Uc);
term4 = (1/sqrt(2)) * abs(Us - Uc);
L1 = min(term1, min(term2, min(term3, term4)));
L2 = term1-term2;
L3 = (L2)./(term1+term2);


V1 = L3;
Vtri = 4/pi*abs(asin(Us))-1;
E1 = V1-Vtri;

C = Vtri./V1;
V2 = C.*V1;

g = term1./term2;

Sg1 = Us>=0;
Sg2 = Uc>=0;
part = ~Sg1;
sele = xor(Sg1, Sg2);
thetaC = (pi/2) * ( ((-1).^sele .* V2 + 1)/2 + 2*part + sele );
displacement = unwrap(thetaC);
displacement = (displacement-displacement(1)) / (2*pi) * 750;
err = disp_t_um - displacement;

%% Error metrics
rmse_val = rms(err);
mae_val = mean(abs(err));
R_matrix = corrcoef(disp_t_um, displacement); 
r_sq_val = R_matrix(1,2).^2;

fprintf('RMSE: %.4f\n', rmse_val);
fprintf('MAE:  %.4f\n', mae_val);
fprintf('R^2:  %.4f\n', r_sq_val);


%% Plots
figure;
yyaxis left
plot(t, disp_t_um, 'b', 'LineWidth', 2, 'DisplayName', 'Displacement from interferometer');
hold on;
plot(t, displacement, '-r', 'LineWidth', 2, 'DisplayName', 'Displacement from ADC')
hold off
ylabel('Displacement [\mum]');

yyaxis right
plot(t, err, 'Color', [0.17, 0.63, 0.17])

ylim('auto');
hold off;

ylabel('Error of displacement calculated from raw sensor data [\mum]', 'Color', [0.17, 0.63, 0.17], 'FontWeight', 'bold');
ax = gca; % Get the handle for the current axes
ax.YColor = [0.17, 0.63, 0.17]; % Change everything related to Y to green

grid on

legend('show', 'Location', 'southeast');


%% Reactivity


tol = 1e-9;   % µm, treat differences smaller than this as "no change"

delta_out = diff(displacement);        % change in Zhao output per sample
same_out  = abs(delta_out) < tol;          % plateau where output stays constant

plateau_spans  = [];
plateauStarts  = [];
plateauEnds    = [];

start_idx = 1;
for i = 2:length(displacement)
    if ~same_out(i-1)
        % [start_idx .. i-1] is one plateau
        span = disp_t_um(i-1) - disp_t_um(start_idx);
        plateau_spans(end+1)  = span;      %#ok<SAGROW>
        plateauStarts(end+1)  = start_idx; %#ok<SAGROW>
        plateauEnds(end+1)    = i-1;       %#ok<SAGROW>
        start_idx = i;
    end
end
% last plateau:
span = disp_t_um(end) - disp_t_um(start_idx);
plateau_spans(end+1) = span;
plateauStarts(end+1) = start_idx;
plateauEnds(end+1)   = length(displacement);

[max_plateau, idx_max] = max(plateau_spans);      % worst-case flat region

% --- Best resolution point: smallest displacement step that changed the output ---
delta_disp_true = diff(disp_t_um);               % true displacement step between samples
mask_change = ~same_out;                         % where output actually changes
[best_step, idx_rel] = min(delta_disp_true(mask_change));  % smallest Δx with Δoutput ≠ 0
change_indices   = find(mask_change);
idx_best_delta   = change_indices(idx_rel);      % index in delta_out / delta_disp_true
idx_best_sample  = idx_best_delta + 1;           % sample index in displacement / disp_t_um

fprintf('Worst-case plateau of Zhao output:         %.4f µm\n', max_plateau);
fprintf('=> Any displacement change > %.4f µm is guaranteed to change the output.\n', max_plateau);
fprintf('Best resolution step (smallest Δx with Δoutput≠0): %.4f µm\n', best_step);

idx_worst_start = plateauStarts(idx_max);
idx_worst_end   = plateauEnds(idx_max);
idx_worst       = idx_worst_start:idx_worst_end;

% for delta_out indexing (plateau s..e  ⇒ deltas s..e-1)
idx_worst_delta = idx_worst_start:(idx_worst_end-1);

figure;

% --- 1) Integer ideal signals ---
subplot(3,1,1);
plot(disp_t_um, sigA_fit_int, '.-', 'DisplayName','A int'); hold on;
plot(disp_t_um, sigB_fit_int, '.-', 'DisplayName','B int');
xlabel('Displacement [\mum]'); ylabel('ADC counts');
title('Integer ideal encoder signals');
legend; grid on;

% --- 2) Zhao output vs true displacement ---
subplot(3,1,2);
plot(disp_t_um, displacement, 'LineWidth',1.2, ...
     'DisplayName','Zhao displacement'); hold on;

% Worst plateau: thick red segment
plot(disp_t_um(idx_worst), displacement(idx_worst), ...
     'r', 'LineWidth', 2, 'DisplayName','Worst plateau');

% Best resolution point: green star
plot(disp_t_um(idx_best_sample), displacement(idx_best_sample), ...
     'g*', 'MarkerSize',10, 'DisplayName','Best resolution point');

xlabel('True displacement [\mum]');
ylabel('Zhao displacement [\mum]');
title('Zhao output vs true displacement (integer signals)');
legend('Location','best');
grid on;

% --- 3) Δ-output per step ---
subplot(3,1,3);
plot(disp_t_um(2:end), delta_out, 'LineWidth',1, ...
     'DisplayName','\Delta output per sample'); hold on;

% Worst plateau steps: red dots at Δoutput≈0
plot(disp_t_um(idx_worst_delta+1), delta_out(idx_worst_delta), ...
     'r.', 'MarkerSize',12, 'DisplayName','Worst plateau steps');

% Best resolution step: green star
plot(disp_t_um(idx_best_delta+1), delta_out(idx_best_delta), ...
     'g*', 'MarkerSize',10, 'DisplayName','Best resolution step');

xlabel('True displacement [\mum]');
ylabel('\Delta Zhao output between samples [\mum]');
title('Output change per sample');
legend('Location','best');
grid on;
hold off;
