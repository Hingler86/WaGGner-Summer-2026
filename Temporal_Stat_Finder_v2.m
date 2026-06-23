%% Temporal Response Identification: Command -> RPM -> Velocity
%
% PROJECT CONTEXT
%   Characterizes the temporal (frequency) response of a single fan block
%   in the WaGGner fan-array wind generation wall. A sinusoidal PWM
%   command is applied to the fan at a sweep of forcing periods; the
%   resulting fan speed (RPM) and downstream flow velocity (measured by a
%   five-hole probe at X = 500 mm, on the jet centerline) are recorded
%   simultaneously on three independent acquisition systems with three
%   independent clocks. This script identifies, for each forcing period:
%     - the gain and phase lag of RPM relative to the command
%     - the gain and phase lag of streamwise velocity (Ux) relative to
%       the command
%   and combines those across all periods into a Bode-style frequency
%   response, from which a first-order-lag-plus-delay transfer function
%   can be identified (see fit_first_order_delay.m).
%
% INPUT FILE TRIO (one trio per forcing period T, all in the same folder
% as this script, or anywhere on the MATLAB path):
%   Commande_Scenario_Block_3_3_Sinus_<T*1000>ms.lvm
%       PWM command log. 50 data channels + a leading X_Value (time)
%       column; only one channel (column 26) carries the real command --
%       the rest are unused/zero. Filename has NO trailing ID.
%   Mesure_Scenario_Block_3_3_Sinus_<T*1000>ms.lvm
%       Fan RPM log, same 50-ish-channel layout. Column 7 is the RPM
%       channel; it idles at -1 between updates (placeholder, not a real
%       reading) and otherwise reports true RPM directly -- already
%       calibrated, no conversion needed. Filename has NO trailing ID.
%   Sinu_x500_P30-90_T<T>_t<trial>
%       Five-hole probe log (no file extension). Tab-delimited, comma
%       decimals, fixed 100 Hz sample rate. Column 13 is Ux (m/s). The
%       trial number in the filename isn't derivable from T alone, so
%       this one file is found by wildcard.
%
% WHY THREE SEPARATE TIME AXES NEED RECONCILING
%   Each file's elapsed-time column starts at 0 relative to that
%   acquisition system's OWN start trigger -- the three systems were not
%   triggered in lockstep. Each LVM file's header also records an
%   absolute wall-clock start time for that channel; reading those three
%   header timestamps and differencing them against the command file's
%   gives the offset needed to bring RPM and velocity onto the command
%   file's time axis (see Section 2, "Header timestamps & offsets").
%
% WHY ONSET/OFFSET OF THE ACTIVE SINUSOID CAN'T BE FOUND FROM SIGNAL
% VALUE ALONE
%   The command idles at PWM = 60 between tests, which is also the
%   sinusoid's own mean -- so a single sample's value (or even a window
%   of values) can't distinguish "idling at the mean" from "passing
%   through the mean mid-cycle." What does distinguish them, verified
%   directly on real data: the command logger is event/change-triggered.
%   It idles with large, irregular gaps between samples (seconds to tens
%   of seconds) and only logs densely and regularly (~0.1 s apart,
%   regardless of T) while the sinusoid is actively running. detect_window
%   (Section 5) finds the start and end of the active region purely from
%   this gap structure in the command file, then that same window (in
%   command-clock time) is applied to the RPM and velocity signals after
%   they've been synchronized onto that clock.
%
% SIGN CONVENTION
%   With each signal fit to y = offset + A*sin(phi + psi), phase lag is
%   defined as wrap180(psi_command - psi_output), which is POSITIVE when
%   the output lags the command (the physically expected case for both
%   RPM and velocity), and the corresponding delay in seconds is
%   phaselag/360 * T.

clear; clc; close all;

%% 1. Experimental Parameter Setup
all_periods = [1, 2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20, 22.5, 25, 27.5, 30, 32.5, 35];

n_skip_cycles          = 1;     % full cycles discarded after onset, to let
                                 % transients settle before fitting
nbins_phase            = 36;    % bins for the phase-averaged-response plots (10 deg/bin)
make_diagnostic_plots  = false;  % set false to suppress the per-period figures
representative_period  = 15;    % which forcing period (s) to capture for the
                                 % normalized RPM/velocity time-domain comparison
                                 % figure in Section 3.5 -- must be one of the
                                 % values in all_periods
representative_n_cycles = 2;     % how many forcing cycles to display in that
                                 % figure -- kept short so the RPM-vs-velocity
                                 % peak lag is visually clear; the full active
                                 % window (~9+ cycles) would compress the lag
                                 % too much to see at a glance

% Captured inside the batch loop when T == representative_period, used to
% build the time-domain comparison figure after the loop ends.
rep_t_rpm = []; rep_rpm = []; rep_t_vel = []; rep_vel = []; rep_T = NaN;

num_runs = length(all_periods);
summary_freq         = nan(num_runs, 1);
summary_gain_rpm     = nan(num_runs, 1);
summary_phaselag_rpm = nan(num_runs, 1);
summary_delay_rpm    = nan(num_runs, 1);
summary_gain_vel     = nan(num_runs, 1);
summary_phaselag_vel = nan(num_runs, 1);
summary_delay_vel    = nan(num_runs, 1);
summary_r2_cmd       = nan(num_runs, 1);
summary_r2_rpm       = nan(num_runs, 1);
summary_r2_vel       = nan(num_runs, 1);

fprintf('===========================================================\n');
fprintf('     STARTING SYSTEM IDENTIFICATION & BATCH PROCESSING     \n');
fprintf('===========================================================\n\n');

%% 2. Batch Processing Loop
% Repeats the full pipeline -- file discovery, clock sync, active-window
% detection, sine-regression fitting -- independently for every forcing
% period, accumulating one row of results per period into the summary_*
% arrays above.
for i = 1:num_runs
    T = all_periods(i);
    summary_freq(i) = 1 / T;

    % --- File discovery ---
    % Command/Measurement filenames are exact (no trailing ID); only the
    % probe file's trial number needs a wildcard.
    ms_val = round(T * 1000);
    cmd_filename   = find_file_matching(sprintf('Commande_Scenario_Block_3_3_Sinus_%dms.lvm', ms_val));
    mes_filename   = find_file_matching(sprintf('Mesure_Scenario_Block_3_3_Sinus_%dms.lvm', ms_val));
    probe_filename = find_file_matching(sprintf('Sinu_x500_P30-90_T%g_t*', T));

    fprintf('--- Processing Forcing Period: T = %g seconds (f = %.4f Hz) ---\n', T, summary_freq(i));

    if isempty(cmd_filename) || isempty(mes_filename) || isempty(probe_filename)
        warning('Missing data files for T = %g s. Skipping run.', T);
        continue;
    end

    % --- Header timestamps & offsets ---
    % Each LVM header has two "Time" lines: a file-level creation time,
    % and (after the first ***End_of_Header*** marker) the per-channel
    % acquisition start time that the data rows are actually relative to.
    % read_lvm_channel_start_time grabs the second one specifically.
    cmd_time_str   = read_lvm_channel_start_time(cmd_filename);
    mes_time_str   = read_lvm_channel_start_time(mes_filename);
    probe_time_str = parse_probe_start_time(probe_filename);

    t0_cmd   = datetime(cmd_time_str, 'InputFormat', 'HH:mm:ss.SSSSSS');
    t0_mes   = datetime(mes_time_str, 'InputFormat', 'HH:mm:ss.SSSSSS');
    t0_probe = datetime(probe_time_str, 'InputFormat', 'HH:mm:ss.SSSSSS');

    % Offsets map the Measurement/Probe clocks onto the Command clock:
    % t_synced = t_raw + offset puts every signal on the same time axis.
    offset_mes   = seconds(t0_mes - t0_cmd);
    offset_probe = seconds(t0_probe - t0_cmd);

    % --- Data extraction ---
    % Both LVM files (.lvm extension) and the probe file (no extension)
    % are tab-delimited text with comma decimals; 'FileType','text' is
    % required because detectImportOptions can't infer a file type from
    % an unrecognized/missing extension. dataStartLine skips the header
    % block specific to each file (23 for Command/Measurement, 6 for the
    % probe file -- see read_european_lvm below).
    cmd_data   = read_european_lvm(cmd_filename, 23);
    mes_data   = read_european_lvm(mes_filename, 23);
    probe_data = read_european_lvm(probe_filename, 6);

    t_cmd_raw   = cmd_data(:, 1);
    t_mes_raw   = mes_data(:, 1);
    % The probe file's own time column is a date/time string, not a
    % plain number, and isn't needed anyway: its sample rate is a
    % verified-uniform 100 Hz, so a synthetic grid is simpler and exact.
    t_probe_raw = (0:size(probe_data,1)-1)' * 0.010;

    pwm_command = cmd_data(:, 26);   % the only nonzero command channel
    rpm         = mes_data(:, 7);    % true RPM (already calibrated);
                                      % idles at -1 between sensor updates
    u_velocity  = probe_data(:, 13); % Ux (m/s)

    t_mes_synced   = t_mes_raw + offset_mes;
    t_probe_synced = t_probe_raw + offset_probe;

    % --- RPM channel: drop the -1 idle/placeholder rows ---
    % The sensor only updates roughly every 5 samples of the measurement
    % file's own clock; -1 marks rows where no new reading was available.
    % No other transformation is applied: this channel is already true,
    % instantaneous RPM, not a cumulative pulse count, so it must NOT be
    % differentiated.
    valid_mask  = rpm > -0.5;
    t_rpm_valid = t_mes_synced(valid_mask);
    rpm_valid   = rpm(valid_mask);

    if numel(rpm_valid) < 3
        warning('Not enough valid RPM samples for T = %g s. Skipping run.', T);
        continue;
    end

    % --- Active window detection ---
    % Finds where the command's sinusoid actually starts and stops
    % (as opposed to the surrounding idle ramp-up/ramp-down), using the
    % command file's own sampling-gap structure -- see detect_window for
    % why this can't be done from signal value alone.
    [t_onset, t_offset] = detect_window(t_cmd_raw);
    if isnan(t_onset) || isnan(t_offset)
        warning('Could not detect active sinusoid window for T = %g s. Skipping run.', T);
        continue;
    end

    t_analysis_start = t_onset + n_skip_cycles * T;  % skip transient cycles
    t_analysis_end   = t_offset;

    if t_analysis_end - t_analysis_start < T
        warning('Less than one full cycle available after trimming for T = %g s. Skipping run.', T);
        continue;
    end

    % --- Per-signal phase tagging ---
    % Each signal keeps its own native (synchronized) timestamps -- no
    % resampling/interpolation onto a common grid -- and is converted
    % to a phase angle within the forcing cycle via mod(t,T)/T*2*pi.
    mask_cmd = t_cmd_raw >= t_analysis_start & t_cmd_raw <= t_analysis_end;
    mask_rpm = t_rpm_valid >= t_analysis_start & t_rpm_valid <= t_analysis_end;
    mask_vel = t_probe_synced >= t_analysis_start & t_probe_synced <= t_analysis_end;

    phi_cmd = mod(t_cmd_raw(mask_cmd)      - t_analysis_start, T) / T * 2*pi;
    phi_rpm = mod(t_rpm_valid(mask_rpm)    - t_analysis_start, T) / T * 2*pi;
    phi_vel = mod(t_probe_synced(mask_vel) - t_analysis_start, T) / T * 2*pi;

    y_cmd = pwm_command(mask_cmd);
    y_rpm = rpm_valid(mask_rpm);
    y_vel = u_velocity(mask_vel);

    % --- Capture for the time-domain RPM/velocity comparison figure ---
    % Keeps the windowed, synced time + signal arrays for whichever period
    % matches representative_period, since the loop overwrites these
    % variables on every iteration and they would otherwise be lost by the
    % time the loop finishes. Limited to representative_n_cycles so the
    % resulting figure shows the RPM/velocity lag clearly rather than
    % compressing many cycles onto one wide axis.
    if T == representative_period
        t_rep_end = t_analysis_start + representative_n_cycles * T;
        rep_mask_rpm = mask_rpm & (t_rpm_valid <= t_rep_end);
        rep_mask_vel = mask_vel & (t_probe_synced <= t_rep_end);
        rep_t_rpm = t_rpm_valid(rep_mask_rpm);
        rep_rpm   = rpm_valid(rep_mask_rpm);
        rep_t_vel = t_probe_synced(rep_mask_vel);
        rep_vel   = u_velocity(rep_mask_vel);
        rep_T     = T;
    end

    if numel(y_cmd) < 3 || numel(y_rpm) < 3 || numel(y_vel) < 3
        warning('Not enough samples inside the active window for T = %g s. Skipping run.', T);
        continue;
    end

    % --- Sine regression fits: y = offset + A*sin(phi + psi) ---
    [A_cmd, psi_cmd, off_cmd] = fit_sine_regression(phi_cmd, y_cmd);
    [A_rpm, psi_rpm, off_rpm] = fit_sine_regression(phi_rpm, y_rpm);
    [A_vel, psi_vel, off_vel] = fit_sine_regression(phi_vel, y_vel);

    % --- Fit quality (R^2) for each signal's sine regression ---
    % Quantifies how much of each signal's variance the fitted sinusoid
    % actually explains, independent of gain/phase/delay. This matters
    % most at the fastest forcing periods: a low R^2 there means the
    % corresponding phase/delay number is statistically unreliable (the
    % fit is largely chasing noise, not a real phase-locked oscillation),
    % even though the script still reports a number for it regardless.
    % See Section 9.6.3 of the report for how this is used to judge which
    % periods' results should be trusted.
    r2_cmd = sine_fit_r2(phi_cmd, y_cmd, A_cmd, psi_cmd, off_cmd);
    r2_rpm = sine_fit_r2(phi_rpm, y_rpm, A_rpm, psi_rpm, off_rpm);
    r2_vel = sine_fit_r2(phi_vel, y_vel, A_vel, psi_vel, off_vel);

    summary_r2_cmd(i) = r2_cmd;
    summary_r2_rpm(i) = r2_rpm;
    summary_r2_vel(i) = r2_vel;

    % --- Gain / Phase lag / Delay ---
    % Phase lag is POSITIVE when the output lags the command (see header
    % docstring for the sign convention).
    gain_rpm     = A_rpm / A_cmd;
    phaselag_rpm = wrap180(rad2deg(psi_cmd - psi_rpm));
    delay_rpm    = phaselag_rpm / 360 * T;

    gain_vel     = A_vel / A_cmd;
    phaselag_vel = wrap180(rad2deg(psi_cmd - psi_vel));
    delay_vel    = phaselag_vel / 360 * T;

    summary_gain_rpm(i)     = gain_rpm;
    summary_phaselag_rpm(i) = phaselag_rpm;
    summary_delay_rpm(i)    = delay_rpm;
    summary_gain_vel(i)     = gain_vel;
    summary_phaselag_vel(i) = phaselag_vel;
    summary_delay_vel(i)    = delay_vel;

    fprintf('  >> RPM (min^-1)       Gain: %8.3f   Phase Lag: %7.2f deg   Delay: %.4f s   R^2: %.4f\n', gain_rpm, phaselag_rpm, delay_rpm, r2_rpm);
    fprintf('  >> Velocity (m/s)     Gain: %8.3f   Phase Lag: %7.2f deg   Delay: %.4f s   R^2: %.4f\n\n', gain_vel, phaselag_vel, delay_vel, r2_vel);

    % --- Diagnostic plot: phase-averaged response + fitted sine overlay ---
    % Bins all in-window samples by phase (regardless of which cycle they
    % came from) so multiple forcing cycles are effectively averaged
    % together into one representative period, with error bars showing
    % the cycle-to-cycle scatter at each phase bin.
    if make_diagnostic_plots
        phi_dense = linspace(0, 2*pi, 200);

        figure('Name', sprintf('Phase-averaged response, T = %g s', T), 'Color', [1 1 1]);

        subplot(3,1,1);
        [bc, bm, bs] = phase_bin_average(phi_cmd, y_cmd, nbins_phase);
        errorbar(bc, bm, bs, 'o', 'Color', [0 0.4470 0.7410]); hold on;
        plot(rad2deg(phi_dense), off_cmd + A_cmd*sin(phi_dense+psi_cmd), 'k-', 'LineWidth', 1.5);
        ylabel('Command (%)'); title(sprintf('T = %g s -- Phase-Averaged Cycle', T)); grid on;

        subplot(3,1,2);
        [bc, bm, bs] = phase_bin_average(phi_rpm, y_rpm, nbins_phase);
        errorbar(bc, bm, bs, 'o', 'Color', [0.4660 0.6740 0.1880]); hold on;
        plot(rad2deg(phi_dense), off_rpm + A_rpm*sin(phi_dense+psi_rpm), 'k-', 'LineWidth', 1.5);
        ylabel('RPM (min^{-1})'); grid on;

        subplot(3,1,3);
        [bc, bm, bs] = phase_bin_average(phi_vel, y_vel, nbins_phase);
        errorbar(bc, bm, bs, 'o', 'Color', [0.8500 0.3250 0.0980]); hold on;
        plot(rad2deg(phi_dense), off_vel + A_vel*sin(phi_dense+psi_vel), 'k-', 'LineWidth', 1.5);
        ylabel('U_x (m/s)'); xlabel('Phase (deg)'); grid on;
    end
end

%% 3. Final Summary Table
fprintf('=========================================================================================\n');
fprintf('                          FINAL CHARACTERIZATION RESULTS MATRIX                          \n');
fprintf('=========================================================================================\n');
fprintf('Period(s)  Freq(Hz)  GainRPM   PhaseRPM(deg)  DelayRPM(s)   GainVel   PhaseVel(deg)  DelayVel(s)\n');
fprintf('-----------------------------------------------------------------------------------------\n');
for i = 1:num_runs
    fprintf('%8.1f  %8.4f  %8.3f  %12.2f  %11.4f   %8.3f  %12.2f  %10.4f\n', ...
        all_periods(i), summary_freq(i), summary_gain_rpm(i), summary_phaselag_rpm(i), summary_delay_rpm(i), ...
        summary_gain_vel(i), summary_phaselag_vel(i), summary_delay_vel(i));
end
fprintf('=========================================================================================\n');

%% 3.1 Sine-Fit Quality (R^2) Table
% R^2 of the fitted sinusoid against the raw (unbinned) samples in the
% active window, for each of the three signals independently. Low R^2 at
% a given period means that period's gain/phase/delay numbers above are
% statistically unreliable -- the fit is mostly tracking noise rather
% than a genuine phase-locked oscillation -- which matters most for
% judging which points to trust in the frequency-response plots and
% transfer function fit that follow. See Section 9.6.3 of the report.
fprintf('\n');
fprintf('=========================================================================\n');
fprintf('                     SINE-FIT QUALITY (R^2) BY PERIOD                    \n');
fprintf('=========================================================================\n');
fprintf('Period(s)  Freq(Hz)  R2_Command  R2_RPM    R2_Velocity\n');
fprintf('-------------------------------------------------------------------------\n');
for i = 1:num_runs
    fprintf('%8.1f  %8.4f  %9.4f  %7.4f  %10.4f\n', ...
        all_periods(i), summary_freq(i), summary_r2_cmd(i), summary_r2_rpm(i), summary_r2_vel(i));
end
fprintf('=========================================================================\n');

%% 3.5 Normalized RPM/Velocity Time-Domain Comparison (representative period)
% Plots the RPM and velocity signals together, over real elapsed time
% (not wrapped into a single cycle), for whichever period was captured
% above as representative_period. Each signal is independently min-max
% normalized to [0, 1] so the two -- which live in completely different
% units and magnitudes (raw RPM counts vs. m/s) -- can be compared
% directly on one set of axes. This makes the time LAG between the RPM
% peak and the later velocity peak visible directly as a horizontal
% offset, which is a more intuitive way to see the delay than reading it
% off a phase number alone.
if ~isempty(rep_rpm) && ~isempty(rep_vel)
    rpm_norm = (rep_rpm - min(rep_rpm)) / (max(rep_rpm) - min(rep_rpm));
    vel_norm = (rep_vel - min(rep_vel)) / (max(rep_vel) - min(rep_vel));

    % Time axis zeroed to the start of this period's active window, purely
    % for a cleaner x-axis -- does not change any relative timing between
    % the two signals, since both were already on the same synced clock.
    t0_rep = min([rep_t_rpm; rep_t_vel]);

    figure('Name', sprintf('Normalized RPM/Velocity Comparison, T = %g s', rep_T), 'Color', [1 1 1]);
    plot(rep_t_rpm - t0_rep, rpm_norm, 'o-', 'Color', [0.4660 0.6740 0.1880], ...
        'MarkerSize', 3, 'LineWidth', 1); hold on;
    plot(rep_t_vel - t0_rep, vel_norm, 's-', 'Color', [0.8500 0.3250 0.0980], ...
        'MarkerSize', 3, 'LineWidth', 1);
    xlabel('Time (s)'); ylabel('Normalized signal (–)');
    title(sprintf('Normalized RPM and Velocity Time Series, T = %g s', rep_T));
    legend('RPM (normalized)', 'Velocity (normalized)', 'Location', 'best');
    grid on;
else
    warning(['representative_period (%g s) was not found among the ' ...
             'periods actually processed -- skipping the time-domain ' ...
             'comparison figure.'], representative_period);
end

%% 4. Bode-Style Gain/Phase vs Frequency Plots (wrapped phase, as fitted)
% Phase lag here is wrap180'd into (-180, 180], matching exactly what's
% in the summary table/arrays above. At the fastest forcing periods this
% can visually "cliff" even though the true underlying lag keeps growing
% smoothly past 180 deg -- see Section 4.1 for the unwrapped version,
% which is the one to use for reporting.
%
% RPM and velocity are plotted as SEPARATE figures rather than overlaid
% on shared axes: RPM gain (~250 min^-1 per % command) and velocity gain
% (~0.17 m/s per % command) differ by three orders of magnitude, so a
% shared y-axis makes the velocity curve visually collapse to a flat
% line near zero even though it has its own meaningful shape.
figure('Name', 'Frequency Response - RPM (Wrapped Phase)', 'Color', [1 1 1], 'Position', [100 100 800 600]);

subplot(2,1,1);
semilogx(summary_freq, summary_gain_rpm, 'o-', 'Color', [0.4660 0.6740 0.1880], 'LineWidth', 1.5);
ylabel('Gain (min^{-1} per % command)');
title('RPM Gain vs. Forcing Frequency'); grid on;

subplot(2,1,2);
semilogx(summary_freq, summary_phaselag_rpm, 'o-', 'Color', [0.4660 0.6740 0.1880], 'LineWidth', 1.5);
ylabel('Phase Lag (deg)'); xlabel('Frequency (Hz)');
title('RPM Phase Lag vs. Forcing Frequency (Wrapped)'); grid on;

figure('Name', 'Frequency Response - Velocity (Wrapped Phase)', 'Color', [1 1 1], 'Position', [950 100 800 600]);

subplot(2,1,1);
semilogx(summary_freq, summary_gain_vel, 's-', 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 1.5);
ylabel('Gain (m/s per % command)');
title('Velocity Gain vs. Forcing Frequency'); grid on;

subplot(2,1,2);
semilogx(summary_freq, summary_phaselag_vel, 's-', 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 1.5);
ylabel('Phase Lag (deg)'); xlabel('Frequency (Hz)');
title('Velocity Phase Lag vs. Forcing Frequency (Wrapped)'); grid on;

%% 4.1 Bode Plot, Unwrapped Phase (canonical version for reporting)
% Reconstructs the true, smoothly-growing phase trend by unwrapping in
% ascending-frequency order -- matching exactly what fit_first_order_delay.m
% does internally, so this figure and the fitted model parameters
% reflect the same underlying curve. Use THIS figure as the primary
% gain/phase-vs-frequency figure for reporting, not the wrapped one above.
% As in Section 4, RPM and velocity get separate figures to avoid the
% gain scale mismatch described above.
[freq_sorted, sort_idx] = sort(summary_freq);

phase_rpm_unwrapped = rad2deg(unwrap(deg2rad(summary_phaselag_rpm(sort_idx))));
phase_vel_unwrapped = rad2deg(unwrap(deg2rad(summary_phaselag_vel(sort_idx))));
gain_rpm_sorted = summary_gain_rpm(sort_idx);
gain_vel_sorted = summary_gain_vel(sort_idx);

figure('Name', 'Frequency Response - RPM (Unwrapped Phase)', 'Color', [1 1 1], 'Position', [100 750 800 600]);

subplot(2,1,1);
semilogx(freq_sorted, gain_rpm_sorted, 'o-', 'Color', [0.4660 0.6740 0.1880], 'LineWidth', 1.5);
ylabel('Gain (min^{-1} per % command)');
title('RPM Gain vs. Forcing Frequency'); grid on;

subplot(2,1,2);
semilogx(freq_sorted, phase_rpm_unwrapped, 'o-', 'Color', [0.4660 0.6740 0.1880], 'LineWidth', 1.5);
ylabel('Phase Lag (deg, unwrapped)'); xlabel('Frequency (Hz)');
title('RPM Phase Lag vs. Forcing Frequency (Unwrapped)'); grid on;

figure('Name', 'Frequency Response - Velocity (Unwrapped Phase)', 'Color', [1 1 1], 'Position', [950 750 800 600]);

subplot(2,1,1);
semilogx(freq_sorted, gain_vel_sorted, 's-', 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 1.5);
ylabel('Gain (m/s per % command)');
title('Velocity Gain vs. Forcing Frequency'); grid on;

subplot(2,1,2);
semilogx(freq_sorted, phase_vel_unwrapped, 's-', 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 1.5);
ylabel('Phase Lag (deg, unwrapped)'); xlabel('Frequency (Hz)');
title('Velocity Phase Lag vs. Forcing Frequency (Unwrapped)'); grid on;

%% 4.5 First-Order + Delay Transfer Function Fit (Command -> Velocity)
% Fits H(s) = K*exp(-s*Td) / (1+tau*s) to the velocity gain/phase data
% above. Requires fit_first_order_delay.m to be on the MATLAB path (same
% folder as this script is sufficient) and the Optimization Toolbox
% (lsqnonlin). The fastest period (T = 1 s) is excluded by default, since
% its underlying sine fit has near-zero R^2 -- the flow physically can't
% track a 1 s oscillation, so its phase number is noise-dominated and
% would skew the fit. Remove the 4th/5th arguments below to include all
% periods instead.
if exist('fit_first_order_delay', 'file')
    fprintf('\n');
    [K_vel, tau_vel, Td_vel, gof_vel] = fit_first_order_delay( ...
        summary_freq, summary_gain_vel, summary_phaselag_vel, 1, all_periods);
else
    warning(['fit_first_order_delay.m not found on the path -- skipping ' ...
             'transfer function fit. Place it in the same folder as this script.']);
end


%% =========================================================================
%% 5. CORE ANALYSIS FUNCTIONS
%% =========================================================================

function fname = find_file_matching(pattern)
% FIND_FILE_MATCHING Return the first filename matching a wildcard
% pattern, or '' if none found. Used for the probe filename, whose trial
% number isn't derivable from the forcing period alone.
    listing = dir(pattern);
    if isempty(listing)
        fname = '';
        return;
    end
    if numel(listing) > 1
        warning('Multiple files match pattern "%s" -- using "%s".', pattern, listing(1).name);
    end
    fname = listing(1).name;
end

function [t_onset, t_offset] = detect_window(t)
% DETECT_WINDOW Find the start and end times of the actively-logged
% sinusoid region in the command signal, sampled at (possibly irregular)
% times `t`.
%
% Why value-based detection doesn't work here: the idle/baseline value of
% the command signal is identical to the sinusoid's own mean (both are
% the PWM center point), so a single sample's value -- or even a window
% of values -- can't reliably tell "idle" apart from "sinusoid passing
% through its center." A naive amplitude-window approach also fails for
% long periods (e.g. T = 15 s): a window wide enough to contain a
% meaningful fraction of one cycle can "see across" the entire idle gap
% and falsely flag idle samples as already active.
%
% What actually distinguishes idle from active, verified on real data:
% the command logger is event/change-triggered. It idles with large,
% irregular gaps between samples (seconds to tens of seconds); while the
% sinusoid is active it logs densely and regularly (~0.1 s apart,
% independent of T). So onset/offset are found purely from the SAMPLING
% GAP structure, not from signal values at all -- this sidesteps the
% idle-equals-mean ambiguity entirely and works the same way regardless
% of period or PWM range.
    t = t(:);
    n = numel(t);
    gap_thresh   = 1.0;  % seconds; comfortably separates dense active
                         % logging (~0.1 s/sample) from sparse idle gaps
                         % (seconds+) for every period in the 1-35 s sweep
    run_len_need = 5;    % consecutive dense steps required to call a
                         % region "active" (rejects isolated stray samples)

    if n < run_len_need + 1
        t_onset = NaN; t_offset = NaN;
        return;
    end

    dt = diff(t);
    dense = dt <= gap_thresh;   % dense(k) true => gap from sample k to k+1 is small

    % --- Onset: first dense step that begins a sustained dense run ---
    t_onset = NaN;
    for i = 1:(n-1)
        if dense(i)
            j = i; count = 1;
            while j+1 < n && dense(j+1) && count < run_len_need
                count = count + 1;
                j = j + 1;
            end
            if count >= run_len_need
                t_onset = t(i+1);
                break;
            end
        end
    end

    % --- Offset: last dense step that ends a sustained dense run ---
    t_offset = NaN;
    for i = (n-1):-1:1
        if dense(i)
            j = i; count = 1;
            while j-1 >= 1 && dense(j-1) && count < run_len_need
                count = count + 1;
                j = j - 1;
            end
            if count >= run_len_need
                t_offset = t(i+1);
                break;
            end
        end
    end
end

function [A, psi, offset] = fit_sine_regression(phi, y)
% FIT_SINE_REGRESSION Linear least-squares fit of y = offset + A*sin(phi+psi)
% via the sin/cos basis trick, so no nonlinear solver is needed: writing
% A*sin(phi+psi) = c1*sin(phi) + c2*cos(phi) makes the model linear in
% (c1, c2, offset), solvable with a single backslash.
    phi = phi(:); y = y(:);
    H = [sin(phi), cos(phi), ones(size(phi))];
    coeffs = H \ y;
    c1 = coeffs(1); c2 = coeffs(2); offset = coeffs(3);
    A   = sqrt(c1^2 + c2^2);
    psi = atan2(c2, c1);
end

function r2 = sine_fit_r2(phi, y, A, psi, offset)
% SINE_FIT_R2 Coefficient of determination (R^2) for a fit_sine_regression
% result: the fraction of y's variance explained by offset + A*sin(phi+psi),
% computed on the same raw, unbinned samples the fit itself used (not the
% phase-binned averages, which would inflate R^2 by averaging out noise).
% R^2 = 1 means the sinusoid perfectly reproduces the data; R^2 near 0
% means the fitted oscillation explains essentially none of the observed
% variation, and the corresponding gain/phase/delay for that fit should be
% treated as unreliable (see Section 9.6.3 of the report).
    phi = phi(:); y = y(:);
    y_pred = offset + A*sin(phi + psi);
    ss_res = sum((y - y_pred).^2);
    ss_tot = sum((y - mean(y)).^2);
    if ss_tot == 0
        r2 = NaN;   % constant signal, R^2 undefined
    else
        r2 = 1 - ss_res/ss_tot;
    end
end

function wrapped = wrap180(deg)
% WRAP180 Wrap an angle in degrees into (-180, 180].
    wrapped = mod(deg + 180, 360) - 180;
end

function [bin_centers_deg, bin_means, bin_stds] = phase_bin_average(phi, y, nbins)
% PHASE_BIN_AVERAGE Bin (phi, y) pairs by phase for the phase-averaged
% plot only. The Gain/Phase numbers used everywhere else come from
% fit_sine_regression on the raw, unbinned samples -- these bins are
% purely a visualization aid (and a way to show cycle-to-cycle scatter
% via the bin standard deviation), not an input to the fit itself.
    phi = mod(phi(:), 2*pi);
    edges = linspace(0, 2*pi, nbins+1);
    bin_centers_deg = rad2deg((edges(1:end-1) + edges(2:end)) / 2)';
    bin_means = nan(nbins, 1);
    bin_stds  = nan(nbins, 1);
    for k = 1:nbins
        mask = phi >= edges(k) & phi < edges(k+1);
        if any(mask)
            bin_means(k) = mean(y(mask));
            bin_stds(k)  = std(y(mask));
        end
    end
end

%% =========================================================================
%% 6. DATA INGESTION UTILITIES
%% =========================================================================

function time_str = read_lvm_channel_start_time(filename)
% READ_LVM_CHANNEL_START_TIME Extract the per-channel acquisition start
% time from an LVM header, independent of time-of-day.
%
% An LVM file has two header blocks separated by '***End_of_Header***':
%   (1) a file-level block with a single global Date/Time (when LabVIEW
%       wrote the file) -- not what we want.
%   (2) a per-channel block, repeated once per column, whose Date/Time is
%       what the X_Value = 0 column is actually relative to -- this is
%       the one needed for synchronizing against the other two systems.
% This walks the file and takes the first "Time" line that appears AFTER
% the first End_of_Header marker, which is always the per-channel one
% regardless of what hour the test happened to run in (an earlier
% approach matched a literal hour string, which only worked by
% coincidence and broke for sessions starting in a different hour).
    fid = fopen(filename, 'r');
    if fid == -1, error('Could not open file: %s', filename); end
    seen_end_of_header = false;
    time_str = '';
    while ~feof(fid)
        line = fgetl(fid);
        if ~ischar(line), continue; end
        if contains(line, '***End_of_Header***')
            seen_end_of_header = true;
            continue;
        end
        if seen_end_of_header && startsWith(strtrim(line), 'Time')
            parts = split(line, char(9));
            time_str = strrep(parts{2}, '"', '');
            time_str = strrep(time_str, ',', '.');  % comma decimal -> period
            break;
        end
    end
    fclose(fid);
end

function probe_time_str = parse_probe_start_time(filename)
% PARSE_PROBE_START_TIME Extract the probe file's "t0" acquisition start
% timestamp (a differently-formatted header than the LVM files: a single
% "t0" line near the top, with a date+time string per probe channel --
% only the time portion is needed here, since synchronization is done by
% time-of-day, not calendar date).
    fid_probe = fopen(filename, 'r');
    if fid_probe == -1, error('Could not open probe file: %s', filename); end
    while ~feof(fid_probe)
        line = fgetl(fid_probe);
        if startsWith(line, 't0')
            parts = split(line, char(9));
            probe_time_str = strrep(parts{2}, '"', '');
            probe_time_str = strrep(probe_time_str, ',', '.');
            if contains(probe_time_str, ' ')
                subparts = split(probe_time_str, ' ');
                probe_time_str = subparts{end};  % drop the date, keep HH:mm:ss.ffffff
            end
            break;
        end
    end
    fclose(fid_probe);
end

function data = read_european_lvm(filename, dataStartLine)
% READ_EUROPEAN_LVM Import a tab-delimited LabVIEW/probe file that uses
% comma decimal separators, skipping header rows up to dataStartLine.
%   filename      - path to the .lvm file, or the extensionless probe file
%   dataStartLine - first 1-indexed file line containing real numeric data
%                   (23 for Command/Measurement files, 6 for probe files)
% 'FileType','text' is required because detectImportOptions can't infer a
% file type from the .lvm extension (or the probe file's missing
% extension) on its own. DecimalSeparator must be set via setvaropts
% AFTER setvartype forces every column numeric -- it isn't a valid
% constructor argument to detectImportOptions/delimitedTextImportOptions
% directly.
    opts = detectImportOptions(filename, 'FileType', 'text', 'Delimiter', '\t');
    opts = setvartype(opts, opts.VariableNames, 'double');
    opts = setvaropts(opts, 'DecimalSeparator', ',');
    opts.DataLines = [dataStartLine, Inf];
    data = readmatrix(filename, opts);
end