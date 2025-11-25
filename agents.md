---
description: A suite of signal processing tools for analyzing interferometer data, modeling stepper motor mechanics, fitting ideal encoder signals, and performing Zhao et al. linearization for high-precision displacement tracking.
globs: ["**/*.m", "**/*.csv"]
---

# Interferometry and Encoder Analysis Tools

This toolset allows for the extraction of displacement data from interferometer ADC logs, the generation of synthetic ideal encoder signals via least-squares fitting, and the evaluation of sensor resolution using the Zhao linearization algorithm.

## `alpha_trim_filter`

Applies an alpha-trimmed mean filter to a signal to remove impulsive noise while preserving edges.

- **description**: A non-linear smoothing filter that sorts data within a window and averages the central portion, discarding the outliers (alpha-trimming).
- **parameters**:
  - `x`: The input signal array.
  - `k`: The window size (integer).
  - `p`: The percentage of data to trim from the cumulative distribution (0 to 100).
- **returns**: Filtered signal `y` of the same length as `x`.

## `compute_interferometer_physics`

Calculates the theoretical fringe frequency and displacement constants based on the physical setup of the laser and stepper motor.

- **description**: Derives conversion factors (meters per fringe, meters per microstep) and expected frequencies using laser wavelength and mechanical screw pitch.
- **parameters**:
  - `fs`: Sampling frequency (Hz).
  - `lambda`: Laser wavelength (meters).
  - `pitch_mm`: Screw pitch (mm/rev).
  - `steps_per_rev`: Motor full steps per revolution.
  - `microsteps`: Microsteps per full step.
  - `f_stepper`: Actual motor step frequency (microsteps/s).
- **returns**: A structure containing `disp_per_period`, `microstep_size`, `f_ideal` (expected fringe frequency), and `microsteps_per_fringe`.

## `estimate_displacement_hilbert`

Estimates the "Ground Truth" displacement from the raw interferometer ADC signal using the Hilbert transform.

- **description**: Band-passes the raw ADC signal around the expected fringe frequency, computes the analytic signal via Hilbert transform, unwraps the phase, and converts it to displacement.
- **parameters**:
  - `adc_signal`: The raw analog-to-digital converter signal.
  - `fs`: Sampling frequency.
  - `f_ideal`: The expected center frequency of the fringes.
  - `bw`: Bandwidth half-width for the bandpass filter (default 0.8 Hz).
  - `disp_per_period`: Displacement corresponding to $2\pi$ phase change ($\lambda/2$).
- **returns**: Time-series array `disp_t` representing continuous displacement in meters.

## `fit_ideal_encoder_signals`

Generates idealized integer-based Quadrature signals (A and B) by fitting sinusoids to noisy measurement data.

- **description**: Performs a least-squares fit on noisy input signals (`sigA`, `sigB`) against a design matrix constructed from the spatial frequency ($k_{spatial} = 2\pi / 0.75mm$). It then scales the fitted sinusoids to a specific integer range (e.g., 300 to 3796 ADC counts).
- **parameters**:
  - `sigA_raw`: Raw Signal A data.
  - `sigB_raw`: Raw Signal B data.
  - `displacement_mm`: The spatial domain vector (displacement in mm) used for the fitting basis.
  - `spatial_period_mm`: The encoder wavelength (default 0.75 mm).
- **returns**: `sigA_fit_int`, `sigB_fit_int` (Clean, idealized, integer-quantized quadrature signals).

## `zhao_linearization`

Demodulates Quadrature signals (A/B) into angular displacement using the Zhao et al. algorithm.

- **description**: A geometric demodulation method that calculates displacement by linearizing the Lissajous figure formed by signals A and B. It handles signal normalization, segmentation, and arctangent calculation to derive `thetaC`.
- **parameters**:
  - `sigA`: Component A of the quadrature signal (clean or fitted).
  - `sigB`: Component B of the quadrature signal (clean or fitted).
  - `scale_factor`: Factor to convert radians to physical units (e.g., $750 / 2\pi$).
- **returns**: `displacement_zhao` (The calculated displacement array).

## `analyze_resolution_reactivity`

Evaluates the effective resolution and reactivity of the calculated displacement.

- **description**: Analyzes the step-wise behavior of the displacement output. It identifies "plateaus" (where output does not change despite physical movement) and "steps" to determine the worst-case resolution and the system's sensitivity.
- **parameters**:
  - `displacement_calculated`: The output from the Zhao algorithm.
  - `displacement_truth`: The reference displacement (e.g., from the Hilbert transform).
  - `tolerance`: Threshold for considering two values identical (e.g., 1e-9).
- **returns**:
  - `max_plateau`: The largest physical displacement occurred without a change in sensor output (Worst-case resolution).
  - `best_step`: The smallest physical displacement that triggered a sensor output change.
  - `rmse`: Root Mean Square Error between calculated and truth displacement.
