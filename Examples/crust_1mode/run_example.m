clear; close all; clc;

% Must point to the src directory
addpath('/Users/dmikesell/GIT/raylee_codes/src');

% -------------------------------------------------------------------------
% This script runs example 1 from Haney & Tsai (2017) with just the
% fundamental mode instead of two modes.
% -------------------------------------------------------------------------

% Make the synthetic data w/ noise
make_synthetic_ex1
% Creates files in input_files/ folder: 
% velocity_values.txt      The measured Rayleigh velocities
% velocity_values_errs.txt Error bars on the measurements
% frequency_values.txt     Frequencies at which the measurements are made
% mode_values.txt          Mode number (Fundamental=1, First overtone=2, etc)
% vtype_values.txt         Vector of velocity type, either phase (0) or group (1)


% Make intial model for inversion
make_initial_model_ex1
% Creates files in input_files/ folder: 
% vp_init.txt               Initial P-wave velocity model in solid
% vs_init.txt               Initial S-wave velocity model in solid
% rho_init.txt              Initial density model in solid
% vpf.txt                   P-wave velocity in fluid layer (if no water layer, this is not used)
% rhof.txt                  Density in fluid layer (if no water layer, this is not used)
% grid_values_solid.txt     Finite element grid in solid layer
% grid_values_fluid.txt     Finite element grid in fluid layer (if no water layer, this is not used)
% input_params.txt          See description below


% Run the inversion 
raylee_invert  % this file is in the src/ directory
% Right now this does not write a solution to file. Instead everything is
% in memory. So call the plotting code now.

plot_results_ex1


