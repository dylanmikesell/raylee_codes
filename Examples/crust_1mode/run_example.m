clear; close all; clc;

% Must point to the src directory
addpath('/Users/dmikesell/GIT/raylee_codes/src');

% -------------------------------------------------------------------------
% This script runs example 1 from Haney & Tsai (2017) with just the
% fundamental mode instead of two modes.
% -------------------------------------------------------------------------

% Make the synthetic data
make_synthetic_ex1

% Creates files in input_files/ folder: 
% velocity_values.txt, 
% velocity_values_errs.txt, 
% frequency_values.txt, 
% mode_values.txt,
% vtype_values.txt

% Make intial model for inversion
make_initial_model_ex1

% Creates files in input_files/ folder: 
% vp_init.txt, 
% vs_init.txt, 
% rho_init.txt, 
% vpf.txt, 
% rhof.txt, 
% grid_values_solid.txt, 
% grid_values_fluid.txt, 
% input_params.txt


% Run the inversion 
raylee_invert  % this file is in the src/ directory
% Right now this does not right the solution to file. Instead everything is
% in memory. So call the plotting code now.

plot_results_ex1



