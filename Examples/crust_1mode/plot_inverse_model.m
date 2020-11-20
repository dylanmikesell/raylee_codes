% This script replaces the original plot_results script. This now uses
% actual output data files from the raylee_invert.m routine. Thus, you do
% not need to keep things in memory in order to plot the results.

clear; close all; clc;


%% Plot the dispersion data

% Load the observed dispersion data
fks = load(fullfile('input_files', 'frequency_values.txt'));
U_data = load(fullfile('input_files', 'velocity_values.txt'));

% Load the initial model dispersion data
fksr_guess = load(fullfile('output_files', 'frequency_guess.txt'));
U_guess = load(fullfile('output_files', 'velocity_guess.txt'));

% Load the final model dispersion data
fks = load(fullfile('output_files', 'frequency_final.txt'));
U = load(fullfile('output_files', 'velocity_final.txt'));

% plot dispersion data comparison
h = figure('Color','w');
plot(fks,U_data,'bo','MarkerSize',6); hold on; % data
plot(fks,U,'ko','MarkerSize',6); % inverted
plot(fksr_guess,U_guess,'ro','MarkerSize',6); % initial
xlabel('Frequency [Hz]'); ylabel('Velocity [m/s]'); 
legend({'data','inverted','initial'},'location','best','box','off');
axis('tight');

% title(' Data (blue), initial guess (red), and final update (black) ');
% fsize = 16;
% axis([.1 .65 1500 3400])
% set(gca,'Fontsize',fsize,'FontWeight','bold');

%% Plot the Vs models 

Nn = 300; % number of solid elements for this example

% parameters
vpvsr = 1.7321;

% make a three layered model
layrth1 = 5; % thickness in elements, layrth1*h = thickness in meters
layrth2 = 10; % thickness in elements, layrth2*h = thickness in meters
layrth3 = 50;

% the true model
vplay1 = 4000; 
vslay1 = vplay1/vpvsr; 
vplay2 = 3396; 
vslay2 = vplay2/vpvsr; 
vplay3 = 4500; 
vslay3 = vplay3/vpvsr; 
vplay4 = 6000; 
vslay4 = vplay4/vpvsr; 

vsv_true = [...
    vslay1*ones(1,layrth1)...
    vslay2*ones(1,layrth2) ...
    vslay3*ones(1,layrth3) ...
    vslay4*ones(1,(Nn-(layrth1+layrth2+layrth3)))...
    ];

vsv_guess = load(fullfile('output_files', 'vs_guess.txt'));
vsv_final = load(fullfile('output_files', 'vs_final.txt'));
hss = load(fullfile('output_files', 'depth_final.txt'));
hss = hss/1000; % [km] convert meters to kilometers

% plot model comparisons
h = figure('Color','w');
plot(vsv_true,hss,'b-','LineWidth',4); hold on; axis ij; % true model
plot(vsv_final,hss,'k--','LineWidth',4); % final model
plot(vsv_guess,hss,'r--','LineWidth',4); % initial model
axis('tight');
xlabel('Shear velocity [m/s]'); ylabel('Depth [km]'); 
legend({'data','inverted','initial'},'location','best','box','off');

% axis([1500 4000 0 30]);
% axis ij; hold on
% fsize = 16;
% title(' True (blue) and initial models (red), and final update (black) ');
% set(gca,'Fontsize',fsize,'FontWeight','bold');

%% Plot the kernels

% Don't worry about this for now.