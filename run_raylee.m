clear; close all; clc;

tic 

% This script will run the Raylee inversion code on a complete 2D group
% velocity tomography data set.

% NOTE: To change the initial model and inversion parameters see the
% function: write_initial_model()

% Must point to the Raylee src directory (this must be Dylan's fork from
% github: https://github.com/dylanmikesell/raylee_codes).
addpath('/Users/dmikesell/GIT/raylee_codes/src');

% This should point to the directory where your tomo*.txt files are located
map_dir = 'TOMO_MAPS';

% Set the frequency limits of your 1D dispersion curve inversion.
fmax = 20; % [Hz] highest frequency to invert
fmin = 6; % [Hz] lowest frequency to invert


% First, we need to get the grid information.
filter = '14';
comp = 'ZZ';
gridfile = fullfile('TOMO_FILES', filter, comp, 'XYgrid.dat');
grid = load(gridfile);

%% Start loading the 2D tomography maps

% set up the grid x and y axes using the .dat file from the tomography.
y = grid(2,1) : grid(3,2) : grid(2,2);
x = grid(1,1) : grid(3,1) : grid(1,2);

% get the tomography maps at each frequency
map_files = dir(fullfile(map_dir, 'tomo*.txt'));

% get periods from filenames
period = zeros(1,numel( map_files ));
for ii = 1 : numel( map_files )
    file = map_files(ii).name;
    [~,file,~] = fileparts(file);
    tmp = split(file,'_');
    period(ii) = str2double(tmp{2}(1:end-1));
end
% [period,p_sort_idx] = sort(period)

% load the maps into a 3D matrix
maps = zeros(numel(y),numel(x),numel(period) ); % allocate the maps
for ii = 1 : numel( map_files )
    maps(:,:,ii) = load( fullfile( map_files(ii).folder, map_files(ii).name ) );
end

%% Now we need to invert each point in space


% Use a double for loop to go through pixels
% parfor ii = 1 : numel(y) % use this to make a parallel loop
for ii = 10 % : numel(y)
    parfor jj = 1 : numel(x)
        node = (ii-1)*numel(y) + jj; % keep track of this pixel/node numbre
        per = period; % make a local copy of the periods
        disp = squeeze( maps(ii,jj,:) ); % pull ouf the dispersion curve at this pixel
        
        % Step 1: remove NaNs
        per(isnan(disp)) = []; % must remove periods
        disp(isnan(disp)) = []; % must remove velocity

        % Step 2: if there are more than 5 points in the dispersion cuve,
        % then invert it
        if numel(disp) > 5
            fprintf('Inverting dispersion curve at node %d (x=%0.0f, y=%0.0f)\n',node,x(jj),y(ii));
            
            % preprocess dispersion curve
            [disp, freq] = prepare_disp_curve( per, disp, fmin, fmax );
            
            % do this in a node-based way for parallel purposes
            write_initial_model( numel(freq), disp, node ); % build initial model for the inversion
            write_raylee_data( freq, disp, node ); % write the pixel data
            raylee_invert( node );  % this file is in the src/ directory
        else
            fprintf('Not enough data to invert dispersion curve at node %d.\n', node);
        end
        
    end
end

% Save the model matrix
% vs_mod = build_3d_Vs_matrix(x, y);

toc

%% Functions

% -------------------------------------------------------------------------
function vs_mod = build_3d_Vs_matrix(x, y)
% load the data into a 3D matrix and save as MAT file

for ii = 1 : numel(y)
    for jj = 1: numel(x)
        node = (ii-1)*numel(y) + jj;
        
        % process the inversion output for this node.
        output_path = 'raylee_output_files';
        output_path = fullfile(output_path,num2str(node,'%05d'));
        
        if isfolder(output_path) % make sure this node folder exists
            depth = load( fullfile(output_path,'depth_final.txt') );
            vs = load( fullfile(output_path,'vs_final.txt') );
            
            if ~exist('vs_mod','var') % allocate this variable after first iteration
                vs_mod = nan( numel(y), numel(x), numel(depth) );
            end
            vs_mod(ii,jj,:) = vs;
        end
    end
end

% write the MAT file
save('vs_3d_model.mat','vs_mod','x','y','depth');

end
% -------------------------------------------------------------------------
function [disp_new, freq] = prepare_disp_curve( per, disp, fmin, fmax )

% Convert units as needed for Raylee
disp = disp.*1000; % [m/s] convert from km/s to m/s

% Step 2: set frequency limits and sort for Raylee
freq = 1 ./ per; % [Hz]
disp(freq>fmax) = []; % remove data above fmax
freq(freq>fmax) = [];
disp(freq<fmin) = []; % remove data below fmin
freq(freq<fmin) = [];
[freq,sort_idx] = sort(freq,'ascend');
disp = disp(sort_idx);

% Step 3: smooth the disperion curve using 5th order polynomial approx.
p_coeff = polyfit(freq, disp, 5);
disp_new = polyval(p_coeff, freq);

end

% -------------------------------------------------------------------------
function write_raylee_data( freq, disp, node )

% User parameters to set
noise = 0.025; % noise estimate (0.05 = 5% of velocity value)
% end user parameters to set

% Creates files in raylee_input_files/ folder:
% velocity_values.txt      The measured Rayleigh velocities
% velocity_values_errs.txt Error bars on the measurements
% frequency_values.txt     Frequencies at which the measurements are made
% mode_values.txt          Mode number (Fundamental=1, First overtone=2, etc)
% vtype_values.txt         Vector of velocity type, either phase (0) or group (1)

Nf = numel(freq);
modnv = ones(1,Nf);  % mode, 1 fundamental
vtypv = ones(1, Nf); % what type of velocity: group=1

file_path = 'raylee_input_files';
if ~isfolder(file_path)
    mkdir(file_path);
end
% append the node number
file_path = fullfile(file_path,num2str(node,'%05d'));
if ~isfolder(file_path)
    mkdir(file_path);
end

% Write the different single column data files
write_txt_file(file_path, Nf, disp, 'velocity');
write_txt_file(file_path, Nf, noise.*disp, 'errors');
write_txt_file(file_path, Nf, freq, 'frequency');
write_txt_file(file_path, Nf, modnv, 'mode_value');
write_txt_file(file_path, Nf, vtypv, 'disp_type');

end
% -------------------------------------------------------------------------
function write_initial_model( Nf, disp, node )

% =========================================================================
% USER parameters to set.

% Inversion parameters
chilo      = 1;     % low end of acceptable chi-squared window
chihi      = 1.5;   % high end of acceptable chi-squared window
nupdats    = 15;    % maximum number of potential updates
mstdfact   = 3;     % model standard deviation factor (was 2)
smscl      = 10;  % smoothing scale (was 1000)
pratioflag = 1;     % flag for if Vp/Vs ratio is fixed (1) or not (0)

% Build the initial solid part of the model

% NOTE: You can read "Appendix D: Optimal Layers for Rayleigh Waves" in
% Haney and Tsai (2017). The finite-element grid needs to be the same for
% all dispersion curves, so you only need to set this once after doing some
% tests.

% construct a grid in solid
n_sol_elements = 300; % number of elements in solid
h_sol          = 0.2; % [m] thickness of elements in the solid

% Elastic solid properties following Gardner's relation
vpvsr = 1.7321; % Vp/Vs ratio
gardc = 309.6; % constant in Gardner relation
powr  = 0.25; % exponent in Gardner relation

% homogeneous initial model --> need to set the Vs velocity in this model
% vslay1  = 400; % [m/s]
% vslay1  = (1/0.92)*mean(disp); % make the initial velocity
vslay1  = max(disp);
vplay1  = vslay1*vpvsr;
rholay1 = gardc*(vplay1^powr);

% END USER parameters to set.
% =========================================================================

% fluid grid properties --> zero fluid elements for now
vp_flu = 1500;
rho_flu = 1000;
% construct a grid in fluid
n_flu_elements = 0; % number of elements in fluid
h_flu          = 100; % thickness of elements in the fluid

% Build the elements
h    = h_sol*ones(1,n_sol_elements); % grid spacing of mesh (meters)
vpv  = vplay1*ones(1,n_sol_elements);
vsv  = vslay1*ones(1,n_sol_elements);
rhov = rholay1*ones(1,n_sol_elements);

% Build the fluid elements
hfv   = h_flu*ones(1,n_flu_elements); % grid spacing of mesh (meters)
vpvf  = vp_flu*ones(1,n_flu_elements);
rhovf = rho_flu*ones(1,n_flu_elements);

% Do Not Edit below this

file_path = 'raylee_input_files';
if ~isfolder(file_path)
    mkdir(file_path);
end
% append the node number
file_path = fullfile(file_path,num2str(node,'%05d'));
if ~isfolder(file_path)
    mkdir(file_path);
end

% Nf = length(load(fullfile(file_path,'velocity_values.txt')));  % number of frequency measurements

% write the parameter file
write_param_file(file_path,pratioflag,smscl,mstdfact,nupdats,...
    Nf, n_sol_elements, n_flu_elements, chilo, chihi)

% write out Vp model in single column format
if (pratioflag == 1)
    write_txt_file(file_path, 1, vpvsr, 'vp');
else
    write_txt_file(file_path, length(vpv), vpv, 'vp');
end

% Write the different single column data files
write_txt_file(file_path, length(vsv), vsv, 'vs');
write_txt_file(file_path, length(rhov), rhov, 'density');
write_txt_file(file_path, length(vpvf), vpvf, 'vp_fluid');
write_txt_file(file_path, length(rhovf), rhovf, 'rho_fluid');
write_txt_file(file_path, n_sol_elements, h, 'solid_grid');
write_txt_file(file_path, n_flu_elements, hfv, 'fluid_grid');

end

% -------------------------------------------------------------------------
function write_param_file(file_path,pratioflag,smscl,mstdfact,nupdats,...
    Nf, Nn, Nnf, chilo, chihi)

filename = fullfile(file_path,'input_params.txt');

% write out input parameters file
fidt = fopen(filename,'w');
fprintf(fidt,'%% input parameters for Rayleigh/Scholte wave inversion\n');
fprintf(fidt,'\n');
fprintf(fidt,'%i  %% flag for fixed poisson''s ratio (0=no,1=yes)\n',pratioflag);
fprintf(fidt,'%10.5f  %% smoothness scale (m)\n',smscl);
fprintf(fidt,'%10.5f  %% a priori model standard deviation factor\n',mstdfact);
fprintf(fidt,'%i  %% maximum number of updates (iterations)\n',nupdats);
fprintf(fidt,'%i  %% number of measurements\n',Nf);
fprintf(fidt,'%i  %% number of elements in solid part of model\n',Nn);
fprintf(fidt,'%i  %% number of elements in fluid part of model\n',Nnf);
fprintf(fidt,'%10.5f  %% lower chi squared window\n',chilo);
fprintf(fidt,'%10.5f  %% higher chi squared window\n',chihi);
fclose(fidt);
fprintf('Done writing %s\n', filename);

end
% -------------------------------------------------------------------------
function write_txt_file(file_path, Nf, data, data_type)
%
% Write single column text files for each type of data.

switch data_type
    case 'velocity'
        filename = fullfile(file_path,'velocity_values.txt');
    case 'frequency'
        filename = fullfile(file_path,'frequency_values.txt');
    case 'mode_value'
        filename = fullfile(file_path,'mode_values.txt');
    case 'disp_type'
        filename = fullfile(file_path,'vtype_values.txt');
    case 'errors'
        filename = fullfile(file_path,'velocity_values_errs.txt');
    case 'vp'
        filename = fullfile(file_path,'vp_init.txt');
    case 'vs'
        filename = fullfile(file_path,'vs_init.txt');
    case 'density'
        filename = fullfile(file_path,'rho_init.txt');
    case 'vp_fluid'
        filename = fullfile(file_path,'vpf.txt');
    case 'rho_fluid'
        filename = fullfile(file_path,'rhof.txt');
    case 'solid_grid'
        filename = fullfile(file_path,'grid_values_solid.txt');
    case 'fluid_grid'
        filename = fullfile(file_path,'grid_values_fluid.txt');
end

% write out data for each frequency
fid = fopen(filename,'w');
for ii = 1 : Nf
    fprintf(fid,'%10.5f\n',data(ii));
end
fclose(fid);
fprintf('Done writing %s\n', filename);

end
% -------------------------------------------------------------------------
