%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% PROGRAM:
% make_initial_model_ex1.m
%
% PROGRAMMERS:
% Matt Haney and Victor Tsai
%
% Last revision date:
% 26 April 2017
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This code is distributed as part of the source-code package
%                   raylee_inversion_codes
% that accompanies Haney and Tsai (2017). The package can be downloaded
% from the Geophysics source-code archive at
%                   http://software.seg.org/2017/0003/index.html
% Use of this code is subject to acceptance of the terms and conditions
% that can be found at http://software.seg.org/disclaimer.txt
% Copyright 2017 by The Society of Exploration Geophysicists (SEG)
% Reference:
% Haney, M. M., Tsai, V. C. (2017) Perturbational and nonperturbational
% inversion of Rayleigh-wave velocities, Geophysics, 82(3), F15-F28.
% doi: 10.1190/geo2016-0397.1
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Program make_initial_model_ex1 is a Matlab script to make an initial
% model for use in iterative Rayleigh/Scholte wave phase or group velocity
% inversion.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Output:
%
% files vp_init.txt, vs_init.txt, rho_init.txt, vpf.txt, rhof.txt,
%       grid_values_solid.txt, grid_values_fluid.txt, input_params.txt
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% this is a script
clear all

% Inversion parameters
chilo      = 1;     % low end of acceptable chi-squared window
chihi      = 1.5;   % high end of acceptable chi-squared window
nupdats    = 16;    % maximum number of potential updates
mstdfact   = 2;     % model standard deviation factor
smscl      = 1000;  % smoothing scale
pratioflag = 1;     % flag for if Vp/Vs ratio is fixed (1) or not (0)

% Build the initial model

% fluid properties
vp_flu = 1500;
rho_flu = 1000;

% Elastic solid properties following Gardner's relation
vpvsr = 1.7321; % Vp/Vs ratio
gardc = 309.6; % constant in Gardner relation
powr  = 0.25; % exponent in Gardner relation

% construct a grid in solid
n_sol_elements = 300; % number of elements in solid
h_sol          = 250; % thickness of elements in the solid
% construct a grid in fluid
n_flu_elements = 0; % number of elements in fluid
h_flu          = 100; % thickness of elements in the fluid

% homogeneous initial model
vslay1  = 3400;
vplay1  = vslay1*vpvsr;
rholay1 = gardc*(vplay1^powr);

% Build the elements
h    = h_sol*ones(1,n_sol_elements); % grid spacing of mesh (meters)
vpv  = vplay1*ones(1,n_sol_elements);
vsv  = vslay1*ones(1,n_sol_elements);
rhov = rholay1*ones(1,n_sol_elements);

% Build the fluid elements
hfv   = h_flu*ones(1,n_flu_elements); % grid spacing of mesh (meters)
vpvf  = vp_flu*ones(1,n_flu_elements);
rhovf = rho_flu*ones(1,n_flu_elements);

%% Do Not Edit below this

file_path = 'input_files';
if ~isfolder(file_path)
    mkdir(file_path);
end

Nf = length(load(fullfile(file_path,'velocity_values.txt')));  % number of frequency measurements

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


%% Functions

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

