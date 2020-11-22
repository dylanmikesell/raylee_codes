%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% PROGRAM:
% make_synthetic_ex1.m
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
% Program make_synthetic_ex1 is a Matlab script to make some synthetic
% phase or group velocity data for Rayleigh/Scholte waves. It writes out
% input files for eventual use by the program raylee_invert.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Output:
%
% files velocity_values.txt, velocity_values_errs.txt,
%       frequency_values.txt, mode_values.txt, vtype_values.txt
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% this is a script
clear all

file_path = 'input_files';
if ~isfolder(file_path)
    mkdir(file_path);
end

%% Make the frequencies

% construct vector of frequencies
Nf   = 56;                   % number of measurements
fmin = 0.10;                    % minimum frequency (Hz)
df   = 0.01;                      % frequency spacing (Hz)
fks  = fmin + (df*[0:(Nf-1)]);     % vector of frequencies (Hz)
% these are the frequencies at which the velocities are measured

% % model 2 modes, 1 fundamental and 1st higher mode
% Nf = 2*Nf;
% % vector of mode numbers
% modnv = [ones(1,Nf/2) 2*ones(1,Nf/2)];
% % try to model over same frequencies
% fks = [fks fks];

% % model 1 mode, 1 fundamental
% % vector of mode numbers
modnv = ones(1,Nf);

% what type of velocity - phase (0) or group (1)?
vtypv = ones(1, Nf);

%% Make the vertial finite-element model

% make the fluid part of the model
n_flu_elements = 0; % number of elements in fluid
h_flu = 100; % thickness of elements in the fluid

% make the solid part of the 4-layer model
vplay1  = 4000;
vplay2  = 3396;
vplay3  = 4500;
vplay4  = 6000;

vpvsr = 1.7321; % Vp/Vs ratio

n_sol_elements = 240; % number of elements in solid
h_sol = 250; % thickness of elements in the solid

% make a three layered model
lay1_n_elements = 5; % thickness in elements, layrth1*h = thickness in meters
lay2_n_elements = 10; % thickness in elements, layrth2*h = thickness in meters
lay3_n_elements = 50;
% Note the remaining lay4 has (n_sol_elments - (lay1+lay2+lay3)) elements


[vpv, vsv, rhov, vpfv, rhofv, h, hfv] = make_model(n_sol_elements,h_sol,...
    lay1_n_elements, lay2_n_elements, lay3_n_elements, vplay1, vplay2,...
    vplay3, vplay4, vpvsr, n_flu_elements, h_flu);

% Make synthetic dispersion curve
vout = make_dispersion_curve(fks,n_sol_elements,vsv,vpv,rhov,h,modnv,n_flu_elements,vpfv,rhofv,hfv,vtypv);

% add 1% Gaussian noise
% set randn to its default initial state so test is repeatable
randn('state', 0);
noise = 0.025; % scale of the noise
vout = ( 1 + noise*randn(1,Nf) ) .* vout;

% take out the NaNs
nan_idx = isnan( vout ); % True if vout==NaN
vout(nan_idx)  = []; % remove NaNs
fks(nan_idx)   = []; % remove NaNs
modnv(nan_idx) = []; % remove NaNs
vtypv(nan_idx) = []; % remove NaNs

Nf = length(fks); % Update number of frequencies after removing NaNs

% Write the different single column data files
write_txt_file(file_path, Nf, vout, 'velocity');
write_txt_file(file_path, Nf, noise.*vout, 'errors');
write_txt_file(file_path, Nf, fks, 'frequency');
write_txt_file(file_path, Nf, modnv, 'mode_value');
write_txt_file(file_path, Nf, vtypv, 'disp_type');


%% Functions

% -------------------------------------------------------------------------
function [vpv, vsv, rhov, vpfv, rhofv, h, hfv] = make_model(n_sol_elements,h_sol,...
    lay1_n_elements, lay2_n_elements, lay3_n_elements, vplay1, vplay2,...
    vplay3, vplay4, vpvsr, n_flu_elements, h_flu)


% Note that the fluid is hardcoded here as water
vplayf = 1500;
rholayf = 1000;

% The model in the fluid
vpfv  = vplayf * ones(1,n_flu_elements);
rhofv = rholayf * ones(1,n_flu_elements);
% grid thicknesses in fluid
hfv = h_flu*ones(1, n_flu_elements); % grid spacing of mesh (meters)


% The model in the solid

% elastic solid parameters are also hardcoded here for Gardner's relation
gardc = 309.6; % constant in Gardner relation
powr  = 0.25; % exponent in Gardner relation

% layer 1 Vs and rho
vslay1  = vplay1/vpvsr;
rholay1 = gardc*(vplay1^powr);
% layer 2 Vs and rho
vslay2  = vplay2/vpvsr;
rholay2 = gardc*(vplay2^powr);
% layer 3 Vs and rho
vslay3  = vplay3/vpvsr;
rholay3 = gardc*(vplay3^powr);
% layer 4 Vs and rho
vslay4  = vplay4/vpvsr;
rholay4 = gardc*(vplay4^powr);

% Remaining elements go the last layer
lay4_n_elements = n_sol_elements - (lay1_n_elements+lay2_n_elements+lay3_n_elements);

% Vp
vpv = [...
    vplay1*ones(1,lay1_n_elements)...
    vplay2*ones(1,lay2_n_elements)...
    vplay3*ones(1,lay3_n_elements)...
    vplay4*ones(1,lay4_n_elements)...
    ];
% Vs
vsv = [...
    vslay1*ones(1,lay1_n_elements)...
    vslay2*ones(1,lay2_n_elements)...
    vslay3*ones(1,lay3_n_elements)...
    vslay4*ones(1,lay4_n_elements)...
    ];
% rho
rhov = [...
    rholay1*ones(1,lay1_n_elements)...
    rholay2*ones(1,lay2_n_elements)...
    rholay3*ones(1,lay3_n_elements)...
    rholay4*ones(1,lay4_n_elements)...
    ];

% grid thicknesses in solid
h = h_sol*ones(1, n_sol_elements); % grid spacing of mesh (meters)

end

% -------------------------------------------------------------------------
function vout = make_dispersion_curve(fks,Nn,vsv,vpv,rhov,h,modnv,Nnf,vpfv,rhofv,hfv,vtypv)

% make some synthetic data
countr = 0;
for f=fks
    
    countr = countr + 1;
    modn = modnv(countr);
    
    [kk, vpk, vgk, ev] = ...
        raylee_lysmer(Nn,vsv,vpv,rhov,f,h,modn,Nnf,vpfv,rhofv,hfv);
    
    if (vtypv(countr) == 0) % phase or group velocity
        vout(countr) = vpk;
    else
        vout(countr) = vgk;
    end
    
end

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
