%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% PROGRAM:
% raylee_invert.m
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
% Program raylee_invert is a Matlab script to invert any collection of
% fundamental mode/higher mode Rayleigh wave group or phase
% velocities measured at a set of frequencies for a shear wave velocity
% depth model.
%
% The program can invert for shear velocity assuming the Vp/Vs ratio is
% fixed in the subsurface or assuming that the original Vp is unchanged
% during the inversion.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Input
%
% 1. Input parameters are read from the file "input_params.txt"
%
% 2. Measured (group or phase) velocities are read from
%    the file "velocity_values.txt"
%
% 3. The frequencies at which the velocities are measured are read
%    from the file "frequency_values.txt"
%
% 4. Error bars on the measured (group or phase) velocities are read from
%    the file "velocity_values_errs.txt"
%
% 5. The types of velocities measured (group or phase) are read from
%    the file "vtype_values.txt"
%
% 6. The mode numbers of the measurements (group or phase) are read from
%    the file "mode_values.txt"
%
% 7. The finite-element grid used for modeling and inversion is read
%    from the files "grid_values_solid.txt" and "grid_values_fluid.txt"
%
% 8. The initial model of the subsurface material properties
%    (Vp, Vs, and density) for the inversion is read from the files
%    "vp_init.txt", "vs_init.txt", and "rho_init.txt".
%
% 9. The subsurface material properties (Vp and density) for a water
%    layer above the solid are read from the files "vpf.txt" and
%    "rhof.txt".
%
% Output
%
% 1. The Vs model updates are available in the matrix vsv_update. This
%    matrix is size (nupdat x Nn) where nupdat is the number of iterations
%    executed before the stopping criterion is met and Nn is the number
%    of nodes for the finite-element grid. The first update is in row 1
%    and the final update is in row nupdat.
%
% 2. The sensitivity kernel matrix for the last iteration is the
%    matrix snsmf_vstot. This matrix has size (Nn x Nf), where Nf is the
%    number of frequencies at which measurements exist.
%
% 3. The computed group or phase velocity for the final update is the
%    vector U, size (1 x Nf).
%
% 4. The RMS error for the initial guess and all the updates is stored
%    in the vector rmserror, size (1 x (nupdat+1)).
%
% 5. The Chi-squared error for the initial guess and all updates is stored
%    in the vector chisqurd, size (1 x (nupdat+1)).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% 
% % this is a script
% clear all
% 
% % time this
% tic

function raylee_invert

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% begin input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

file_path = 'raylee_input_files';

% Get the data "PATH/files"
[ifil, efil, ffil, mfil, vfil] = get_data_files(file_path);

% Get the initial model "PATH/files"
[gfil, gffil, nprm, vpinit, vsinit, rhoinit, vpf, rhof] = get_model_files(file_path);

% Read the input parameter file
[pratioflag, lsmth, msigmaf, nupds, Nf, Nn, Nnf, chilo, chihi] = load_input_params(nprm);

% Read the data-related text files
fks = read_txt_file(ffil, Nf); % vector of frequencies at which the velocity is measured [Hz]
U_data = read_txt_file(ifil, Nf); % data vector of measured velocities [m/s]
U_data_errs = read_txt_file(efil, Nf); % data vector of errors [m/s]
modn = read_txt_file(mfil, Nf); % data vector of mode numbers [integer]
vflg = read_txt_file(vfil, Nf); % data vector of velocity types [integer]


% Read the initial model input text files

% Fluid
hfv = read_txt_file(gffil, Nnf); % grid spacing of mesh (m)
vpfv = read_txt_file(vpf, Nnf); % vector of fluid Vp model [m/s]
rhofv = read_txt_file(rhof, Nnf); % vector of fluid density model

% Solid
h = read_txt_file(gfil, Nn); % grid spacing of mesh (m)
vsv = read_txt_file(vsinit, Nn); % vector of initial Vs model [m/s]
rhov = read_txt_file(rhoinit, Nn); % vector of initial density model

if (pratioflag == 0) % load initial Vp model
    vpv = read_txt_file(vpinit, Nn); % vector of initial Vp model
elseif (pratioflag == 1)
    vpvsratio = read_txt_file(vpinit, 1); % Vp/Vs ratio
    vpv = vpvsratio*vsv; % compute the Vp vector
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize velocity update vector
vsv_update = zeros(nupds,Nn);

% make a vector of depths at nodes by a running sum of the grid spacings
hs(1) = 0;
for ii=2:length(h)
    hs(ii) = sum(h(1:(ii-1)));
end

% make a vector of depths at center of elements by a running sum
hss(1) = h(1)/2;
for ii=2:length(h)
    hss(ii) = sum(h(1:(ii-1))) + h(ii)/2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check whether input parameters are physically possible
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% density greater than zero
if (sum(rhov <= 0) > 0)
    error('Negative density values exist in initial guess');
else
end
% shear velocity greater than zero
if (sum(vsv <= 0) > 0)
    error('Negative shear velocity values exist in initial guess');
else
end
% poisson's ratio between two bounds
pratio = (vpv.^2 - 2*(vsv.^2))./(2*(vpv.^2 - vsv.^2));
if ((sum(pratio <= -1) > 0) || (sum(pratio >= 0.5) > 0))
    error('Impossible Poisson ratio values exist in initial guess');
else
end
% density greater than zero in fluid
if (sum(rhofv <= 0) > 0)
    error('Negative density values exist in initial guess');
else
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prepare for initial inversion step
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

output_path = 'raylee_output_files';
if ~isfolder(output_path)
    mkdir(output_path);
end

% compute sensitivity kernel using initial guess
[U, snsmf_vstot] = raylee_sensitivity(Nn,vsv,vpv,...
    rhov,fks,h,modn,vflg,Nnf,vpfv,rhofv,hfv,pratioflag);

% find the measurements for which both data and model are not NaN
[Ur, U_datar, fksr, fksri, modnr, vflgr, snsmf_vstotr] = ...
    check_nans(U, U_data, fks, modn, vflg, snsmf_vstot);
Nfr = length(fksr);

% save the S-wave velocity guess and the resulting data
vsv_guess = vsv;
U_guess = Ur;
fksr_guess = fksr;

% Save the initial model
save_disp_model(U_guess,fksr_guess,output_path,'initial');
save_vs_model(vsv_guess, hss, output_path, 'initial')

% calculate the a priori model covariance matrix and inverse square root
msigma = mean(U_data_errs(fksri))*msigmaf;
mcm = (msigma^2)*exp(-abs(repmat(hs,Nn,1)-repmat(hs',1,Nn))/lsmth);
mcmisr = sqrtm(inv(mcm));

% calculate the a priori data covariance matrix and inverse square root
dcm = diag(U_data_errs(fksri).^2);
dcmisr = diag(1./U_data_errs(fksri));

% rms error of the initial guess
rmserror(1) = sqrt(mean(((U_guess-U_datar)./1).^2));
chisqurd(1) = (U_guess-U_datar)*dcmisr*dcmisr*transpose(U_guess-U_datar);
Nfrv(1) = Nfr;

% check to see if initial guess has chi^2 less than 1
if ((chisqurd(1)/Nfr) < chilo)
    error('Initial model fits data to less than 1 chi-squared');
elseif ((chisqurd(1)/Nfr) < chihi)
    error('Initial model fits data within acceptable chi-squared window');
else
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % invert using damped least squares method of Tarantola and Valette (1982)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('\nRunning the damped least-squares inversion.\n');
    
    % initial damped linear inversion
    dvs = linvers(U_datar,Ur,snsmf_vstotr,mcmisr,dcmisr,Nn,vsv,vsv_guess);
    
    % add to the initial model
    vsv = dvs' + vsv_guess;
    if (pratioflag == 1)
        vpv = vpvsratio*vsv;
    else
    end
    
    % compute new sensitivity kernel
    [U, snsmf_vstot] = raylee_sensitivity(Nn,vsv,vpv,...
        rhov,fksr,h,modnr,vflgr,...
        Nnf,vpfv,rhofv,hfv,pratioflag);
    
    % find NaNs
    [U, U_datar, fksr, fksri, modnr, vflgr, snsmf_vstot] = ...
        check_nans(U, U_datar, fksr, modnr, vflgr, snsmf_vstot);
    
    % if number of NaNs changed, recompute data and model covariances
    if (length(fksr) ~= Nfr)
        Nfr = length(fksr);
        msigma = mean(U_data_errs(fksri))*msigmaf;
        mcm = (msigma^2)*exp(-abs(repmat(hs,Nn,1)-repmat(hs',1,Nn))/lsmth);
        mcmisr = sqrtm(inv(mcm));
        dcm = diag(U_data_errs(fksri).^2);
        dcmisr = diag(1./U_data_errs(fksri));
    else
    end
    
    % compute RMS error and chi-squared
    rmserrorp = sqrt(mean(((U-U_datar)./1).^2));
    chisqurdp = (U-U_datar)*dcmisr*dcmisr*transpose(U-U_datar);
    
    % a reduced line search if chi^2 of update is not lower
    nreds = 0;
    while ((chisqurdp >= chisqurd(1) && nreds < nupds) || ...
            ((chisqurdp/Nfr) < 1 && nreds < nupds))
        
        fprintf('Chi2 increased. Reducing step size by factor of 2.\n');
        nreds = nreds + 1;
        fprintf('Reduction number: %d\n', nreds);
        
        % reduce step by a factor of 2, and add it in
        dvs = dvs/2;
        vsv = vsv_guess + dvs';
        if (pratioflag == 1)
            vpv = vpvsratio*vsv;
        else
        end
        
        % call the sensitivity function to compute U
        [U, snsmf_vstot] = raylee_sensitivity(Nn,vsv,vpv,...
            rhov,fksr,h,modnr,vflgr,...
            Nnf,vpfv,rhofv,hfv,pratioflag);
        
        % check for NaNs
        [U, U_datar, fksr, fksri, modnr, vflgr, snsmf_vstot] = ...
            check_nans(U, U_datar, fksr, modnr, vflgr, snsmf_vstot);
        
        % if number of NaNs changed recompute data and model covariances
        if (length(fksr) ~= Nfr)
            Nfr = length(fksr);
            msigma = mean(U_data_errs(fksri))*msigmaf;
            mcm = (msigma^2)*exp(-abs(repmat(hs,Nn,1)-repmat(hs',1,Nn))/lsmth);
            mcmisr = sqrtm(inv(mcm));
            dcm = diag(U_data_errs(fksri).^2);
            dcmisr = diag(1./U_data_errs(fksri));
        else
        end
        
        % the rms of this potential update
        rmserrorp = sqrt(mean(((U-U_datar)./1).^2));
        % the chi^2 of this potential update
        chisqurdp = (U-U_datar)*dcmisr*dcmisr*transpose(U-U_datar);
        
    end
    
    % shear velocity must be greater than zero
    if (sum(vsv <= 0) > 0)
        error('Negative shear velocity values encountered in inversion');
    else
    end
    % poisson's ratio between two bounds
    pratio = (vpv.^2 - 2*(vsv.^2))./(2*(vpv.^2 - vsv.^2));
    if ((sum(pratio <= -1) > 0) || (sum(pratio >= 0.5) > 0))
        error('Impossible Poisson ratio values encountered in inversion');
    else
    end
    
    % the updated model, print number of update to screen
    nupdat = 1;
    fprintf('\nModel update (nupdat): %d\n', nupdat);
    vsv_update(nupdat,:) = vsv;
    
    % the rms of this update
    rmserror(nupdat+1) = rmserrorp;
    % the chi^2 of this update
    chisqurd(nupdat+1) = chisqurdp;
    
end

% now full modeling
[U, snsmf_vstot] = raylee_sensitivity(Nn,vsv,vpv,...
    rhov,fks,h,modn,vflg,...
    Nnf,vpfv,rhofv,hfv,pratioflag);

% check for NaNs
[Ur, U_datar, fksr, fksri, modnr, vflgr, snsmf_vstotr] = ...
    check_nans(U, U_data, fks, modn, vflg, snsmf_vstot);

% if number of NaNs changed recompute data and model covariances
if (length(fksr) ~= Nfr)
    Nfr = length(fksr);
    msigma = mean(U_data_errs(fksri))*msigmaf;
    mcm = (msigma^2)*exp(-abs(repmat(hs,Nn,1)-repmat(hs',1,Nn))/lsmth);
    mcmisr = sqrtm(inv(mcm));
    dcm = diag(U_data_errs(fksri).^2);
    dcmisr = diag(1./U_data_errs(fksri));
else
end

% compute RMS and chi-squared
rmserrorp = sqrt(mean(((Ur-U_datar)./1).^2));
chisqurdp = (Ur-U_datar)*dcmisr*dcmisr*transpose(Ur-U_datar);

% the rms of this update
rmserror(nupdat+1) = rmserrorp;
% the chi^2 of this update
chisqurd(nupdat+1) = chisqurdp;
Nfrv(nupdat+1) = Nfr;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% now an iterative loop, updating the initial guess
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% while the stopping criterion and the maximum
% allowed number of iterations has not been met, continue updating
fprintf('\nNow running iterative inversion step.\n');
while ((chisqurdp/Nfr) > chihi && nupdat < nupds )
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % invert again as in Tarantola and Valette (1982)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % linear inverse
    dvs = linvers(U_datar,Ur,snsmf_vstotr,mcmisr,dcmisr,Nn,vsv,vsv_guess);
    
    % add to the initial model
    vsv = dvs' + vsv_guess;
    % if fixed vpvs ratio, adjust vp model
    if (pratioflag == 1)
        vpv = vpvsratio*vsv;
    else
    end
    
    % call the sensitivity function to model
    [U, snsmf_vstot] = raylee_sensitivity(Nn,vsv,vpv,...
        rhov,fksr,h,modnr,vflgr,...
        Nnf,vpfv,rhofv,hfv,pratioflag);
    
    % check for NaNs
    [U, U_datar, fksr, fksri, modnr, vflgr, snsmf_vstot] = ...
        check_nans(U, U_datar, fksr, modnr, vflgr, snsmf_vstot);
    
    % if number of data changed, recompute data and model covariances
    if (length(fksr) ~= Nfr)
        Nfr = length(fksr);
        msigma = mean(U_data_errs(fksri))*msigmaf;
        mcm = (msigma^2)*exp(-abs(repmat(hs,Nn,1)-repmat(hs',1,Nn))/lsmth);
        mcmisr = sqrtm(inv(mcm));
        dcm = diag(U_data_errs(fksri).^2);
        dcmisr = diag(1./U_data_errs(fksri));
    else
    end
    
    % compute rms and chi of this potential update
    rmserrorp = sqrt(mean(((U-U_datar)./1).^2));
    chisqurdp = (U-U_datar)*dcmisr*dcmisr*transpose(U-U_datar);
    
    % a reduced line search if chi^2 of update is not lower
    nreds = 0;
    % the gradient - difference between the current update and previous
    dvs = (vsv' - transpose(vsv_update(nupdat,:)));
    
    % Check chi2
    while ((chisqurdp >= 1.01*chisqurd(nupdat+1) && nreds < nupds) || ...
            ((chisqurdp/Nfr) < chilo && nreds < nupds))
        
        fprintf('Chi2 increased. Reducing step size by factor of 2.\n');
        nreds = nreds + 1;
        fprintf('Reduction number: %d\n', nreds);
        
        % reduce step by a factor of 2, and add it in
        dvs = dvs/2;
        vsv = vsv_update(nupdat,:) + dvs';
        % if vpvs ratio fixed, adjust vp model
        if (pratioflag == 1)
            vpv = vpvsratio*vsv;
        else
        end
        
        % call the sensitivity function to compute U
        [U, snsmf_vstot] = raylee_sensitivity(Nn,vsv,vpv,...
            rhov,fksr,h,modnr,vflgr,...
            Nnf,vpfv,rhofv,hfv,pratioflag);
        
        % check for NaNs
        [U, U_datar, fksr, fksri, modnr, vflgr, snsmf_vstot] = ...
            check_nans(U, U_datar, fksr, modnr, vflgr, snsmf_vstot);
        %Nfr = length(fksr);
        
        % if number of data changed, adjust model and data covariances
        if (length(fksr) ~= Nfr)
            Nfr = length(fksr);
            msigma = mean(U_data_errs(fksri))*msigmaf;
            mcm = (msigma^2)*exp(-abs(repmat(hs,Nn,1)-repmat(hs',1,Nn))/lsmth);
            mcmisr = sqrtm(inv(mcm));
            dcm = diag(U_data_errs(fksri).^2);
            dcmisr = diag(1./U_data_errs(fksri));
        else
        end
        
        % the rms of this potential update
        rmserrorp = sqrt(mean(((U-U_datar)./1).^2));
        % the chi^2 of this potential update
        chisqurdp = (U-U_datar)*dcmisr*dcmisr*transpose(U-U_datar);
        
    end
    
    % shear velocity must be greater than zero
    if (sum(vsv <= 0) > 0)
        error('Negative shear velocity values encountered in inversion');
    else
    end
    % poisson's ratio between two bounds
    pratio = (vpv.^2 - 2*(vsv.^2))./(2*(vpv.^2 - vsv.^2));
    if ((sum(pratio <= -1) > 0) || (sum(pratio >= 0.5) > 0))
        error('Impossible Poisson ratio values encountered in inversion');
    else
    end
    
    % the next updated model, print number of update to screen
    nupdat = nupdat + 1;
    fprintf('\nModel update (nupdat): %d\n', nupdat);
    vsv_update(nupdat,:) = vsv;
    
    % the rms of this update
    rmserror(nupdat+1) = rmserrorp;
    % the chi^2 of this update
    chisqurd(nupdat+1) = chisqurdp;
    
    % now full modeling
    [U, snsmf_vstot] = raylee_sensitivity(Nn,vsv,vpv,...
        rhov,fks,h,modn,vflg,...
        Nnf,vpfv,rhofv,hfv,pratioflag);
    
    % check for NaNs
    [Ur, U_datar, fksr, fksri, modnr, vflgr, snsmf_vstotr] = ...
        check_nans(U, U_data, fks, modn, vflg, snsmf_vstot);
    %Nfr = length(fksr);
    
    % if number of data changed, recompute data and model covariances
    if (length(fksr) ~= Nfr)
        Nfr = length(fksr);
        msigma = mean(U_data_errs(fksri))*msigmaf;
        mcm = (msigma^2)*exp(-abs(repmat(hs,Nn,1)-repmat(hs',1,Nn))/lsmth);
        mcmisr = sqrtm(inv(mcm));
        dcm = diag(U_data_errs(fksri).^2);
        dcmisr = diag(1./U_data_errs(fksri));
    else
    end
    
    % compute rms and chi^2
    rmserrorp = sqrt(mean(((Ur-U_datar)./1).^2));
    chisqurdp = (Ur-U_datar)*dcmisr*dcmisr*transpose(Ur-U_datar);
    
    % the rms of this update
    rmserror(nupdat+1) = rmserrorp;
    % the chi^2 of this update
    chisqurd(nupdat+1) = chisqurdp;
    Nfrv(nupdat+1) = Nfr;
    
end

% Save the final dispersion data model
save_disp_model(U,fks,output_path,'final');
% Save the final Vs model
save_vs_model(vsv_update(nupdat,:), hss, output_path, 'final');

% end the timer
% toc

sprintf('%d of %d measurements used',Nfr,Nf-sum(isnan(U_data)))

if ((chisqurd(nupdat+1)/Nfr) > chihi)
    sprintf('WARNING: Inversion did not converge to stopping criterion and underfitted data. Increase number of updates.')
else
end

if ((chisqurd(nupdat+1)/Nfr) < chilo)
    sprintf('WARNING: Inversion did not converge to stopping criterion and overfitted data. Increase number of reduction steps.')
else
end

end


%% Functions

% -------------------------------------------------------------------------
function [ifil, efil, ffil, mfil, vfil] = get_data_files(file_path)

% data file of velocities at each frequency
ifil = sprintf(fullfile(file_path,'velocity_values.txt'));
% data file of velocity error bars at each frequency
efil = sprintf(fullfile(file_path,'velocity_values_errs.txt'));
% data file of frequencies
ffil = sprintf(fullfile(file_path,'frequency_values.txt'));
% data file of mode number at each frequency
mfil = sprintf(fullfile(file_path,'mode_values.txt'));
% data file of mode number at each frequency
vfil = sprintf(fullfile(file_path,'vtype_values.txt'));

end
% -------------------------------------------------------------------------
function [gfil, gffil, nprm, vpinit, vsinit, rhoinit, vpf, rhof] = get_model_files(file_path)

% data file of element thicknesses
gfil = sprintf(fullfile(file_path,'grid_values_solid.txt'));

% data file of element thicknesses
gffil = sprintf(fullfile(file_path,'grid_values_fluid.txt'));

% data file of input parameters
nprm = sprintf(fullfile(file_path,'input_params.txt'));

% data files of initial model guess
vpinit = sprintf(fullfile(file_path,'vp_init.txt'));
vsinit = sprintf(fullfile(file_path,'vs_init.txt'));
rhoinit = sprintf(fullfile(file_path,'rho_init.txt'));

% data files of water layer
vpf = sprintf(fullfile(file_path,'vpf.txt'));
rhof = sprintf(fullfile(file_path,'rhof.txt'));

end

% -------------------------------------------------------------------------
function [pratioflag, lsmth, msigmaf, nupds, Nf, Nn, Nnf, chilo, chihi] = load_input_params(nprm)

% load input parameters
inp = load(nprm);

% flag to indicate if vp/vs ratio should be fixed
pratioflag = inp(1);
% inversion parameters
lsmth = inp(2);                     % model smoothness scale (m)
msigmaf = inp(3);                   % model standard deviation factor
% a factor times the mean
% data standard deviation
% stopping criteria
nupds = inp(4);                     % max number of updates (iterations)

% data and model size
Nf = inp(5);                        % number of measurements
Nn = inp(6);                        % number of elements/nodes for solid
Nnf = inp(7);                       % number of elements for fluid

% chi squared window
chilo = inp(8);
chihi = inp(9);

end


% -------------------------------------------------------------------------
function d = read_txt_file(filename, N)

d = zeros(1,N);
fid = fopen(filename,'r');
for ii=1:N
    d(ii) = fscanf(fid,'%f',1);
end
fclose(fid);
fprintf('Done reading %s\n', filename);

end

% -------------------------------------------------------------------------
function save_disp_model(U, fks, output_path, model_type)

switch model_type
    case 'initial'
        vel_file = fullfile(output_path,'velocity_guess.txt');
        freq_file = fullfile(output_path,'frequency_guess.txt');
    case 'final'
        vel_file = fullfile(output_path,'velocity_final.txt');
        freq_file = fullfile(output_path,'frequency_final.txt');
end

% Save the initial model
fid = fopen(vel_file,'w');
for ii = 1 : numel(U)
    fprintf(fid,'%10.5f\n',U(ii)); % Initial dispersion data   
end
fclose(fid);
fprintf('Done writing %s\n', vel_file);

fid = fopen(freq_file,'w');
for ii = 1 : numel(fks)
    fprintf(fid,'%10.5f\n',fks(ii)); % Initial frequencies  
end
fclose(fid);
fprintf('Done writing %s\n', freq_file);

end

% -------------------------------------------------------------------------
function save_vs_model(vs, hss, output_path, model_type)

switch model_type
    case 'initial'
        vel_file = fullfile(output_path,'vs_guess.txt');
        depth_file = fullfile(output_path,'depth_guess.txt');
    case 'final'
        vel_file = fullfile(output_path,'vs_final.txt');
        depth_file = fullfile(output_path,'depth_final.txt');
end

% Save the initial model
fid = fopen(vel_file,'w');
for ii = 1 : numel(vs)
    fprintf(fid,'%10.5f\n',vs(ii)); % Initial dispersion data   
end
fclose(fid);
fprintf('Done writing %s\n', vel_file);

fid = fopen(depth_file,'w');
for ii = 1 : numel(hss)
    fprintf(fid,'%10.5f\n',hss(ii)); % Initial frequencies  
end
fclose(fid);
fprintf('Done writing %s\n', depth_file);

end