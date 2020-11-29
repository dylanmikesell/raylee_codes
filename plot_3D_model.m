close all; clc;

load('vs_3d_model.mat');

%%

[X,Y] = meshgrid(x,y);

figure;
plot(X,Y,'k*'); xlabel('x-coordinate [m]'); ylabel('y-coordinate [m]');
grid on; axis('tight');

z_idx = 1:5:100;

vmax = max(vs_mod(:));
vmin = min(vs_mod(:));

for ii = 1 : numel(z_idx)
   depth_slice = squeeze( vs_mod(:,:,z_idx(ii)) );
   figure;
   imagesc(x,y,depth_slice); colormap('jet'); grid on;
   xlabel('x-coordinate [m]'); ylabel('y-coordinate [m]');
   caxis([vmin vmax]);
   title( sprintf('Depth: %0.f [m]', depth(z_idx(ii))) );
   c = colorbar; ylabel(c,'V_s [m/s]');
   
end




