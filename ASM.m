%%
%parameters
clear all; clc;

c = 299792458;
lambda0 = 1.55e-6;
f = c/lambda0;
omega = 2*pi*f;
k0 = 2*pi/lambda0;

n = 1 ; %free space
lambda = lambda0/n;
k = k0*n;
%%
%Define input source field

%load file
S=load('f_50.mat');
Ey=S.Ey;
E_i = Ey(:,:,1,2);

%Define 
E_i = E_i(2:end,:)';
source_size = size(E_i);
%% Input Field window

% nx and ny have to be odd, so that no artifacts are introduced when calculating the discrete Fourier transform
% If they are even, the propagated beam will be shifted from its center, for example
nx=source_size(2); % Number of grid points in the x direction
ny=source_size(1); % Number of grid points in the y direction 
x_width = max(S.x)-min(S.x); % Computation window width in the x direction
y_width = max(S.y)-min(S.y); % Computation window width in the y direction
xA = S.x'; % Array of window points in the x direction
yA = S.y'; % Array of window points in the y direction
[XA,YA]=meshgrid(xA,yA); % Aperture grid


%%

%Simulation window
dx=0.02e-6; % Spacing between sucessive points in the x direction
dy=0.02e-6; % Spacing between sucessive points in the y direction

x_width_s = 3*x_width;
y_width_s = 3*y_width;
xA_s = cat(2,[xA(1)-1*x_width:dx:xA(1)],xA(2:end-1),[xA(end):dx:xA(end)+1*x_width]);
yA_s = cat(2,[yA(1)-1*y_width:dy:yA(1)],yA(2:end-1),[yA(end):dy:yA(end)+1*y_width]);
[XA_s,YA_s] = meshgrid(xA_s,yA_s);

obj_size=size(XA_s); 

%% Wave-vector domain parameters

Fs_x = obj_size(2)/x_width_s; % Fx axis range
Fs_y = obj_size(1)/y_width_s; % Fy axis range
dFx = Fs_x/(obj_size(2)-1); % Sampling spacing in the Fx axis
dFy = Fs_y/(obj_size(1)-1); % Sampling spacing in the Fy axis
Fx=-Fs_x/2:dFx:Fs_x/2; % Points in Fx direction
Fy=-Fs_y/2:dFy:Fs_y/2; % Points in the Fy direction
[FX,FY]=meshgrid(Fx,Fy); % Grid in wave-vector space
% alpha and beta (wavenumber components) 
alpha = lambda.*FX;
beta2 = lambda.*FY; 
gamma_cust = sqrt(1 - alpha.^2 - beta2.^2);

%%
%Combine source field into Computational window

E0 = zeros(obj_size(1),obj_size(2));

nx_start = ceil(obj_size(2)/2)-floor(nx/2) ;
nx_end = ceil(obj_size(2)/2)+floor(nx/2) ;
ny_start = ceil(obj_size(1)/2)-floor(ny/2) ;
ny_end = ceil(obj_size(1)/2)+floor(ny/2);

E0(ny_start:ny_end,nx_start:nx_end)=E_i;


figure(1)
imagesc(abs(E0))
title('Source field')

%%
%Agular spectrum method
A0 = fftshift(fft2(E0));

propagation_distance = 100e-6;
nz = 41;

z = linspace(0,propagation_distance,nz);
E = zeros(obj_size(1),obj_size(2),nz);

for i=1:nz
    %transfer function for free space
    T = exp(1i*k.*gamma_cust.*(z(1,i)));
    
    %band limiting
    ulx = 1/(((2*dFx*z(1,i))^2+1)^0.5*lambda);
    uly = 1/(((2*dFy*z(1,i))^2+1)^0.5*lambda);
    [~, nx_l] = min(abs(abs(Fx)-ulx));
    [~, ny_l] = min(abs(abs(Fy)-uly));
    if -ulx > Fx(1,nx_l)
        nx_l = nx_l+1;
    end    
    if -uly > Fx(1,ny_l)
        ny_l = ny_l+1;
    end
    
    Band_Limit_Matrix = zeros(obj_size(1),obj_size(2));
    Band_Limit_Matrix(nx_l:obj_size(2)-nx_l+1,ny_l:obj_size(1)-ny_l+1)=1;
    
    T = T.*Band_Limit_Matrix;
    
    %Calculation
    E(:,:,i) = ifft2(ifftshift(A0.*T));
    %figure(i+1)
    %imagesc(abs(E(nx_start:nx_end,ny_start:ny_end,i)))
end
%EE = abs(E(nx_start:nx_end,ny_start:ny_end,:));


%%
h = figure;
axis tight manual % this ensures that getframe() returns a consistent size
filename = 'testAnimated.gif';
for n = 1:nz
 
    
    imagesc(abs(E(:,:,n)))
    title(['z = ',num2str(z(1,n))])
    drawnow 
    % Capture the plot as an image 
    frame = getframe(h); 
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256); 
    % Write to the GIF File 
    if n == 1 
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
    else 
        imwrite(imind,cm,filename,'gif','WriteMode','append'); 
    end 
end

%}