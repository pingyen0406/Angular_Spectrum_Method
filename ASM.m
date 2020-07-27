% parameters
clear all;
lambda0 = 1.55e-6;
k0 = 2*pi/lambda0;
n = 1 ; %free space
lambda = lambda0/n;
k = k0*n;
f = 50e-6;

% Define input source field

% load input E-field file
S=load('80nmSiO2_f_50.mat');
Ey=S.Ey;
E_i = Ey(:,:,1,2);
E_i = E_i(2:end,:)';
source_size = size(E_i);
%%%%%%%%%%%%%%%%%%%%%%%%%%% Input Field window %%%%%%%%%%%%%%%%%%%%%%%%%%%
dx=0.02e-6; % Spacing outside the source field in the x direction
dy=0.02e-6; % Spacing outside the source field in the y direction
% nx and ny have to be odd, so that no artifacts are introduced when calculating the discrete Fourier transform
% If they are even, the propagated beam will be shifted from its center, for example
nx=source_size(2); % Number of grid points in the x direction
ny=source_size(1); % Number of grid points in the y direction 
x_width = max(S.x)-min(S.x); % Computation window width in the x direction
y_width = max(S.y)-min(S.y); % Computation window width in the y direction
xA = S.x'; % Array of window points in the x direction
yA = S.y'; % Array of window points in the y direction
%xA_n = (min(xA)):dx:(max(xA));
%yA_n = (min(yA)):dy:(max(yA));
%[XA_n,YA_n] = ndgrid(xA_n,yA_n);
%Ei_n = griddedInterpolant(XA_n,YA_n,E_i);
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Simulation window %%%%%%%%%%%%%%%%%%%%%%%%%%%%

dx_mid = zeros(1,length(xA));% spacing of the source field in x direction
dy_mid = zeros(1,length(yA));% spacing of the source field in y direction
for i=1:length(xA)-1
    dx_mid(i) = xA(i+1)-xA(i);
end
for i=1:length(yA)-1
    dy_mid(i) = yA(i+1)-yA(i);
end
padding_factor = 7;
x_width_s = padding_factor*x_width;
y_width_s = padding_factor*y_width;
xA_s = cat(2,[xA(1)-(padding_factor-1)/2*x_width:dx:xA(1)],xA(2:end-1),...
    [xA(end):dx:xA(end)+(padding_factor-1)/2*x_width]);
yA_s = cat(2,[yA(1)-(padding_factor-1)/2*y_width:dy:yA(1)],yA(2:end-1),...
    [yA(end):dy:yA(end)+(padding_factor-1)/2*y_width]);
obj_size=size(meshgrid(xA_s,yA_s)); 

% Combine source field into Computational window
E0 = zeros(obj_size(1),obj_size(2));
nx_start = ceil(obj_size(2)/2)-floor(nx/2) ;
nx_end = ceil(obj_size(2)/2)+floor(nx/2) ;
ny_start = ceil(obj_size(1)/2)-floor(ny/2) ;
ny_end = ceil(obj_size(1)/2)+floor(ny/2);
E0(ny_start:ny_end,nx_start:nx_end)=E_i;
% plot input source field
figure(1);
imagesc(xA_s,yA_s,abs(E0));
title('Source field');

%%%%%%%%%%%%%%%%%%%%% Wave-vector domain parameters %%%%%%%%%%%%%%%%%%%%%%%

Fs_x = obj_size(2)/x_width_s; % Fx axis range
Fs_y = obj_size(1)/y_width_s; % Fy axis range 
dFx = Fs_x/(obj_size(2)-1); % Sampling spacing in the Fx axis
dFy = Fs_y/(obj_size(1)-1); % Sampling spacing in the Fy axis
Fx=-Fs_x/2:dFx:Fs_x/2; % Points in Fx direction
Fy=-Fs_y/2:dFy:Fs_y/2; % Points in the Fy direction
[FX,FY]=meshgrid(Fx,Fy); % Grid in wave-vector space
% alpha and beta (wavenumber components) 
alpha = lambda.*FX;
beta = lambda.*FY; 
gamma_cust = sqrt(1 - alpha.^2 - beta.^2);
clear FX FY;
F_struct = struct('Fx',Fx,'Fy',Fy,'dFx',dFx,'dFy',dFy);

%% Agular spectrum method
tic;
prop_L = 100e-6;
nz = 101;
z_list = linspace(0,prop_L,nz);
% slice through the focal spot(x,y v.s. propgation distance)
E_row = zeros(obj_size(2),nz);
E_col = zeros(obj_size(1),nz);
plot_gif = false;
gif_filename = "TEST.gif";

for i=1:nz
    T = exp(1i*k.*gamma_cust.*(z_list(1,i)));
    E = AS_method(E0,T,F_struct,z_list(1,i),1.55e-6);
    E_row(:,i) = abs(E(4304,:))^2;
    E_col(:,i) = abs(E(:,4018))^2;            
    % plot gif
    if plot_gif==true
        h = figure(2);
        imagesc(xA_s,yA_s,abs(E)^2);
        hold on
        line([xA(1),xA(end)],[yA(end),yA(end)],'Color','black');
        line([xA(end),xA(end)],[yA(end),yA(1)],'Color','black');
        line([xA(end),xA(1)],[yA(1),yA(1)],'Color','black');
        line([xA(1),xA(1)],[yA(1),yA(end)],'Color','black');
        hold off
        title(['z = ',num2str(z_list(i))]);
        % Capture the plot as an image 
        frame = getframe(h); 
        im = frame2im(frame); 
        [imind,cm] = rgb2ind(im,256); 
        % Write to the GIF File 
        if i == 1 
            imwrite(imind,cm,gif_filename,'gif', 'Loopcount',Inf); 
        else 
            imwrite(imind,cm,gif_filename,'gif','WriteMode','append'); 
        end
    end
    
    %find the actual focal spot
    if i==41
        max_val = max(max(abs(E)));
        f_real = z_list(i);
    elseif z_list(i)>40e-6
        tmp_val = max(max(abs(E)));
        if tmp_val>max_val
            max_val = tmp_val;
            f_real = z_list(i);
        end
    end
    
end
%}
toc;
% plot focal plane
T = exp(1i*k.*gamma_cust.*f_real);
E = AS_method(E0,T,F_struct,f_real,1.55e-6);
plot_range_x = 2*nx_start-nx_end:2*nx_end-nx_start;
plot_range_y = 2*ny_start-ny_end:2*ny_end-ny_start;
figure;
imagesc(xA_s(plot_range_x),yA_s(plot_range_y),abs(E(plot_range_x,plot_range_y)));
hold on
% source field area
line([xA(1),xA(end)],[yA(end),yA(end)],'Color','black');
line([xA(end),xA(end)],[yA(end),yA(1)],'Color','black');
line([xA(end),xA(1)],[yA(1),yA(1)],'Color','black');
line([xA(1),xA(1)],[yA(1),yA(end)],'Color','black');
hold off
title(["focal plane_",num2str(f_real)]);
xlabel("x");ylabel("y");

figure;imagesc(z_list,xA_s,E_row);title("rowSlice");ylabel("x");
figure;imagesc(z_list,yA_s,E_col);title("columnSlice");ylabel("y");
% release some space 
%clear E0 T alpha beta gamma_cust


function E = AS_method(Ein, T, F, z, lambda)
A0 = fftshift(fft2(Ein));
%Create band-limiting matrix
ulx = 1/(((2*F.dFx*z)^2+1)^0.5*lambda);
uly = 1/(((2*F.dFy*z)^2+1)^0.5*lambda);
nx_min = find((abs(F.Fx)-ulx)<0,1);
nx_max = find((abs(F.Fx)-ulx)<0,1,'last');
ny_min = find((abs(F.Fy)-uly)<0,1);
ny_max = find((abs(F.Fy)-uly)<0,1,'last');

Band_Limit_Matrix = zeros(size(Ein));
Band_Limit_Matrix(nx_min:nx_max,ny_min:ny_max)=1;
T = T.*Band_Limit_Matrix;
%Calculation
E = ifft2(ifftshift(A0.*T));

end
