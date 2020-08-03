% parameters
clear all;
lambda0 = 1.55e-6;
k0 = 2*pi/lambda0;
n = 1 ; %free space
lambda = lambda0/n;
k = k0*n;
f = 50e-6;
z0 = 5e-6; % distance between monitor and meta-atoms
%%%%%%%%%%%%%%%%%%%%%%%%%%% Input Field window %%%%%%%%%%%%%%%%%%%%%%%%%%%
% load input E-field file
S=load('monitorTest/f_50_monitor1.mat');
Ey=S.Ey;
Ei_tmp = squeeze(Ey(:,18,1,2))';

% nx and ny have to be odd, so that no artifacts are introduced when calculating the discrete Fourier transform
% If they are even, the propagated beam will be shifted from its center, for example
xA = S.x'; % Array of window points in the x direction
yA = linspace(-(max(xA)-min(xA))/2,max(xA)+min(xA),length(xA)); % Array of window points in the y direction
Ei = zeros(length(yA),length(xA));
for i=1:length(yA)
    Ei(i,:) = Ei_tmp;
end

[XA,YA] = meshgrid(xA,yA);



% Interpolate the non-uniform mesh grids to uniform grids
dx=0.015e-6; % Spacing
dy=0.015e-6; % Spacing
xA_n = (min(xA)):dx:(max(xA));
yA_n = (min(yA)):dy:(max(yA));

if mod(length(xA_n),2)==0
    xA_n(end)=[];
end
if mod(length(yA_n),2)==0
    yA_n(end)=[];
end

[XA_n,YA_n] = meshgrid(xA_n,yA_n);
Ei_n = interp2(XA,YA,Ei,XA_n,YA_n); % new input field
source_size = size(Ei_n);
nx=source_size(2); % Number of grid points in the x direction
ny=source_size(1); % Number of grid points in the y direction 
x_width = max(xA_n)-min(xA_n); % Computation window width in the x direction
y_width = max(yA_n)-min(yA_n); % Computation window width in the y direction
clear XA YA XA_n YA_n;
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Simulation window %%%%%%%%%%%%%%%%%%%%%%%%%%%%

padding_factor_x = 2;
padding_factor_y = 2;
x_width_s = padding_factor_x*x_width;
y_width_s = padding_factor_y*y_width;
xA_s = cat(2,[xA_n(1)-(padding_factor_x-1)/2*x_width:dx:xA_n(1)],xA_n(2:end-1),...
    [xA_n(end):dx:xA_n(end)+(padding_factor_x-1)/2*x_width]);
yA_s = cat(2,[yA_n(1)-(padding_factor_y-1)/2*y_width:dy:yA_n(1)],yA_n(2:end-1),...
    [yA_n(end):dy:yA_n(end)+(padding_factor_y-1)/2*y_width]);
obj_size=size(meshgrid(xA_s,yA_s)); 

% Combine source field into Computational window
E0 = zeros(obj_size);
nx_start = ceil(obj_size(2)/2)-floor(nx/2) ;
nx_end = ceil(obj_size(2)/2)+floor(nx/2) ;
ny_start = ceil(obj_size(1)/2)-floor(ny/2) ;
ny_end = ceil(obj_size(1)/2)+floor(ny/2);
[~,maxcolumn] = max(max(abs(real(Ei_n))));
[~,maxrow]=max(abs(real(Ei_n(:,maxcolumn))));
maxfield = abs(real(Ei_n(maxrow,maxcolumn)));
Ei_phase = atan2(imag(Ei_n),real(Ei_n));
dphase = wrapToPi(2*pi*(sqrt((xA_n-30e-6).^2+(f-z0)^2)-(f-z0))/lambda);
figure;
plot(xA_n,dphase,xA_n,Ei_phase(ceil(source_size(1)/2),:));
%Ei_fake = real(Ei_n)*maxfield./abs(real(Ei_n))...
%    +imag(Ei_n)*1i*maxfield./abs(real(Ei_n));
E0(ny_start:ny_end,nx_start:nx_end)=Ei_n;

% plot input source field
figure;
imagesc(xA,yA,abs(Ei_n).^2);
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
clear FX FY alpha beta;
F_struct = struct('Fx',Fx,'Fy',Fy,'dFx',dFx,'dFy',dFy);

 %% Multiple planes calculation
tic;
prop_L =60e-6;
nz = 101;
z_list = linspace(z0,prop_L,nz);
% slice through the focal spot(x,y v.s. propgation distance)
E_row = zeros(obj_size(2),nz);
%E_col = zeros(obj_size(1),nz);
outfname = "monitor1.txt";
plot_gif = false;
gif_filename1 = "monitor1.gif";
gif_filename2 = "monitor1_slice.gif";
for i=1:nz
    T = exp(1i*k.*gamma_cust.*(z_list(1,i)));
    E = AS_method(E0,T,F_struct,z_list(1,i),1.55e-6);
    E_row(:,i) = abs(E(ceil(obj_size(1)/2),:)).^2;
    E_col(:,i) = abs(E(:,ceil(obj_size(2)/2))).^2;           
    % plot gif
    if plot_gif==true
        h = figure(3);
        imagesc(xA_s(nx_start:nx_end),yA_s(ny_start:ny_end),...
            abs(E(ny_start:ny_end,nx_start:nx_end)).^2);
        truesize(h,[200,1000]);
        colorbar;
        hold on
        line([xA_n(1),xA_n(end)],[yA_n(end),yA_n(end)],'Color','black');
        line([xA_n(end),xA_n(end)],[yA_n(end),yA_n(1)],'Color','black');
        line([xA_n(end),xA_n(1)],[yA_n(1),yA_n(1)],'Color','black');
        line([xA_n(1),xA_n(1)],[yA_n(1),yA_n(end)],'Color','black');
        hold off
        title(['z = ',num2str(z_list(i))]);
        
        w = figure(4);
        plot(xA_s(nx_start:nx_end),E_row(nx_start:nx_end,i));
        title(['z = ',num2str(z_list(i))]);
        % Capture the plot as an image 
        frame1 = getframe(h); 
        im1 = frame2im(frame1); 
        [imind1,cm1] = rgb2ind(im1,256);
        frame2 = getframe(w); 
        im2 = frame2im(frame2); 
        [imind2,cm2] = rgb2ind(im2,256);
        % Write to the GIF File 
        if i == 1 
            imwrite(imind1,cm1,gif_filename1,'gif', 'Loopcount',Inf);
            imwrite(imind2,cm2,gif_filename2,'gif', 'Loopcount',Inf); 
        else 
            imwrite(imind1,cm1,gif_filename1,'gif','WriteMode','append');
            imwrite(imind2,cm2,gif_filename2,'gif','WriteMode','append');
        end
        
    end
    
    %find the actual focal spot
    %{
    if i==18
        max_val = max(max(abs(E)));
        f_real = z_list(i);
    elseif z_list(i)>90e-6
        tmp_val = max(max(abs(E)));
        if tmp_val>max_val
            max_val = tmp_val;
            f_real = z_list(i);
        end
    end
    %}
end
%}
toc;
figure;imagesc(z_list,xA_s(nx_start:nx_end),...
    E_row(nx_start:nx_end,:));title("rowSlice");ylabel("x");
%figure;imagesc(z_list,yA_s(ny_start:ny_end),...
%    E_col(:,ny_start:ny_end));title("columnSlice");ylabel("x");
writematrix(E_row,outfname,"Delimiter",'tab');

%% Single plane calculation
% plot focal plane
f_real = 50e-6;
T = sparse(exp(1i*k.*gamma_cust.*f_real));
E = AS_method(E0,T,F_struct,f_real,1.55e-6);
plot_range_x = nx_start:nx_end;
plot_range_y = ny_start:ny_end;
figure;
imagesc(xA_s(plot_range_x),yA_s(plot_range_y),abs(E(plot_range_y,plot_range_x)).^2);
hold on
% source field area
line([xA(1),xA(end)],[yA(end),yA(end)],'Color','black');
line([xA(end),xA(end)],[yA(end),yA(1)],'Color','black');
line([xA(end),xA(1)],[yA(1),yA(1)],'Color','black');
line([xA(1),xA(1)],[yA(1),yA(end)],'Color','black');
hold off
title(["focal plane at",num2str(f_real/1e-6)]);
xlabel("x");ylabel("y");


%}
% release some space 
%clear E0 T alpha beta gamma_cust

%% Angular spectrum method
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
Band_Limit_Matrix(ny_min:ny_max,nx_min:nx_max)=1;
T = T.*Band_Limit_Matrix;
%Calculation
E = ifft2(ifftshift(A0.*T));

end
