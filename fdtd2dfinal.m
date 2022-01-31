function [] = fdtd2dedit(Nx,Ny,dx,dy,Nt,df,shape,r,l,b,NPML,ur, er, nbc, freq, epssrc, musrc,src)
% Usage: fdtd2d(Nx,Ny,dx,dy,Nt,df,r,NPML,ur, er, nbc, freq, epssrc, musrc,src)
% For a quick test, run the following command in the Command Window: fdtd2d(1000,1000,1,1,10,1e6,200,[20 20 20 20],2,-1+ 6i, 1, 1e8, 1, 1,200)
% Function to calculate E field profile for a sphere/rectangle next to a Total-field/Scattered-field (TF/SF) source using the FDTD method in the 2D domain. The function plots a movie of the transient E field as the simulation runs.
% Input arguments:
%                   Nx = width of simulation region
%                   Ny = height of simulation region
%                   dx = number of grids Nx must be split into
%                   dy = number of grids Ny must be split into
%                   Nt = used for computing the number of iterations (time stepps), a minimum value of 10 is needed
%                   df = frequency step (resolution)
%                   shape = 'sphere'/'rectangle'
%                   r = radius of the sphere, if sphere is chosen
%                   l = Length of rectangle, if rectangle is chosen
%                   b = Breadth of rectangle, if rectangle is chosen
%                   NPML = an array containing the thicknesses of the PML; eg. [20,20,20,20]
%                   ur = Relative permeability of the sphere material
%                   er = Relative permittivity of the sphere material
%                   nbc = Refractive index of the background material (1 for air, 1.33 for water, etc.)
%                   freq = Maximum frequency of interest in Hz

if strcmp(shape,'sphere')
    l = "None";
    b = "None";
end

if strcmp(shape,'rectangle')
    r = 'None';
end

Nx = ceil(Nx/dx);
Ny = ceil(Ny/dy);
e0 = 8.85e-12;
Nx2 = 2*Nx;
Ny2 = 2*Ny;
src = ceil(src/dx);
dx2 = dx/2;
dy2 = dy/2;
epsrc = 1;
musrc = 1;
%Build geometry
if strcmp(shape,'sphere')
    xa2 = [0:Nx2-1]*dx2;
    ya2 = [0:Ny2-1]*dy2;
    xa2 = xa2 - mean(xa2);
    ya2 = ya2 - mean(ya2);
    [Y2,X2] = meshgrid(ya2,xa2);
    ER2 = (X2.^2 + Y2.^2)<= r.^2;
    ER2 = epssrc*(1-ER2) + er*ER2;
    UR2 = (X2.^2 + Y2.^2) <= r.^2;
    UR2 = musrc*(1-UR2) + ur*UR2;
    %Exx, Eyy, Ezz
    URxx = UR2(1:2:Nx2,2:2:Ny2);
    URyy = UR2(2:2:Nx2,1:2:Ny2);
    ERzz = ER2(1:2:Nx2,1:2:Ny2);
    imagesc(ERzz);
end
if strcmp(shape,'sphere2') 
    
    xa2 = [0:Nx2-1]*dx2;
    ya2 = [0:Ny2-1]*dy2;
    xa2 = xa2 - mean(xa2) + r + 1e-9;
    ya2 = ya2 - mean(ya2);
    [Y2,X2] = meshgrid(ya2,xa2);
    ER2 = (X2.^2 + Y2.^2) <= r^2;
    UR2 = (X2.^2 + Y2.^2) <= r^2;
    imagesc(ER2);
    xa2 = [0:Nx2-1]*dx2;
    ya2 = [0:Ny2-1]*dy2;
    xa2 = xa2 - mean(xa2) - r - 1e-9;
    ya2 = ya2 - mean(ya2);
    [Y2,X2] = meshgrid(ya2,xa2);
    ER22 = (X2.^2 + Y2.^2)<= r.^2;
    ER2 = ER2 + ER22;
    ER2 = epssrc*(1-ER2) + er*ER2;
    UR22 = (X2.^2 + Y2.^2)<= r.^2;
    UR2 = UR2 + UR22;
    
    UR2 = musrc*(1-UR2) + ur*UR2;
    %Exx, Eyy, Ezz
    URxx = UR2(1:2:Nx2,2:2:Ny2);
    URyy = UR2(2:2:Nx2,1:2:Ny2);
    ERzz = ER2(1:2:Nx2,1:2:Ny2);
    imagesc(real(ERzz));
end
if strcmp(shape,'rectangle')
    
    nx = round(l/dx);
    nx1 = 1 + floor((Nx - nx)/2);
    nx2 = nx1 + nx - 1;
    
    ny = round(b/dy);
    ny1 = 1 + floor((Ny - ny)/2);
    ny2 = ny1 + ny -1;
    ny1 = ny1 - 25 - 1;
    ny2 = ny2 - 25 - 1;
    
    ER2 = zeros(Nx,Ny);
    UR2 = zeros(Nx,Ny);
    ER2(nx1:nx2,ny1:ny2) = 1;
    UR2(nx1:nx2,ny1:ny2) = 1;
    ny1 = ny1 + 26 + 25 + 1;
    ny2 = ny2 + 26 + 25 + 1;
    ER2(nx1:nx2,ny1:ny2) = 1;
    ER2 = epssrc*(1-ER2) + er*ER2;
    UR2(nx1:nx2,ny1:ny2) = 1;
    UR2 = musrc*(1-UR2) + ur*UR2;
    %Exx, Eyy, Ezz
    URxx = UR2;
    URyy = UR2;
    ERzz = ER2; 
    o = nx1;
    clf
    imagesc(real(ER2));
end

imagesc(real(URyy));

%Source stuff
dt = 1/(Nt*freq);
tau = 0.5/freq;
t0 = 10*tau;
steps = 1/(dt*df);
steps = 10000;
t = [0:steps-1]*dt;
Ezsrc = exp(-((t-t0)/tau).^2);
nsrc = sqrt(musrc*epssrc);
c0 = 299792458;
diract = (nsrc*dy/(2*c0)) + dt/2;
Hxsrc = sqrt(epssrc/musrc)*exp(-((t+diract-t0)/tau).^2);
% A = readmatrix("Au.txt");
% omega = 2*pi*299792458./(A(:,1));
% n = A(:,2) + i*A(:,3);
% gamma = (omegapsquared./((1-n)*i*2*pi*299792458./A(:,1))) - (2*pi*299792458./A(:,1));
% gamma = interp1d(A(:,1), gamma, t);
% %plot(linspace(1,steps,steps),Ezsrc)

%PML
sigx = zeros(Nx2,Ny2);
for nx = 1 : 2*NPML(1)
    nx1 = 2*NPML(1) - nx + 1;
    sigx(nx1,:) = (0.5*e0/dt)*(nx/2/NPML(1))^3;
end
for nx = 1 : 2*NPML(2)
    nx1 = Nx2 - 2*NPML(2) + nx;
    sigx(nx1,:) = (0.5*e0/dt)*(nx/2/NPML(2))^3;
end
sigy = zeros(Nx2,Ny2);
for ny = 1 : 2*NPML(3)
    ny1 = 2*NPML(3) - ny + 1;
    sigy(:,ny1) = (0.5*e0/dt)*(ny/2/NPML(3))^3;
end
for ny = 1 : 2*NPML(4)
    ny1 = Ny2 - 2*NPML(4) + ny;
    sigy(:,ny1) = (0.5*e0/dt)*(ny/2/NPML(4))^3;
end

%PML update coefficients
sigHx = sigx(1:2:Nx2,2:2:Ny2);
sigHy = sigy(1:2:Nx2,2:2:Ny2);
mHx0 = (1/dt) + sigHy/(2*e0);
mHx1 = ((1/dt) - sigHy/(2*e0))./mHx0;
mHx2 = - c0./URxx./mHx0;
mHx3 = - (c0*dt/e0) * sigHx./URxx ./ mHx0;
sigHx = sigx(2:2:Nx2,1:2:Ny2);
sigHy = sigy(2:2:Nx2,1:2:Ny2);
mHy0 = (1/dt) + sigHx/(2*e0);
mHy1 = ((1/dt) - sigHx/(2*e0))./mHy0;
mHy2 = - c0./URyy./mHy0;
mHy3 = - (c0*dt/e0) * sigHy./URyy ./ mHy0;
sigDx = sigx(1:2:Nx2,1:2:Ny2);
sigDy = sigy(1:2:Nx2,1:2:Ny2);
mDz0 = (1/dt) + (sigDx + sigDy)/(2*e0) + sigDx.*sigDy*(dt/4/e0^2);
mDz1 = (1/dt) - (sigDx + sigDy)/(2*e0) - sigDx.*sigDy*(dt/4/e0^2);
mDz1 = mDz1 ./ mDz0;
mDz2 = c0./mDz0;
mDz4 = - (dt/e0^2)*sigDx.*sigDy./mDz0;

%Main FDTD loop
CEx = zeros(Nx,Ny);
CEy = zeros(Nx,Ny);
Ez = zeros(Nx,Ny);
Hx = zeros(Nx,Ny);
Hy = zeros(Nx,Ny);
CHz = zeros(Nx,Ny);
Dz = zeros(Nx,Ny);
ICEx = zeros(Nx,Ny);
ICEy = zeros(Nx,Ny);
IDz = zeros(Nx,Ny);
Px = zeros(Nx,Ny);
Py = zeros(Nx,Ny);
Pz = zeros(Nx,Ny);
densx = zeros(Nx,Ny);
densy = zeros(Nx,Ny);
densz = zeros(Nx,Ny);

% %initialise Fourier transforms
% NFREQ = 1000;
% FREQ = linspace(0,freq, NFREQ);
% K = exp(-1i*2*pi*dt*FREQ);
% REF = zeros(1, NFREQ);
% TRN = zeros(1, NFREQ);
% SRC = zeros(1, NFREQ);
% ABS = zeros(1, NFREQ);

% ET = [];
% ER = [];
% ESRC = [];
% EA = [];
axis tight
myVideo = VideoWriter('test1232.gif'); %open video file
myVideo.FrameRate = 30;  %can adjust this, 5 - 10 works well for me
open(myVideo)
dx = dx*1e9;
dy = dy*1e9;
fig = figure('visible','off');
tic
for T = 1:steps
    
    for nx = 1 : Nx
        for ny = 1 : Ny
            Pz(nx,ny) = Pz(nx,ny) + dt*densz(nx,ny);
        end
    end
    
%     for nx = 1 : Nx
%         for ny = 1 : Ny
%             densz(nx,ny) = ((2-gamma(find(t==((T-1)*dt)))*dt)/(2+gamma(find(t==((T-1)*dt)))*dt))*densz(nx.ny) - ...
%                             ((2*((2*pi*299792458/freq)^2)
%     
%     % Find Curl of Ex and Ey
    for nx = 1 : Nx
        for ny = 1 : Ny-1
            CEx(nx,ny) = (Ez(nx,ny+1) - Ez(nx,ny))/dy;
        end
        CEx(nx,Ny) = (0 - Ez(nx, Ny))/dy;
    end
    %Inject source to the curl of E
    for i = 2:Nx
        CEx(i,src-1) = (Ez(i,src) - Ez(i,src-1))/dy - Ezsrc(T)/dy;
    end
    for ny = 1 : Ny
        for nx = 1 : Nx-1
            CEy(nx,ny) = - (Ez(nx+1,ny) - Ez(nx,ny))/dx;
        end
    CEy(Nx,ny) = - (0 - Ez(Nx,ny))/dx;
    end
    
    %Update H integrations
    ICEx = ICEx + CEx;
    ICEy = ICEy + CEy;

    %Update H field
    Hx = mHx1.*Hx + mHx2.*CEx + mHx3.*ICEx;
    Hy = mHy1.*Hy + mHy2.*CEy + mHx3.*ICEy;
    %fprintf(num2str(size(Hx)))
    %Find curl of H
    CHz(1,1) = (Hy(1,1) - 0)/dx - (Hx(1,1) - 0)/dy;
    for nx = 2 : Nx
        CHz(nx,1) = (Hy(nx,1) - Hy(nx-1,1))/dx - (Hx(nx,1) - 0)/dy;
    end
    for ny = 2 : Ny
        CHz(1,ny) = (Hy(1,ny) - 0)/dx - (Hx(1,ny) - Hx(1,ny-1))/dy;
        for nx = 2 : Nx
            CHz(nx,ny) = (Hy(nx,ny) - Hy(nx-1,ny))/dx - (Hx(nx,ny) - Hx(nx,ny-1))/dy;
        end
    end
    
    %Inject source to the curl of H
    for i = 2:Nx
        CHz(i,src) = (Hy(i,src) - Hy(i-1,src))/dx - (Hx(i,src) - Hx(i,src-1))/dy + Hxsrc(T)/dy;
    end
%     for i = 2:Nx
%         CHz(i,src) = (Hy(i,src) - Hy(i-1,src))/dx - (Hx(i,src) - Hx(i,src-1))/dy + Hxsrc(1,T)/dy;
%     end 
    
    %Update D integrations
    IDz = IDz + Dz;
    
    %Update Dz
    Dz = mDz1.*Dz + mDz2.*CHz + mDz4.*IDz;

    %Update Ez
    Ez = mDz1.*Dz;
    
    clf;
    %fig = imagesc(real(Ez),[-0.3 0.3])
    
    draw2d(linspace(1e-9,Nx*1e-9,Nx),linspace(1e-9,Ny*1e-9,Ny),ERzz',real(Ez'),NPML*1e9);
    hold on
    if strcmp(shape,'circle')
        viscircles([Nx/2,Ny/2], r,'Color','k');
    end
    if strcmp(shape,'rectangle')
        rectangle('Position', [o-25-1,o,l,b])
        rectangle('Position', [o+25+1,o,l,b])
    end
    colormap(redblue)
    title("Step " + num2str(T) + " of " + num2str(steps))
    fprintf("Step " + num2str(T) + " of " + num2str(steps) + '\n')
    colorbar
    drawnow
    
    writeVideo(myVideo, getframe)
end
toc
end