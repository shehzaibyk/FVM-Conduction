close all; clear all; clc;

% Two Dimensional Transient Heat Conduction
% Finite Volume Method
% Constant Wall Temperature Boundaries
%     ___________________________
%    |           Top             |
%    |                           |
% Ly | Left                      | Right
%    |           Bottom          |
%    |___________________________|
%                Lx

% PROPERTIES
Lx = 0.2;             % [m]        Length of Plate
Ly = 0.1;             % [m]        Height of Plate
C = 900;              % [J/(kg.K)  Specific Heat
rho = 2700;           % [kg/(m^3)] Density
K = 250;              % [W/(m.K)]  Thermal conductivity
alpha = K/(rho*C);    % [m^2/s]    Thermal diffusivity

% DOMAIN
dx	= 0.005;                 % Node spacing in x
dy	= 0.005;                 % Node spacing in y
nx	= Lx/dx;                 % Number of Control Volumes in x
ny	= Ly/dy;                 % Number of Control Volumes in y
x	= dx/2:dx:Lx+(dx/2);     % Nodal Positions in x
y	= dy/2:dy:Ly+(dy/2);     % Nodal Positions in y
dV  = dx*dy*1;
Ae = dy; Aw = Ae;
As = dx; An = As;

% TIME INITIALIZATION
time = 60;
dt = 0.1;
Cx = alpha*dt/(dx^2);
Cy = alpha*dt/(dy^2);
timesteps = time/dt;

% BOUNDARY CONDITIONS
Tt = 350;   % [degC]    Top Wall Temperature
Tb = 300;   % [degC]    Bottom Wall Temperature
Tl = 300;   % [degC]    Left Wall Temperature
Tr = 300;   % [degC]    Right Wall Temperature

% INITIALIZATION
Nx = nx+1;
Ny = ny+1;
M = zeros(Nx,Ny);
Tini = 300;
T = Tini*(1+M);
error = 1000;
counter = 0;
Told = T;
Tp = Told;

for k = 1:timesteps

while (error>1e-4)

for i = 1:Nx
for j = 1:Ny
    
    aP0 = rho*C*dV/dt;

if (i>1) && (i<Nx) && (j==Ny) %Top wall
    aE = K*Ae/dx;   aW = K*Aw/dx;   aN = 0;   aS = K*As/dy;
    sP = -2*K*An/dy;
    sU = 2*K*An*Tt/dy;
    aP = aP0 + aE + aW + aN + aS -sP;
    H = aE*Told(i+1,j) + aW*Told(i-1,j);
    V = aN*Tt + aS*Told(i,j-1);
    b = sU + (aP0*Tp(i,j));

elseif (i>1) && (i<Nx) && (j==1) %Bottom wall
    aE = K*Ae/dx;   aW = K*Aw/dx;   aN = K*An/dy;   aS = 0;
    sP = -2*K*As/dy;
    sU = 2*K*As*Tb/dy;
    aP = aP0 + aE + aW + aN + aS -sP;
    H = aE*Told(i+1,j) + aW*Told(i-1,j);
    V = aN*Told(i,j+1) + aS*Tb;
    b = sU + (aP0*Tp(i,j));

elseif (i==1) && (j>1) && (j<Ny) % Left wall
    aE = K*Ae/dx;   aW = 0;   aN = K*An/dy;   aS = K*As/dy;
    sP = -2*K*Aw/dx;
    sU = 2*K*Aw*Tl/dx;
    aP = aP0 + aE + aW + aN + aS -sP;
    H = aE*Told(i+1,j) + aW*Tl;
    V = aN*Told(i,j+1) + aS*Told(i,j-1);
    b = sU + (aP0*Tp(i,j));

elseif (i==Nx) && (j>1) && (j<Ny) % Right wall
    aE = 0;   aW = K*Aw/dx;   aN = K*An/dy;   aS = K*As/dy;
    sP = -2*K*Ae/dx;
    sU = 2*K*Ae*Tr/dx;
    aP = aP0 + aE + aW + aN + aS -sP;
    H = aE*Tr + aW*Told(i-1,j);
    V = aN*Told(i,j+1) + aS*Told(i,j-1);
    b = sU + (aP0*Tp(i,j));

elseif (i==1) && (j==1) % Bottom Left Corner
    aE = K*Ae/dx;   aW = 0;   aN = K*An/dy;   aS = 0;
    sP = -2*(K*As/dy + K*Aw/dx);
    sU = 2*(K*As*Tb/dy + K*Aw*Tl/dx);
    aP = aP0 + aE + aW + aN + aS -sP;
    H = aE*Told(i+1,j) + aW*Tl;
    V = aN*Told(i,j+1) + aS*Tb;
    b = sU + (aP0*Tp(i,j));

elseif (i==Nx) && (j==1)    % Bottom Right Corner
    aE = 0;   aW = K*Aw/dx;   aN = K*An/dy;   aS = 0;
    sP = -2*(K*As/dy + K*Ae/dx);
    sU = 2*(K*As*Tb/dy + K*Ae*Tr/dx);
    aP = aP0 + aE + aW + aN + aS -sP;
    H = aE*Tr + aW*Told(i-1,j);
    V = aN*Told(i,j+1) + aS*Tb;
    b = sU + (aP0*Tp(i,j));

elseif (i==1) && (j==Ny)    % Top Left Corner
    aE = K*Ae/dx;   aW = 0;   aN = 0;   aS = K*As/dy;
    sP = -2*(K*An/dy + K*Aw/dx);
    sU = 2*(K*An*Tt/dy + K*Aw*Tl/dx);
    aP = aP0 + aE + aW + aN + aS -sP;
    H = aE*Told(i+1,j) + aW*Tl;
    V = aN*Tt + aS*Told(i,j-1);
    b = sU + (aP0*Tp(i,j));

elseif (i==Nx) && (j==Ny)    % Top Right Corner
    aE = 0;   aW = K*Aw/dx;   aN = 0;   aS = K*As/dy;
    sP = -2*(K*An/dy + K*Ae/dx);
    sU = 2*(K*An*Tt/dy + K*Ae*Tr/dx);
    aP = aP0 + aE + aW + aN + aS -sP;
    H = aE*Tr + aW*Told(i-1,j);
    V = aN*Tt + aS*Told(i,j-1);
    b = sU + (aP0*Tp(i,j));

else
    aE = K*Ae/dx;   aW = K*Aw/dx;   aN = K*An/dy;   aS = K*As/dy;
    sP = 0;         sU = 0;
    aP = aP0 + aE + aW + aN + aS -sP;
    H = aE*Told(i+1,j) + aW*Told(i-1,j);
    V = aN*Told(i,j+1) + aS*Told(i,j-1);
    b = sU + (aP0*Tp(i,j));

end
T(i,j)=(1/aP)*(H + V + b);
end
end
    error = max(max(abs(T - Told))) ;
    Told = T;
    counter = counter + 1;
end
    error = 1000;
    error1 = max(max(abs(T - Tp)));
    Tp = T ;


    subplot(2,1,1);
    semilogy(k,error1,'marker','o','color','b');
    hold on;
    xlabel('number of time steps (nt)','Fontsize',15,'Fontweight','bold','color','k');
    ylabel('error(residue)','Fontsize',15,'Fontweight','bold','color','k');
    title({['2D  Heat Conduction in Transient State.'],['Number of iterations = ',num2str(counter)]});
    pause(0.001)
   
    subplot(2,1,2);
    contourf(x,y,T');
    colorbar;
    caxis([min([Tl Tr Tb Tt]) max([Tl Tr Tb Tt])]);
    colormap(jet)
    pbaspect([Nx/Ny 1 1]) 
    xlabel('X Axis','Fontsize',15,'Fontweight','bold','color','k');
    ylabel('Y Axis','Fontsize',15,'Fontweight','bold','color','k');

    J(k)=T(32,15);
end

    figure(2)
    plot(dt:dt:time,J)
    xlabel('Time (s)')
    ylabel('Temperature (\circC)')
    title('Time-temperature Profile at (x,y)=(0.16,0.075)');
