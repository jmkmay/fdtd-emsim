function FDTD2D
%% 2D FDTD modified for TM polarization
%% ENGR 459 Midterm Question 3
%% March 10, 2016

%% Fundamental constants
epsilon_0 = 8.85e-12;
mu_0 = 4*pi*1e-7;
c = 1/sqrt(epsilon_0*mu_0);
eta_0 = sqrt(mu_0/epsilon_0);

%% Units
nanometers = 1e-9;


%% Grid parameters
delta_x = 4 * nanometers;
delta_z = delta_x;
delta_t = delta_z/(2*c);  %% Courant condition
Nz = 600; 
Nx = 600;
Nt = 5000;

x = [1:Nx] * delta_x;
z = [1:Nz] * delta_z;
t = [1:Nt] * delta_t;

[X Z] = meshgrid(x,z);

%% Source parameters
source_z = 50;
wavelength = 600 * nanometers;
omega = 2 * pi * c / wavelength;
beam_waist = 10 * delta_x;
beam_centre = Nx/2*delta_x;


%% Field initialization - TE polarization
Ex(1:Nx,1:Nz) = 0.0;
Ez(1:Nx,1:Nz) = 0.0;
Hy(1:Nx,1:Nz) = 0.0;

Hy_prev(1:Nx,1:Nz)=0.0;
Hy_prev_prev(1:Nx,1:Nz)=0.0;

Ex_prev(1:Nx,1:Nz)=0.0;
Ex_prev_prev(1:Nx,1:Nz)=0.0;

Ez_prev(1:Nx,1:Nz)=0.0;
Ez_prev_prev(1:Nx,1:Nz)=0.0;


eps_r(1:Nx,1:Nz) = 1.0;


eps_r((X-Nx/2*delta_x).^2+(Z-Nz/2*delta_z).^2<=(Nx/10*delta_x)^2) = 4.0;

figure
surf(eps_r), shading flat


reduction = 15;
figure
%% FDTD kernel for TM Polarization
for n = 1:Nt
   Ex(1:Nx,1:Nz-1) = Ex(1:Nx,1:Nz-1) + delta_t/(epsilon_0*delta_z) * (Hy(1:Nx,2:Nz)  - Hy(1:Nx,1:Nz-1)); 
   
   Ex(:,Nz) = Ex_prev_prev(:,Nz);
   Ex_prev_prev = Ex_prev;
   Ex_prev(:,Nz)=Ex(:,Nz-1);
   
   Ez(1:Nx-1,1:Nz) = Ez(1:Nx-1,1:Nz) - delta_t/(epsilon_0*delta_x) * (Hy(2:Nx,1:Nz)   - Hy(1:Nx-1,1:Nz));

   Ez(Nx,:)= Ez_prev_prev(Nx,:);
   Ez_prev_prev = Ez_prev;
   Ez_prev(Nx,:) = Ez(Nx-1,:);
     
   Hy(2:Nx,2:Nz) = Hy(2:Nx,2:Nz) - delta_t./(mu_0 * eps_r(2:Nx,2:Nz)) .* ( (Ez(2:Nx,2:Nz)   -  Ez(1:Nx-1,2:Nz))/delta_x  -  (Ex(2:Nx,2:Nz)  -  Ex(2:Nx,1:Nz-1))/delta_z);
   
   Hy(1,:)= Hy_prev_prev(1,:);
   Hy_prev_prev = Hy_prev;
   Hy_prev(1,:) = Hy(2,:);
   
   Hy(:,1)= Hy_prev_prev(:,1);
   Hy_prev_prev = Hy_prev;
   Hy_prev(:,1) = Hy(:,2);

   % One way sinusoidal source
   Hy(:,source_z) = (sin(omega*t(n)) *...
       exp( -((x-beam_centre)/beam_waist).^2))' + Hy(:,source_z);
   Ex(:,source_z) = (-sin(omega*t(n)) * ...
       exp( -((x-beam_centre)/beam_waist).^2)/eta_0)' + Ex(:,source_z);
   
   if mod(n,reduction)==0
   surf(X,Z,Hy), shading flat, axis tight, axis equal, view([0 90])
   xlabel('z'), ylabel('x'), title('TM polarized EM wave propagation in free space')
   drawnow
   end
end











































