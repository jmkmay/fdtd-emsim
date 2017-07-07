function fdtd_v2
%% January 28, 2016: Second iteration of the FDTD algorithm
%% Assumptions: free-space, plane-wave propagating along z
%% Added additive source, absorbing boundary conditions

%% Fundamental constants
epsilon_0 = 8.85e-12;
mu_0 = 4 * pi * 1e-7;
c = 1/sqrt(epsilon_0*mu_0);
eta_0 = sqrt(mu_0/epsilon_0);
%% Units
nanometers = 1e-9;

%% Grid parameters
delta_z = 1 * nanometers;
delta_t = delta_z/(2 * c);
Nz = 2000;
Nt = 7000;
z = [1:Nz]*delta_z;
t = [1:Nt]*delta_t;

%% Variable initialization
Ex(1:Nz) = 0.0;
Hy(1:Nz) = 0.0;
Ex_prev = 0.0;
Ex_prev_prev = 0.0;
Hy_prev = 0.0;
Hy_prev_prev = 0.0;


%% Source
source = 325;
pulse_width = 50*delta_t;
pulse_centre = 200*delta_t;
wavelength = 100*nanometers;
omega = 2 * pi * c / wavelength;


%% FDTD kernel
reduction = 5;
figure
for n = 2:Nt
    
    Ex(2:Nz) = Ex(2:Nz) - delta_t/(epsilon_0*delta_z) * ...
        (Hy(2:Nz) - Hy(1:Nz-1));
    %% Absorbing boundary condition for Ex
    Ex(1) = Ex_prev_prev;
    Ex_prev_prev = Ex_prev;
    Ex_prev = Ex(2);
    
 %    Ex(source) = exp(-0.5*( (t(n)-pulse_centre)/pulse_width)^2 ) + Ex(source); % Gaussian source
	Ex(source) = sin(omega*t(n)) + Ex(source);
    
    Hy(1:Nz-1) = Hy(1:Nz-1) - delta_t/(mu_0*delta_z) * ...
        (Ex(2:Nz) - Ex(1:Nz-1));
    
    %% Absorbing boundary condition for Hy
    Hy(Nz) = Hy_prev_prev;
    Hy_prev_prev = Hy_prev;
    Hy_prev = Hy(Nz-1);
    
%     Hy(source) = 1/eta_0*exp(-0.5*( (t(n)-pulse_centre)/pulse_width)^2 ) + Hy(source); 
    
   
    if mod(n,reduction)==0 %reduces number of times it plots
    subplot(2,1,1), plot(z,Ex), axis([z(1) z(Nz) -1.5 1.5]), ylabel('Ex')
    subplot(2,1,2), plot(z,Hy), axis([z(1) z(Nz) -1.5/eta_0 1.5/eta_0]), ylabel('Hy')
    drawnow
    end
end
    







