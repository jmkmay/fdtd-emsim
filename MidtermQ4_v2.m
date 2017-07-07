function fdtd_v2
%% February 2, 2016: Third iteration of the FDTD algorithm
%% Assumptions: free-space, plane-wave propagating along z
%% Added additive source, absorbing boundary conditions
%% Added one-way source and simple media

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
Nt = 10000;
z = [1:Nz]*delta_z;
t = [1:Nt]*delta_t;

%% Variable initialization
Ex(1:Nz) = 0.0;
Hy(1:Nz) = 0.0;
er(1:Nz) = 1.0;
mr(1:Nz) = 1.0;
d1 = 200;
d2 = 168;
d3 = 140;
d4 = 58;
d5 = 120;
d6 = 72;
start = 800;
er(start:start+d1) = 1.13^2;
er(start+d1+1:start+d1+d2) = 1.19^2;
er(start+d1+d2+1:start+d1+d2+d3) = 1.25^2;
er(start+d1+d2+d3+1:start+d1+d2+d3+d4) = 1.73^2;
er(start+d1+d2+d3+d4+1:start+d1+d2+d3+d4+d5) = 1.25^2;
er(start+d1+d2+d3+d4+d5+1:start+d1+d2+d3+d4+d5+d6) = 1.73^2;
er(start+d1+d2+d3+d4+d5+d6+1:Nz) = 4.2^2;
reflection = 0;
reflection_ave = 0;
incident = 0;
incident_ave = 0;


Ex_prev = 0.0;
Ex_prev_prev = 0.0;
Hy_prev = 0.0;
Hy_prev_prev = 0.0;


%% Source
source = 325;hj
wavelength = 900*nanometers;
omega = 2 * pi * c / wavelength;


%% FDTD kernel
reduction = 5;
figure
for n = 2:Nt
    
    Ex(2:Nz) = Ex(2:Nz) - delta_t./(er(2:Nz)*epsilon_0*delta_z) .* (Hy(2:Nz) - Hy(1:Nz-1));
    %% Absorbing boundary condition for Ex
    Ex(1) = Ex_prev_prev;
    Ex_prev_prev = Ex_prev;
    Ex_prev = Ex(2);
    
 %    Ex(source) = exp(-0.5*( (t(n)-pulse_centre)/pulse_width)^2 ) + Ex(source); % Gaussian source
	Ex(source) = sin(omega*t(n))/2 + Ex(source);
    Hy(source) = 1/eta_0*sin(omega*t(n))/2 + Hy(source);
    
    Hy(1:Nz-1) = Hy(1:Nz-1) - delta_t./(mr(1:Nz-1)*mu_0*delta_z) .*(Ex(2:Nz) - Ex(1:Nz-1));
    
    %% Absorbing boundary condition for Hy
    Hy(Nz) = Hy_prev_prev;
    Hy_prev_prev = Hy_prev;
    Hy_prev = Hy(Nz-1);
    
    reflection = reflection_ave + abs(Ex(source-50));
    reflection_ave = reflection/n;
    incident = incident + abs(Ex(source+50));
    incident_ave = incident/n;
  
    if mod(n,reduction)==0 %reduces number of times it plots
    subplot(2,1,1), plot(z,Ex), axis([z(1) z(Nz) -1.5 1.5]), ylabel('Ex')
         line([z(start) z(start)], [-1.5 1.5])
         line([z(start+d1) z(start+d1)], [-1.5 1.5])
         line([z(start+d1+d2) z(start+d1+d2)], [-1.5 1.5])
         line([z(start+d1+d2+d3) z(start+d1+d2+d3)], [-1.5 1.5])
         line([z(start+d1+d2+d3+d4) z(start+d1+d2+d3+d4)], [-1.5 1.5])
         line([z(start+d1+d2+d3+d4+d5) z(start+d1+d2+d3+d4+d5)], [-1.5 1.5])
         line([z(start+d1+d2+d3+d4+d5+d6) z(start+d1+d2+d3+d4+d5+d6)], [-1.5 1.5])
        
       
        
%     subplot(2,1,2), plot(z,Hy), axis([z(1) z(Nz) -1.5/eta_0 1.5/eta_0]), ylabel('Hy')
%         line([z(left_edge) z(left_edge)], [-1.5/eta_0 1.5/eta_0])
%         line([z(right_edge) z(right_edge)], [-1.5/eta_0 1.5/eta_0])
    drawnow
    end
end

r_coeff = reflection_ave/incident_ave;
save('rcoeffs','wavelength','r_coeff');
    







