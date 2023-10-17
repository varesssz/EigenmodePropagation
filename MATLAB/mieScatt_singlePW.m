clear
close all

addpath("./Mie_cylinder")
addpath("./Mie_cylinder/@MieSolver")

% Wavenumber
c = 299792458;
freq = 1e9; % Hz
kwave = 2*pi*freq/c;
la0 = 2*pi/kwave;

figure(1)
% Create scatterer object
scatt1 = scatterer(-1.0+0i,0.2,'soft'); % Soft = PEC
scatt2 = scatterer(+0.0+0i,0.2,'soft');
scatt3 = scatterer(+0.8+0i,0.2,'soft');
scatt1.show()
scatt2.show()
scatt3.show()
xlim([-4 4])
ylim([-3 3])
axis equal

%% setup an incident plane wave
direction_rad = 2; %deg2rad(90);
inc = plane_wave(direction_rad,kwave);

% setup the solver using the incident wave
p = MieSolver(inc);

% Add the scatterer
p.addScatterer(scatt1)
p.addScatterer(scatt2)
p.addScatterer(scatt3)

% configure for TE transmission conditions (Transverse Electric) = Electric
% field vector is perpendicular to the plane of propagation
p.transmissionTE()

% solve the scattering problem
p.solve()

% Create space points to evaluate the field
eval_x = linspace(-4,4,600);
eval_y = linspace(-3,3,600);
[X_var,Y_var] = meshgrid(eval_x,eval_y);
points_to_evaluate = X_var + 1i*Y_var;

% Calculate total field
E_field = p.getTotalField(points_to_evaluate);

% Plot
figure(2)
imagesc(eval_x,eval_y,real(E_field))
axis equal
colorbar
set(gca,'YDir','normal')

