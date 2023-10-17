clear
close all

addpath("./Mie_cylinder")

% From Python script load
% kwave, sample_points_x, incident_angles_rad
load("data/pw_set_w120_x5605_k801.mat")
cyl = load("data/cylinders_structA.txt");

sample_points_y = 0;

% Select PWs for validation simulation
excitation_pw_angles = [
    incident_angles_rad(470), ...
    incident_angles_rad(401), ...
    incident_angles_rad(94)
];

figure(1)
% Create scatterer object (Soft = PEC)
scatt1 = scatterer(cyl(1,1) + cyl(1,2)*1i, cyl(1,3), 'soft');
scatt2 = scatterer(cyl(2,1) + cyl(2,2)*1i, cyl(2,3), 'soft');
scatt3 = scatterer(cyl(3,1) + cyl(3,2)*1i, cyl(3,3), 'soft');
scatt1.show()
scatt2.show()
scatt3.show()
xlim([-3 3])
ylim([-2 2])
axis equal

%% Simulation
% Allocate memory for results
e_field = zeros( ...
    size(sample_points_y, 2), ...
    size(sample_points_x, 2) ...
);

% Create space points to evaluate the field
eval_x = sample_points_x;
eval_y = sample_points_y;
[X_var,Y_var] = meshgrid(eval_x,eval_y);
points_to_evaluate = X_var + 1i*Y_var;

for dir = excitation_pw_angles
    % setup an incident plane wave
    direction_rad = dir; % deg2rad(dir);
    inc = plane_wave(direction_rad,kwave);
    
    % setup the solver using the incident wave
    p = MieSolver(inc);
    
    % Add the scatterer
    p.addScatterer(scatt1)
    p.addScatterer(scatt2)
    p.addScatterer(scatt3)
    
    % configure for TE transmission conditions
    % (Transverse Electric) = Electric
    % field vector is perpendicular to the plane of propagation
    p.transmissionTE()
    
    % solve the scattering problem
    p.solve()
    
    % Calculate total field
    e_field = e_field + p.getTotalField(points_to_evaluate);
end

% Write out results in file
writematrix(e_field, "data/e_field_y0_w120_x5605_k801_structA.txt")

%% Plot
figure(2)
imagesc(eval_x,eval_y,real(e_field))
xlim([-4 4])
ylim([-4 3])
axis equal
colorbar
set(gca,'YDir','normal')

