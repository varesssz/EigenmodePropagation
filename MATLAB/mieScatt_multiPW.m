clear
close all

addpath("./Mie_cylinder")

% From Python script load
% kwave, sample_points_x, incident_angles_rad
load("data/pw_set_w120_x5605_k801.mat")

sample_points_y = 0;

% Cylinder position_x, position_y, radius
cyl = [
    [-1.0, -2.0, 0.2]
    [+0.0, -2.0, 0.2]
    [+0.8, -2.0, 0.2]
];
writematrix(cyl, "data/cylinders_structA.txt")

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
    size(sample_points_x, 2), ...
    size(incident_angles_rad, 2) ...
);
i = 0;

% Create space points to evaluate the field
eval_x = sample_points_x;
eval_y = sample_points_y;
[X_var,Y_var] = meshgrid(eval_x,eval_y);
points_to_evaluate = X_var + 1i*Y_var;

for dir = incident_angles_rad
    i = i + 1;
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
    e_field(:, i) = p.getTotalField(points_to_evaluate);
end

% Write out results in file
writematrix(e_field, "data/efield_transfer_mat_w120_x5605_k801_structA.txt")

% Plot
figure(2)
imagesc(sample_points_x,incident_angles_rad,real(e_field))
axis equal
colorbar
set(gca,'YDir','normal')

