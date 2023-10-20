%% Initialization
clear
close all

addpath("./Mie_cylinder")

% From Python script load
% kwave, sample_points_x, incident_angles_rad
load("data/pw_set.mat")
% structure description (cylinder positions and size)
load("data/cylinder_struct.mat");

sample_points_y = 0;

% Create scatterer object (Soft = PEC)
figure(1)
xlim([-3 3])
ylim([-2 2])
axis equal
cyl_list = [];
for cyl = transpose(clyinders)
    cyl_list=[cyl_list, scatterer(cyl(1) + cyl(2)*1i, cyl(3), 'soft')];
    cyl_list(end).show()
end

% Create space points to evaluate the field
eval_x = sample_points_x;
eval_y = sample_points_y;
[X_var,Y_var] = meshgrid(eval_x,eval_y);
points_to_evaluate = X_var + 1i*Y_var;

%% Simulation
% Allocate memory for results
e_field = zeros(size(eval_x, 2), size(incident_angles_rad, 2));

i = 0;
for dir = incident_angles_rad
    i = i + 1;
    % setup an incident plane wave
    direction_rad = dir; % deg2rad(dir);
    inc = plane_wave(direction_rad,kwave);
    
    % setup the solver using the incident wave
    p = MieSolver(inc);
    
    % Add the scatterer
    for scatterer = cyl_list
        p.addScatterer(scatterer)
    end
    
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
writematrix(e_field, "data/efield_transfer_mat.txt")
