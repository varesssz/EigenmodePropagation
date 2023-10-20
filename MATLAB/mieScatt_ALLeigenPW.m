%% Initialization
clear
close all

addpath("./Mie_cylinder")

% From Python script load
% kwave, sample_points_x, incident_angles_rad
load("data/pw_set.mat")
% eigen mode excitation vector
load("data/eigen_vectors.mat")
% structure description (cylinder positions and size)
load("data/cylinders_structB.txt");

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
eval_x = linspace(-4.5, 4.5, 61);
eval_y = -2;
[X_var,Y_var] = meshgrid(eval_x,eval_y);
points_to_evaluate = X_var + 1i*Y_var;

%% Simulation
% Allocate memory
e_fields = zeros(size(eigen_vectors, 2), size(eval_x, 2));

% Progress and time estimation in waitbar
f = waitbar(0, 'Starting');
n = size(eigen_vectors, 2);
start = datetime;

j = 0;
for eigen_vector = eigen_vectors
    j = j + 1;
    % Allocate memory for results
    e_field = zeros(size(eval_y, 2), size(eval_x, 2));
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
        e_field_i = p.getTotalField(points_to_evaluate) * eigen_vector(i);
        e_field = e_field + e_field_i;
        
        % Progress and time estimation in waitbar
        completed = i + (j-1) * size(incident_angles_rad, 2);
        remaining = n * size(incident_angles_rad, 2) - completed;
        waitbar(j/n, f, ...
            sprintf('Eigenmodes: %d out of %d  |  Remaining: ', j, n) + ...
            string((datetime - start)/completed*remaining));
    end
    e_fields(j, :) = e_field;
end
close(f)

%% Plot
figure(2)
imagesc(eval_x,0,abs(e_fields))
title('All eigenmode response on the structure, sampled at y=-2')
xlabel('x')
ylabel('eigenmode number (by descending eigen value)')
colorbar
set(gca,'YDir','normal')
savefig('eigenmode_y=-2_all.fig')
saveas(gcf, 'eigenmode_y=-2_all.png')
