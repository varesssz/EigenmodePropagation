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
load("data/cylinder_struct.mat");

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

% Allocate memory
e_fields = zeros(size(eigen_vectors, 2), size(eval_x, 2));

%% Simulation
starting_eigenmode = 750;
% Progress and time estimation in waitbar
f = waitbar(0, 'Starting');
j_max = size(eigen_vectors, 2);
i_max = size(incident_angles_rad, 2);
start = datetime;

% number of eigenmode simulated in this run
k = 0;
% next eigenmode in line to simulate
j = starting_eigenmode;
for eigen_vector = eigen_vectors(:, starting_eigenmode:end)
    % Allocate memory for results
    e_field = zeros(size(eval_y, 2), size(eval_x, 2));
    % next plane wave in line to excite
    i = 1;
    for direction = incident_angles_rad
        % setup an incident plane wave
        inc = plane_wave(direction,kwave);

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
        completed = i + k * i_max;
        remaining = (j_max - starting_eigenmode) * i_max - completed;
        waitbar(j/j_max, f, ...
            sprintf('Eigenmodes: %d out of %d  |  Remaining: ', j, j_max) + ...
            string((datetime - start)/completed*remaining));
        i = i + 1;
    end
    e_fields(j, :) = e_field;
    j = j + 1;
    k = k + 1;
end
close(f)

%% Write out results in file
writematrix(e_fields, "eigenmode_y=-1_5_all.txt")

%% Plot
figure(2)
imagesc(eval_x,0,abs(e_fields))
title(sprintf('All eigenmode response on the structure, sampled at y=%.1f', eval_y))
xlabel('x')
ylabel('eigenmode index (by descending eigen value)')
colorbar
set(gca,'YDir','normal')
savefig('eigenmode_y=-2_all.fig')
saveas(gcf, 'eigenmode_y=-2_all.png')
