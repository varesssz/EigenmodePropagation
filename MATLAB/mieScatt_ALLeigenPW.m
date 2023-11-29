%% Initialization
clear
close all

addpath("./Mie_cylinder")

% From Python script load
% what pw model to use, eigenmode excitation vectors, evaluation points and saving path
conf = load("data/config_ALLeigenPW.mat");
% kwave, sample_points_x, inc_angles and model parameters
pw_set = load(sprintf("data/pw_set_from_%s.mat", conf.pw_set_model));
% structure description (cylinder positions and size)
structure = load("data/cylinder_struct.mat");

% Create scatterer object (Soft = PEC)
figure(1)
xlim([-4 4])
ylim([-3 -0.5])
daspect([1 1 1])
cyl_list = [];
for cyl = transpose(structure.cylinders)
    cyl_list=[cyl_list, scatterer(cyl(1) + cyl(2)*1i, cyl(3), 'soft')];
    cyl_list(end).show()
end

% Create space points to evaluate the field
[X_var,Y_var] = meshgrid(conf.eval_x,conf.eval_y);
points_to_evaluate = X_var + 1i*Y_var;
% Plot evaluating points on the scatterer figure
hold on
scatter(real(points_to_evaluate), imag(points_to_evaluate), 6, "filled");
xlabel('x [m]')
ylabel('y [m]')
saveas(gcf, 'output/all_eigenmodes_sampling_and_scatterers.pdf')

%% Simulation
% Progress and time estimation in waitbar
start = datetime;

% Phased array model simulate plane waves
% by the result of multiple point source
if pw_set.model == "points_as_phased_array"
    % Allocate source field array
    source_fields = zeros(...
        length(conf.eval_y),...
        length(conf.eval_x),...
        length(pw_set.sources_pos)...
    );
    % Generate field of each antenna array element
    for i = 1 : length(pw_set.sources_pos)
        inc = point_source(pw_set.sources_pos(i), pw_set.kwave);
        p = MieSolver(inc);
        for scatterer = cyl_list
            warning('off')
            p.addScatterer(scatterer)
        end
        p.transmissionTE()
        p.solve()
        source_fields(:, :, i) = p.getTotalField(points_to_evaluate) * pw_set.sources_w(i);
        
        % Progress and time estimation in waitbar
        fprintf(sprintf('\rArray elements simulated: %3d out of %3d', i, length(pw_set.sources_pos)))
        remaining = length(pw_set.sources_pos) - i;
        fprintf("  |  Time remaining: " + string((datetime - start) / i * remaining))
    end

% Other model simulate plane waves
% by simple one-to-one matching to a point sources
else
    % Allocate source field array
    source_fields = zeros(...
        length(conf.eval_y),...
        length(conf.eval_x),...
        length(pw_set.inc_angles)...
    );
    % next plane wave in line to excite
    for i = 1 : length(pw_set.inc_angles)
        % setup an incident plane wave
        if pw_set.model == "plane_wave"
            inc = plane_wave(pw_set.inc_angles(i), pw_set.kwave);
        elseif pw_set.model == "points_along_line"
            inc = point_source(pw_set.sources_pos(i), pw_set.kwave);
        elseif pw_set.model == "points_along_circle"
            inc = point_source(pw_set.sources_pos(i), pw_set.kwave);
        end

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
        source_fields(:, :, i) = p.getTotalField(points_to_evaluate);
        
        % Progress and time estimation in waitbar

        % Progress and time estimation in waitbar
        fprintf(sprintf('\rWave sources simulated: %3d out of %3d', i, length(pw_set.inc_angles)))
        remaining = length(pw_set.inc_angles) - i;
        fprintf("  |  Time remaining: " + string((datetime - start) / i * remaining))
    end
end
fprintf("\n")

%% Eigenmodes calculation
% Allocate memory for each eigenmode result to store row-by-row
eigenmode_e_fields = zeros(size(conf.eig_vectors, 2), length(conf.eval_x));

% Progress and time estimation in waitbar
start = datetime;

for k = 1 : size(conf.eig_vectors, 2)
    % Allocate memory for one eigenmode results
    e_field = zeros(length(conf.eval_y), length(conf.eval_x));

    % Phased array model simulate plane waves
    % by the result of multiple point source
    if pw_set.model == "points_as_phased_array"
        for i = 1 : length(pw_set.inc_angles)
            e_field_i = zeros(length(conf.eval_y), length(conf.eval_x));

            % Element distances from linspace positions
            steps = abs([real(pw_set.sources_pos), 0] - [0, real(pw_set.sources_pos)]);
            element_distance = mean(steps(2:end-1));

            % Phase shift for each element
            phase_shift = pw_set.kwave * element_distance * cos(pw_set.inc_angles(i));
            multiplier_array = -(length(pw_set.sources_pos) - 1)/2 : 1 : (length(pw_set.sources_pos) - 1)/2;
            phases = exp(1i * phase_shift * multiplier_array);

            for j = 1 : length(phases)
                % Multiply each element's generated field with the phase shift
                e_field_j = source_fields(:, :, j) * phases(j);
                % Multiply with eigenvector component add accumulate to create the steered wave
                e_field_i = e_field_i + e_field_j * conf.eig_vectors(i, k);
            end
            % Add each steered wave to the resultant field
            e_field = e_field + e_field_i;
        end
        
    % Other model simulate plane waves
    % by simple one-to-one matching to a point sources
    else
        % next plane wave in line to excite
        for i = 1 : length(pw_set.inc_angles)
            e_field_i = source_fields(:, :, i) * conf.eig_vectors(i, k);
            e_field = e_field + e_field_i;
        end
    end

    % Save each eigenmode result as a row of a matrix
    eigenmode_e_fields(k, :) = e_field;

    % Progress and time estimation in waitbar
    fprintf(sprintf('\rEigenmodes excited: %3d out of %3d', k, size(conf.eig_vectors, 2)))
    remaining = size(conf.eig_vectors, 2) - k;
    fprintf("  |  Time remaining: " + string((datetime - start) / k * remaining))
end
fprintf("\n")

%% Write out results in file
fprintf("Writing results in file...\n")
writematrix(eigenmode_e_fields, sprintf(conf.result_saving_path, ".txt"))

%% Plot
figure(2)
imagesc(conf.eval_x,0,abs(eigenmode_e_fields))
title(sprintf(...
    'All eigenmode response on the structure, sampled at y=%.2f',...
    conf.eval_y...
))
xlabel('x')
ylabel('eigenmode index (by descending eigen value)')
colorbar
set(gca,'YDir','normal')
% savefig(sprintf(conf.result_saving_path, ".fig"))
% saveas(gcf, sprintf(conf.result_saving_path, ".png"))
