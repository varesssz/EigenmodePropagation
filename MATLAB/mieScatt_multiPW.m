%% Initialization
clear
close all

addpath("./Mie_cylinder")

% From Python script load
% kwave, sample_points_x, inc_angles and model parameters
pw_set = load("data/pw_set.mat");
% structure description (cylinder positions and size)
structure = load("data/cylinder_struct.mat");

% Create list scatterer objects
figure(1)
xlim([-3 3])
ylim([-2 2])
axis equal
cyl_list = [];
for cyl = transpose(structure.cylinders)
    cyl_list=[cyl_list, scatterer(cyl(1) + cyl(2)*1i, cyl(3), 'soft')];
    cyl_list(end).show()
end

% Create space points to evaluate the field
eval_x = pw_set.sample_points_x;
eval_y = 0;
[X_var,Y_var] = meshgrid(eval_x,eval_y);
points_to_evaluate = X_var + 1i*Y_var;

%% Simulation
% Allocate memory for results
e_field = zeros(length(eval_x), length(pw_set.inc_angles));

% Phased array model simulate plane waves
% by the result of multiple point source
if model_select == "point_sources_as_phased_array"
    % Allocate array
    source_fields = zeros(length(eval_x), length(pw_set.sources_pos));
    % Generate field of each antenna array element
    for i = 1 : length(pw_set.sources_pos)
        inc = point_source(pw_set.sources_pos(i), pw_set.kwave);
        p = MieSolver(inc);
        for scatterer = cyl_list
            p.addScatterer(scatterer)
        end
        p.transmissionTE()
        p.solve()
        source_fields(:, i) = p.getTotalField(points_to_evaluate);
    end
    for i = 1 : length(pw_set.inc_angles)
        % Element distances from linspace positions
        steps = abs([real(pw_set.sources_pos), 0] - [0, real(pw_set.sources_pos)]);
        element_distance = mean(steps(2:end-1));
        % Phase shift for each element
        phase_shift = pw_set.kwave * element_distance * cos(pw_set.inc_angles(i));
        phases = exp(1i * phase_shift * (0 : 1 : length(pw_set.sources_pos) - 1));
        % Multiply each element's generated field with the phase shift
        % Also add together, creating a resultant field
        e_field(:, i) = source_fields * phases;
    end
% Other model simulate plane waves
% by simple one-to-one matching to a point sources
else
    % Iterate over every plane wave direction
    for i = 1 : length(pw_set.inc_angles)
        % setup an incident plane wave
        if pw_set.model_select == "plane_wave"
            inc = plane_wave(pw_set.inc_angles(i), pw_set.kwave);
        elseif pw_set.model_select == "point_sources_along_line"
            inc = point_source(pw_set.sources_pos(i), pw_set.kwave);
        elseif pw_set.model_select == "point_sources_along_circle"
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
        e_field(:, i) = p.getTotalField(points_to_evaluate);
    end
end

% Write out results in file
writematrix(e_field, "data/efield_transfer_mat.txt")
