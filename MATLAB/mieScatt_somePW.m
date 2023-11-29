%% Initialization
clear
close all

addpath("./Mie_cylinder")

% From Python script load
% what pw model to use, list of incident pw indices to excite the structure
conf = load("data/config_somePW.mat");
% kwave, sample_points_x, inc_angles and model parameters
pw_set = load(sprintf("data/pw_set_from_%s.mat", conf.pw_set_model));
% structure description (cylinder positions and size)
structure = load("data/cylinder_struct.mat");

% Create list of scatterer objects
figure(1)
xlim([-3 3])
ylim([-2 2])
axis equal
if conf.with_structure
    cyl_list = [];
    for cyl = transpose(structure.cylinders)
       cyl_list=[cyl_list, scatterer(cyl(1) + cyl(2)*1i, cyl(3), 'soft')];
       cyl_list(end).show()
    end
else
    cyl_list=[scatterer(10, 1, 'dielectric', 1)];
end

% Create space points to evaluate the field
[X_var,Y_var] = meshgrid(conf.eval_x,conf.eval_y);
points_to_evaluate = X_var + 1i*Y_var;

%% Simulation
% Allocate memory for results
e_field = zeros(...
    length(conf.eval_y),...
    length(conf.eval_x)...
);
% Progress and time estimation in waitbar
start = datetime;

% Phased array model simulate plane waves
% by the result of multiple point source
if pw_set.model == "points_as_phased_array"
    % Allocate array
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
        % Progress tracker
        fprintf(sprintf('\rArray elements simulated: %3d out of %3d', i, length(pw_set.sources_pos)))
        remaining = length(pw_set.sources_pos) - i;
        fprintf("  |  Time remaining: " + string((datetime - start) / i * remaining))
    end
    fprintf("\n")
    fprintf(sprintf("Exciting with %d testing incident waves...\n", length(conf.pw_indices)))
    for i = conf.pw_indices + 1
        % Element distances from linspace positions
        steps = abs([real(pw_set.sources_pos), 0] - [0, real(pw_set.sources_pos)]);
        element_distance = mean(steps(2:end-1));
        % Phase shift for each element
        phase_shift = pw_set.kwave * element_distance * cos(pw_set.inc_angles(i));
        multiplier_array = -(length(pw_set.sources_pos) - 1)/2 : 1 : (length(pw_set.sources_pos) - 1)/2;
        phases = exp(1i * phase_shift * multiplier_array);
        % Multiply each element's generated field with the phase shift
        % Also add together, creating a resultant field
        for j = 1 : length(phases)
            e_field = e_field + source_fields(:, :, j) * phases(j);
        end
    end
% Other model simulate plane waves
% by simple one-to-one matching to a point sources
else
    % Iterate over the given plane wave directions
    fprintf(sprintf("Simulating %d testing incident waves...\n", length(conf.pw_indices)))
    for i = conf.pw_indices + 1
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
            warning('off')
            p.addScatterer(scatterer)
        end

        % confure for TE transmission conditions
        % (Transverse Electric) = Electric
        % field vector is perpendicular to the plane of propagation
        p.transmissionTE()

        % solve the scattering problem
        p.solve()

        % Calculate total field
        e_field = e_field + p.getTotalField(points_to_evaluate);
    end
end

% Write out results in file
fprintf("Writing results in file...\n")
writematrix(e_field, conf.result_saving_path)

%% Plot
figure(1)
imagesc(conf.eval_x,conf.eval_y,real(e_field))
xlim([-4.5 4.5])
ylim([-5 1])
axis equal
colorbar
colormap(jet(256));
caxis([-1 1])
set(gca,'YDir','normal')
%print(gcf,sprintf('output/inc120deg_by_model_%s.png', pw_set.model),'-dpng','-r200');
