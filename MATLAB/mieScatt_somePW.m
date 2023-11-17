%% Initialization
clear
close all

addpath("./Mie_cylinder")

% From Python script load
% kwave, sample_points_x, inc_angles and model parameters
pw_set = load("data/pw_set.mat");
% structure description (cylinder positions and size)
structure = load("data/cylinder_struct.mat");
% list of incident pw indices to excite the structure
sim_config = load("data/somePW_sim_config.mat");

% Create list of scatterer objects
figure(1)
xlim([-3 3])
ylim([-2 2])
axis equal
if sim_config.with_structure
    cyl_list = [];
    for cyl = transpose(structure.cylinders)
       cyl_list=[cyl_list, scatterer(cyl(1) + cyl(2)*1i, cyl(3), 'soft')];
       cyl_list(end).show()
    end
else
    cyl_list=[scatterer(0, 1, 'dielectric', 1)];
end

% Create space points to evaluate the field
[X_var,Y_var] = meshgrid(sim_config.eval_x,sim_config.eval_y);
points_to_evaluate = X_var + 1i*Y_var;

%% Simulation
% Allocate memory for results
e_field = zeros(...
    length(sim_config.eval_y),...
    length(sim_config.eval_x)...
);

% Phased array model simulate plane waves
% by the result of multiple point source
if pw_set.model_select == "points_as_phased_array"
    % Allocate array
    source_fields = zeros(...
        length(sim_config.eval_y),...
        length(sim_config.eval_x),...
        length(pw_set.sources_pos)...
    );
    % Generate field of each antenna array element
    fprintf("Source points simulated:\n")
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
        fprintf(sprintf("\r %3d out of %3d", i, length(pw_set.sources_pos)))
    end
    fprintf("\n")
    fprintf(sprintf("Exciting with %2d testing incident waves...\n", length(sim_config.pw_indices)))
    for i = sim_config.pw_indices + 1
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
    fprintf(sprintf("Simulating %2d testing incident waves...\n", length(sim_config.pw_indices)))
    for i = sim_config.pw_indices + 1
        % setup an incident plane wave
        if pw_set.model_select == "plane_wave"
            inc = plane_wave(pw_set.inc_angles(i), pw_set.kwave);
        elseif pw_set.model_select == "points_along_line"
            inc = point_source(pw_set.sources_pos(i), pw_set.kwave);
        elseif pw_set.model_select == "points_along_circle"
            inc = point_source(pw_set.sources_pos(i), pw_set.kwave);
        end

        % setup the solver using the incident wave
        p = MieSolver(inc);

        % Add the scatterer
        for scatterer = cyl_list
            warning('off')
            p.addScatterer(scatterer)
        end

        % configure for TE transmission conditions
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
writematrix(e_field, sim_config.result_saving_path)

%% Plot
figure(1)
imagesc(sim_config.eval_x,sim_config.eval_y,real(e_field))
xlim([-4.5 4.5])
ylim([-5 1])
axis equal
colorbar
colormap(jet(256));
caxis([-1 1])
set(gca,'YDir','normal')
%print(gcf,sprintf('output/inc120deg_by_model_%s.png', pw_set.model_select),'-dpng','-r200');
