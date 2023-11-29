%% Initialization
clear
close all

addpath("./Mie_cylinder")

% From Python script load
% what pw model to use, eigen mode excitation vector
conf = load("data/config_eigenPW.mat");
% kwave, sample_points_x, inc_angles and model parameters
pw_set = load(sprintf("data/pw_set_from_%s.mat", conf.pw_set_model));
% structure description (cylinder positions and size)
structure = load("data/cylinder_struct.mat");

% Create scatterer object (Soft = PEC)
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
    struct_string = "pec";
else
    cyl_list=[scatterer(10, 1, 'dielectric', 1)];
    struct_string = "vacuum";
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
                e_field_i = e_field_i + e_field_j * conf.eigen_vector(i);
            end
            % Add each steered wave to the resultant field
            e_field = e_field + e_field_i;

        % Progress tracker
        fprintf(sprintf('\rWaves steered: %3d out of %3d', i, length(pw_set.inc_angles)))
        remaining = length(pw_set.inc_angles) - i;
        fprintf("  |  Time remaining: " + string((datetime - start) / i * remaining))
    end
    fprintf("\n")

% Other model simulate plane waves
% by simple one-to-one matching to a point sources
else
    % Iterate over the given plane wave directions
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
        e_field_i = p.getTotalField(points_to_evaluate) * conf.eigen_vector(i);
        e_field = e_field + e_field_i;

        % Progress tracker in command line
        fprintf(sprintf('\rWave sources simulated: %3d out of %3d', i, length(pw_set.inc_angles)))
        remaining = length(pw_set.inc_angles) - i;
        fprintf("  |  Time remaining: " + string((datetime - start) / i * remaining))
    end
    fprintf("\n")
end

% Write out results in file
fprintf("Writing results in file...\n")
writematrix(e_field, sprintf(...
    'output/eigenmode_by_%s_e_field_%d_%s.txt',...
    pw_set.model, conf.eigenmode_number, struct_string)...
)

%% Plot
figure(2)
imagesc(conf.eval_x,conf.eval_y,abs(e_field))
xlim([-4.5 4.5])
ylim([-5 1])
axis equal
colorbar
caxis([0 15])
set(gca,'YDir','normal')
%print(gcf,sprintf('output/eigenmode_png_abs_%d_%s.png', conf.eigenmode_number, struct_string),'-dpng','-r200');

figure(3)
imagesc(conf.eval_x,conf.eval_y,real(e_field))
xlim([-4.5 4.5])
ylim([-5 1])
axis equal
colorbar
colormap(jet(256));
caxis([-15 15])
set(gca,'YDir','normal')
%print(gcf,sprintf('output/eigenmode_png_real_%d_%s.png', conf.eigenmode_number, struct_string),'-dpng','-r200');

