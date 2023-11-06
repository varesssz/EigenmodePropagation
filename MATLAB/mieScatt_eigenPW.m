%% Initialization
clear
close all

addpath("./Mie_cylinder")

% From Python script load
% kwave, sample_points_x, incident_angles_rad
load("data/pw_set.mat")
% structure description (cylinder positions and size)
load("data/cylinder_struct.mat");
% eigen mode excitation vector
load("data/eigenvector.mat")

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
eval_x = linspace(-4.5, 4.5, 241);
eval_y = linspace(-5.0, 1.0, 161);
[X_var,Y_var] = meshgrid(eval_x,eval_y);
points_to_evaluate = X_var + 1i*Y_var;

%% Simulation
% Allocate memory for results
e_field = zeros(size(eval_y, 2), size(eval_x, 2));

% Progress and time estimation in waitbar
f = waitbar(0, 'Starting');
n = size(incident_angles_rad, 2);
start = datetime;

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
    waitbar(i/n, f, ...
        sprintf('Plane Waves: %d out of %d  |  Remaining: ', i, n) + ...
        string((datetime - start)/i*(n-i)));
end
close(f)

%% Plot
figure(2)
imagesc(eval_x,eval_y,abs(e_field))
xlim([-4.5 4.5])
ylim([-5 1])
axis equal
colorbar
caxis([0 15])
set(gca,'YDir','normal')
print(gcf,sprintf('eigenmode_abs_%d.png', eigenmode_number),'-dpng','-r200');

figure(3)
imagesc(eval_x,eval_y,real(e_field))
xlim([-4.5 4.5])
ylim([-5 1])
axis equal
colorbar
colormap(jet(256));
caxis([-15 15])
set(gca,'YDir','normal')
print(gcf,sprintf('eigenmode_real_%d.png', eigenmode_number),'-dpng','-r200');
