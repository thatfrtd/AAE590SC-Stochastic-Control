%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AAE 590ACA
% HW1 Q2
% Author: Travis Hastreiter 
% Created On: 28 January, 2024
% Description: Simulate 3 body orbits
% Most Recent Change: 28 January, 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Assumed dynamical parameter values
R_E = 6378.1; % [km] Earth radius
mu_E = 398600; % [km3 / s2] Earth gravitational parameter
R_mn = 1738.2; % [km] Moon radius
d_E_mn = 3.8475e5; % [km] Earth-Moon distance
mu_E_mn = 4.0350e5; % [km3 / s2] Earth-Moon barycenter gravitational parameter

% Initial conditions for s/c (non-dimensionalized) 
% Formated as [x, y, z, xdot, ydot, zdot, propagation_time]
IC_1 = [1.2; 0; 0; 0; -1.06110124; 0; 6.20628]; 
IC_2 = [0.85; 0; 0.17546505; 0; 0.2628980369; 0; 2.5543991];
IC_3 = [0.05; -0.05; 0; 4; 2.6; 0; 15];

default_tolerance = 1e-12;

%% a) Simulate orbit under CR3BP
mu_E_mn_mass_ratio = 1.2151e-2; % Moon/Earth mass ratio

d_E_barycenter = d_E_mn * mu_E_mn_mass_ratio;
d_mn_barycenter = d_E_mn - d_E_barycenter;

[p_star] = nondimensionalized_quantities_CR3BP(d_E_barycenter, d_mn_barycenter, [], [], mu_E_mn);

% Simulate
[x_cartesian_N_1, x_star_R_1, t_CR3BP_1] = propagate_orbit_CR3BP(IC_1(1:6), p_star, [0, IC_1(7)], mu_E_mn_mass_ratio, mu_E_mn, d_E_mn, default_tolerance);

[x_cartesian_N_2, x_star_R_2, t_CR3BP_2] = propagate_orbit_CR3BP(IC_2(1:6), p_star, [0, IC_2(7)], mu_E_mn_mass_ratio, mu_E_mn, d_E_mn, default_tolerance);

[x_cartesian_N_3, x_star_R_3, t_CR3BP_3] = propagate_orbit_CR3BP(IC_3(1:6), p_star, [0, IC_3(7)], mu_E_mn_mass_ratio, mu_E_mn, d_E_mn, default_tolerance);

% Assume Moon in circular orbit
mu_mn =  4.9028e3; % [km3 / s2] Moon gravitational parameter
n_mn = sqrt(mu_E_mn / d_E_mn ^ 3); % [rad / s] Moon mean motion

moon_position = @(t) d_E_mn * [cos(n_mn * t); sin(n_mn * t); zeros([1, numel(t)])];

moon_position_CR3BP_1 = moon_position(t_CR3BP_1');
moon_position_CR3BP_2 = moon_position(t_CR3BP_2');
moon_position_CR3BP_3 = moon_position(t_CR3BP_3');

%% a) Plot Synodic Frame
tiledlayout(3,1);

ax1 = nexttile;
earthy(R_E, [-d_E_barycenter; 0; 0], "Earth", 1); hold on;
earthy(R_E, [d_mn_barycenter; 0; 0], "Luna", 1); hold on;
axis equal

plot3(x_star_R_1(:, 1) * p_star(1), x_star_R_1(:, 2) * p_star(1), x_star_R_1(:, 3) * p_star(1), Color = "magenta");
title("IC1")
xlabel("X [km]")
ylabel("Y [km]")
zlabel("Z [km]")
legend("Earth","Moon","Spacecraft")

ax2 = nexttile;
earthy(R_E, [-d_E_barycenter; 0; 0], "Earth", 1); hold on;
earthy(R_E, [d_mn_barycenter; 0; 0], "Luna", 1); hold on;
axis equal

plot3(x_star_R_2(:, 1) * p_star(1), x_star_R_2(:, 2) * p_star(1), x_star_R_2(:, 3) * p_star(1), Color = "magenta");
title("IC2")
xlabel("X [km]")
ylabel("Y [km]")
zlabel("Z [km]")
legend("Earth","Moon","Spacecraft")

ax3 = nexttile;
earthy(R_E, [-d_E_barycenter; 0; 0], "Earth", 1); hold on;
earthy(R_E, [d_mn_barycenter; 0; 0], "Luna", 1); hold on;
axis equal

plot3(x_star_R_3(:, 1) * p_star(1), x_star_R_3(:, 2) * p_star(1), x_star_R_3(:, 3) * p_star(1), Color = "magenta");
title("IC3")
xlabel("X [km]")
ylabel("Y [km]")
zlabel("Z [km]")
legend("Earth","Moon","Spacecraft")

sgtitle(sprintf("CR3BP Orbits in Synodic Frame \n"))


%% a) Plot Inertial Frame
figure
earthy(R_E, [0; 0; 0], "Earth", 1); hold on;
%earthy(R_mn, [d_mn_barycenter; 0; 0], "Luna", 1); hold on;

axis equal

plot_cartesian_orbit(x_cartesian_N_1, "b", 0.3, 1); hold on;
plot_cartesian_orbit(moon_position_CR3BP_1', "r", 0.3, 1)
title("IC1 Inertial Frame CR3BP")

figure
earthy(R_E, [0; 0; 0], "Earth", 1); hold on;
%earthy(R_mn, [d_mn_barycenter; 0; 0], "Luna", 1); hold on;

axis equal

plot_cartesian_orbit(x_cartesian_N_2, "b", 0.3, 1); hold on;
plot_cartesian_orbit(moon_position_CR3BP_2', "r", 0.3, 1)
title("IC2 Inertial Frame CR3BP")

figure
earthy(R_E, [0; 0; 0], "Earth", 1); hold on;
%earthy(R_mn, [d_mn_barycenter; 0; 0], "Luna", 1); hold on;

axis equal

plot_cartesian_orbit(x_cartesian_N_3, "b", 0.3, 1); hold on;
plot_cartesian_orbit(moon_position_CR3BP_3', "r", 0.3, 1)
title("IC3 Inertial Frame CR3BP")

%%
earthy(R_E, [-d_E_barycenter; 0; 0], "Earth", 1); hold on;

axis equal

ax = gca;
ax.NextPlot = 'replaceChildren';
loops = 20;
t_loop = linspace(0, t_CR3BP_1(end), loops);
[~, t_loop_indices] = min(abs(t_loop-t_CR3BP_1));
F(loops) = struct('cdata',[],'colormap',[]);
for j = 1:loops
    
    omega_NR = -sqrt(mu_E_mn / d_E_mn ^ 3); % [rad / s] Earth-Moon mean motion about barycenter
    C_RN = make_R(omega_NR * t_loop(j), 3);
    
    x_cartesian_R_1(:, 1:3) = pagemtimes(C_RN, x_cartesian_N_1(:, 1:3)')';
    moon_position_CR3BP_R_1(1:3, :) = pagemtimes(C_RN, moon_position_CR3BP_1(1:3, :));


    earthy(R_E, C_RN * [-d_E_barycenter; 0; 0], "Earth", 1); hold on;
    earthy(R_E, moon_position_CR3BP_R_1(:, t_loop_indices(j)), "Earth", 1); hold on;
    earthy(R_E * 1.5, x_cartesian_R_1(t_loop_indices(j), :), "Earth", 1); hold on;
    axis equal

    plot_cartesian_orbit(x_cartesian_R_1(1:t_loop_indices(j) , :), "b", 0.3, 1); hold on;
    plot_cartesian_orbit(moon_position_CR3BP_R_1(:, 1:t_loop_indices(j))', "r", 0.3, 1); hold off;

    drawnow
    F(j) = getframe(gcf);
end

fig = figure;
movie(fig,F,1,3)

%% b) Simulate orbit under 3rd body perturbed 2BP
ECI_center_0 = [-d_E_barycenter; 0; 0];

% Convert IC to cartesian
IC_1_2BP = dimensionalize_to_cartesian(IC_1(1:6)', 0, p_star, mu_E_mn, d_E_mn) - [ECI_center_0; zeros([3, 1])];
IC_2_2BP = dimensionalize_to_cartesian(IC_2(1:6)', 0, p_star, mu_E_mn, d_E_mn) - [ECI_center_0; zeros([3, 1])];
IC_3_2BP = dimensionalize_to_cartesian(IC_3(1:6)', 0, p_star, mu_E_mn, d_E_mn) - [ECI_center_0; zeros([3, 1])];

a_d_mn = @(t,x) third_body_perturbation(x(1:3), moon_position(t), mu_mn);

% Simulate 
tolerances = odeset(RelTol=default_tolerance, AbsTol=default_tolerance); 

[t_cartesian_2BP_1,x_cartesian_2BP_1] = ode45(@(t,x) gauss_planetary_eqn(f0_cartesian(x, mu_E), B_cartesian(x, mu_E), a_d_mn(t, x)), [0, IC_1(7) * p_star(3)], IC_1_2BP .* [ones([3, 1]); ones([3, 1])], tolerances);
[t_cartesian_2BP_2,x_cartesian_2BP_2] = ode45(@(t,x) gauss_planetary_eqn(f0_cartesian(x, mu_E), B_cartesian(x, mu_E), a_d_mn(t, x)), [0, IC_2(7) * p_star(3)], IC_2_2BP .* [ones([3, 1]); ones([3, 1])], tolerances);
[t_cartesian_2BP_3,x_cartesian_2BP_3] = ode45(@(t,x) gauss_planetary_eqn(f0_cartesian(x, mu_E), B_cartesian(x, mu_E), a_d_mn(t, x)), [0, IC_3(7) * p_star(3)], IC_3_2BP .* [ones([3, 1]); ones([3, 1])], tolerances);


x_star_2BP_1 = nondimensionalize_from_cartesian_array(x_cartesian_2BP_1, t_cartesian_2BP_1, p_star, mu_E_mn, d_E_mn);
x_star_2BP_2 = nondimensionalize_from_cartesian_array(x_cartesian_2BP_2, t_cartesian_2BP_2, p_star, mu_E_mn, d_E_mn);
x_star_2BP_3 = nondimensionalize_from_cartesian_array(x_cartesian_2BP_3, t_cartesian_2BP_3, p_star, mu_E_mn, d_E_mn);

moon_position_2BP_1 = moon_position(t_cartesian_2BP_1');
moon_position_2BP_2 = moon_position(t_cartesian_2BP_2');
moon_position_2BP_3 = moon_position(t_cartesian_2BP_3');

%% Validate dimensionalize_to_cartesian 
subplot(1, 3, 1)
bar([IC_1_2BP + [ECI_center_0; zeros([3, 1])] x_cartesian_N_1(1, 1:6)'] ./ (p_star(1) * [ones([3, 1]); ones([3, 1]) / p_star(3)]))
title("IC 1")

subplot(1, 3, 2)
bar([IC_2_2BP + [ECI_center_0; zeros([3, 1])] x_cartesian_N_2(1, 1:6)'] ./ (p_star(1) * [ones([3, 1]); ones([3, 1]) / p_star(3)]))
title("IC 2")

subplot(1, 3, 3)
bar([IC_3_2BP + [ECI_center_0; zeros([3, 1])] x_cartesian_N_3(1, 1:6)'] ./ (p_star(1) * [ones([3, 1]); ones([3, 1]) / p_star(3)]))
title("IC 3")

%% Validate nondimensionalize_from_cartesian
IC_1_test = nondimensionalize_from_cartesian(IC_1_2BP' + [ECI_center_0; zeros([3, 1])]', 0, p_star, mu_E_mn, d_E_mn);
IC_2_test = nondimensionalize_from_cartesian(IC_2_2BP' + [ECI_center_0; zeros([3, 1])]', 0, p_star, mu_E_mn, d_E_mn);
IC_3_test = nondimensionalize_from_cartesian(IC_3_2BP' + [ECI_center_0; zeros([3, 1])]', 0, p_star, mu_E_mn, d_E_mn);

subplot(1, 3, 1)
bar([IC_1_test IC_1(1:6)])
title("IC 1")

subplot(1, 3, 2)
bar([IC_2_test IC_2(1:6)])
title("IC 2")

subplot(1, 3, 3)
bar([IC_3_test IC_3(1:6)])
title("IC 3")

%% b) Plot Synodic Frame
% Synodic Frame

tiledlayout(3,1);

ax1 = nexttile;
earthy(R_E, [-d_E_barycenter; 0; 0], "Earth", 1); hold on;
earthy(R_E, [d_mn_barycenter; 0; 0], "Luna", 1); hold on;
axis equal

plot3(x_star_2BP_1(:, 1) * p_star(1), x_star_2BP_1(:, 2) * p_star(1), x_star_2BP_1(:, 3) * p_star(1), Color = "magenta");
title("IC1")
xlabel("X [km]")
ylabel("Y [km]")
zlabel("Z [km]")
legend("Earth","Moon","Spacecraft")

ax2 = nexttile;
earthy(R_E, [-d_E_barycenter; 0; 0], "Earth", 1); hold on;
earthy(R_E, [d_mn_barycenter; 0; 0], "Luna", 1); hold on;
axis equal

plot3(x_star_2BP_2(:, 1) * p_star(1), x_star_2BP_2(:, 2) * p_star(1), x_star_2BP_2(:, 3) * p_star(1), Color = "magenta");
title("IC2")
xlabel("X [km]")
ylabel("Y [km]")
zlabel("Z [km]")
legend("Earth","Moon","Spacecraft")

ax3 = nexttile;
earthy(R_E, [-d_E_barycenter; 0; 0], "Earth", 1); hold on;
earthy(R_E, [d_mn_barycenter; 0; 0], "Luna", 1); hold on;
axis equal

plot3(x_star_2BP_3(:, 1) * p_star(1), x_star_2BP_3(:, 2) * p_star(1), x_star_2BP_3(:, 3) * p_star(1), Color = "magenta");
title("IC3")
xlabel("X [km]")
ylabel("Y [km]")
zlabel("Z [km]")
legend("Earth","Moon","Spacecraft")

sgtitle(sprintf("3rd-Body Perturbed 2BP Orbits in Synodic Frame \n"))


%% b) Plot Inertial Frame

tiledlayout(3,1);

ax1 = nexttile;
earthy(R_E, [0; 0; 0], "Earth", 1); hold on;
axis equal

plot_cartesian_orbit(x_cartesian_2BP_1, "b", 0.3, 1); hold on;
plot_cartesian_orbit(moon_position_2BP_1', "r", 0.3, 1)
title("IC1")
xlabel("X [km]")
ylabel("Y [km]")
zlabel("Z [km]")
legend("Earth","","Spacecraft","Moon")

ax2 = nexttile;
earthy(R_E, [0; 0; 0], "Earth", 1); hold on;
axis equal

plot_cartesian_orbit(x_cartesian_2BP_2, "b", 0.3, 1); hold on;
plot_cartesian_orbit(moon_position_2BP_2', "r", 0.3, 1)
title("IC2")
xlabel("X [km]")
ylabel("Y [km]")
zlabel("Z [km]")
legend("Earth","","Spacecraft","Moon")

ax3 = nexttile;
earthy(R_E, [0; 0; 0], "Earth", 1); hold on;
axis equal

plot_cartesian_orbit(x_cartesian_2BP_3, "b", 0.3, 1); hold on;
plot_cartesian_orbit(moon_position_2BP_3', "r", 0.3, 1)
title("IC3")
xlabel("X [km]")
ylabel("Y [km]")
zlabel("Z [km]")
legend("Earth","","Spacecraft","Moon")

sgtitle(sprintf("3rd-Body Perturbed 2BP Orbits in ECI Frame \n"))

%%

earthy(R_E, [0; 0; 0], "Earth", 1); hold on;

axis equal

ax = gca;
ax.NextPlot = 'replaceChildren';
loops = 20;
t_loop = linspace(0, t_cartesian_2BP_1(end), loops);
[~, t_loop_indices] = min(abs(t_loop-t_cartesian_2BP_1));
F(loops) = struct('cdata',[],'colormap',[]);
for j = 1:loops
    omega_NR = -sqrt(mu_E_mn / d_E_mn ^ 3); % [rad / s] Earth-Moon mean motion about barycenter
    C_RN = make_R(omega_NR * t_loop(j), 3);
    
    x_cartesian_2BP_R_1(:, 1:3) = pagemtimes(C_RN, x_cartesian_2BP_1(:, 1:3)')';
    moon_position_2BP_R_1(1:3, :) = pagemtimes(C_RN, moon_position_2BP_1(1:3, :));


    earthy(R_E, C_RN * [0; 0; 0], "Earth", 1); hold on;
    earthy(R_E, moon_position_2BP_R_1(:, t_loop_indices(j)), "Earth", 1); hold on;
    earthy(R_E * 1.5, x_cartesian_2BP_R_1(t_loop_indices(j), :), "Earth", 1); hold on;
    axis equal

    plot_cartesian_orbit(x_cartesian_2BP_R_1(1:t_loop_indices(j) , :), "b", 0.3, 1); hold on;
    plot_cartesian_orbit(moon_position_2BP_R_1(:, 1:t_loop_indices(j))', "r", 0.3, 1); hold off;


    drawnow
    F(j) = getframe(gcf);
end

fig = figure;
movie(fig,F,1,2)

%% c) Discuss consistency between a) and b)

%% c) Synodic plot

tiledlayout(3,1);

ax1 = nexttile;
earthy(R_E, [-d_E_barycenter; 0; 0], "Earth", 1); hold on;
earthy(R_E, [d_mn_barycenter; 0; 0], "Luna", 1); hold on;
axis equal

plot3(x_star_R_1(:, 1) * p_star(1), x_star_R_1(:, 2) * p_star(1), x_star_R_1(:, 3) * p_star(1), Color = "r"); hold on;
plot3(x_star_2BP_1(:, 1) * p_star(1), x_star_2BP_1(:, 2) * p_star(1), x_star_2BP_1(:, 3) * p_star(1), Color = "b");
title("IC1")
xlabel("X [km]")
ylabel("Y [km]")
zlabel("Z [km]")
legend("Earth","Moon","S/C CR3BP","S/C 2BP")

ax2 = nexttile;
earthy(R_E, [-d_E_barycenter; 0; 0], "Earth", 1); hold on;
earthy(R_E, [d_mn_barycenter; 0; 0], "Luna", 1); hold on;
axis equal

plot3(x_star_R_2(:, 1) * p_star(1), x_star_R_2(:, 2) * p_star(1), x_star_R_2(:, 3) * p_star(1), Color = "r"); hold on;
plot3(x_star_2BP_2(:, 1) * p_star(1), x_star_2BP_2(:, 2) * p_star(1), x_star_2BP_2(:, 3) * p_star(1), Color = "b");
title("IC2")
xlabel("X [km]")
ylabel("Y [km]")
zlabel("Z [km]")
legend("Earth","Moon","S/C CR3BP","S/C 2BP")

ax3 = nexttile;
earthy(R_E, [-d_E_barycenter; 0; 0], "Earth", 1); hold on;
earthy(R_E, [d_mn_barycenter; 0; 0], "Luna", 1); hold on;
axis equal

plot3(x_star_R_3(:, 1) * p_star(1), x_star_R_3(:, 2) * p_star(1), x_star_R_3(:, 3) * p_star(1), Color = "r"); hold on;
plot3(x_star_2BP_3(:, 1) * p_star(1), x_star_2BP_3(:, 2) * p_star(1), x_star_2BP_3(:, 3) * p_star(1), Color = "b");
title("IC3")
xlabel("X [km]")
ylabel("Y [km]")
zlabel("Z [km]")
legend("Earth","Moon","S/C CR3BP","S/C 2BP")

sgtitle(sprintf("CR3BP and 3rd-Body Perturbed 2BP Orbits in Synodic Frame \n"))

%% c) Inertial plot

tiledlayout(3,1);

ax1 = nexttile;
earthy(R_E, [0; 0; 0], "Earth", 1); hold on;
axis equal

plot_cartesian_orbit(moon_position_2BP_1', "g", 0.3, 1); hold on
plot_cartesian_orbit(x_cartesian_N_1 - [ECI_center_0', 0, 0, 0], "r", 0.3, 1); hold on;
plot_cartesian_orbit(x_cartesian_2BP_1, "b", 0.3, 1);
title("IC1")
xlabel("X [km]")
ylabel("Y [km]")
zlabel("Z [km]")
legend("Earth","","Moon","","S/C CR3BP","","S/C 2BP")

ax2 = nexttile;
earthy(R_E, [0; 0; 0], "Earth", 1); hold on;
axis equal

plot_cartesian_orbit(moon_position_2BP_2', "g", 0.3, 1); hold on;
plot_cartesian_orbit(x_cartesian_N_2 - [ECI_center_0', 0, 0, 0], "r", 0.3, 1); hold on;
plot_cartesian_orbit(x_cartesian_2BP_2, "b", 0.3, 1); 
title("IC2")
xlabel("X [km]")
ylabel("Y [km]")
zlabel("Z [km]")
legend("Earth","","Moon","","S/C CR3BP","","S/C 2BP")

ax3 = nexttile;
earthy(R_E, [0; 0; 0], "Earth", 1); hold on;
axis equal

plot_cartesian_orbit(moon_position_2BP_3', "g", 0.3, 1); hold on;
plot_cartesian_orbit(x_cartesian_N_3 - [ECI_center_0', 0, 0, 0], "r", 0.3, 1); hold on;
plot_cartesian_orbit(x_cartesian_2BP_3, "b", 0.3, 1);
title("IC3")
xlabel("X [km]")
ylabel("Y [km]")
zlabel("Z [km]")
legend("Earth","","Moon","","S/C CR3BP","","S/C 2BP")

sgtitle(sprintf("CR3BP and 3rd-Body Perturbed 2BP Orbits in ECI Frame \n"))

%% d) Discuss a perturbation of your choice and add it into a 3BP simulation
C_D = 2.1; % [] drag coefficient
A_over_m = 5.4e-6; % [km2 / kg] spacecraft area to mass ratio
G_0 = 1.02e14; % [kg km / s2] solar flux constant
AU = 149597898; % [km] astronautical unit, Earth-Sun difference
mu_sun = 132712440017.99; % [km3 / s2] Sun gravitational parameter
n_E = sqrt(mu_sun / AU ^ 3); % [rad / s] Earth mean motion

a_d_mn_drag = @(t,x) third_body_perturbation(x(1:3), moon_position(t), mu_mn) ...
                     + SRP_perturbation(t, x(1:3), A_over_m, G_0, AU, n_E);

[t_cartesian_2BP_1_SRP,x_cartesian_2BP_1_SRP] = ode45(@(t,x) gauss_planetary_eqn(f0_cartesian(x, mu_E), B_cartesian(x, mu_E), a_d_mn_drag(t, x)), [0, IC_1(7) * p_star(3)], IC_1_2BP .* [ones([3, 1]); ones([3, 1])], tolerances);
[t_cartesian_2BP_2_SRP,x_cartesian_2BP_2_SRP] = ode45(@(t,x) gauss_planetary_eqn(f0_cartesian(x, mu_E), B_cartesian(x, mu_E), a_d_mn_drag(t, x)), [0, IC_2(7) * p_star(3)], IC_2_2BP .* [ones([3, 1]); ones([3, 1])], tolerances);
[t_cartesian_2BP_3_SRP,x_cartesian_2BP_3_SRP] = ode45(@(t,x) gauss_planetary_eqn(f0_cartesian(x, mu_E), B_cartesian(x, mu_E), a_d_mn_drag(t, x)), [0, IC_3(7) * p_star(3)], IC_3_2BP .* [ones([3, 1]); ones([3, 1])], tolerances);


x_star_2BP_1_SRP = nondimensionalize_from_cartesian_array(x_cartesian_2BP_1_SRP, t_cartesian_2BP_1_SRP, p_star, mu_E_mn, d_E_mn);
x_star_2BP_2_SRP = nondimensionalize_from_cartesian_array(x_cartesian_2BP_2_SRP, t_cartesian_2BP_2_SRP, p_star, mu_E_mn, d_E_mn);
x_star_2BP_3_SRP = nondimensionalize_from_cartesian_array(x_cartesian_2BP_3_SRP, t_cartesian_2BP_3_SRP, p_star, mu_E_mn, d_E_mn);

%% d) Plots

tiledlayout(3,1);

ax1 = nexttile;
earthy(R_E, [-d_E_barycenter; 0; 0], "Earth", 1); hold on;
earthy(R_E, [d_mn_barycenter; 0; 0], "Luna", 1); hold on;
axis equal

plot3(x_star_R_1(:, 1) * p_star(1), x_star_R_1(:, 2) * p_star(1), x_star_R_1(:, 3) * p_star(1), Color = "r"); hold on;
plot3(x_star_2BP_1(:, 1) * p_star(1), x_star_2BP_1(:, 2) * p_star(1), x_star_2BP_1(:, 3) * p_star(1), Color = "b"); hold on;
plot3(x_star_2BP_1_SRP(:, 1) * p_star(1), x_star_2BP_1_SRP(:, 2) * p_star(1), x_star_2BP_1_SRP(:, 3) * p_star(1), Color = "g");
title("IC1")
xlabel("X [km]")
ylabel("Y [km]")
zlabel("Z [km]")
legend("Earth","Moon","S/C CR3BP","S/C 2BP","S/C 2BP + SRP")

ax2 = nexttile;
earthy(R_E, [-d_E_barycenter; 0; 0], "Earth", 1); hold on;
earthy(R_E, [d_mn_barycenter; 0; 0], "Luna", 1); hold on;
axis equal

plot3(x_star_R_2(:, 1) * p_star(1), x_star_R_2(:, 2) * p_star(1), x_star_R_2(:, 3) * p_star(1), Color = "r"); hold on;
plot3(x_star_2BP_2(:, 1) * p_star(1), x_star_2BP_2(:, 2) * p_star(1), x_star_2BP_2(:, 3) * p_star(1), Color = "b"); hold on;
plot3(x_star_2BP_2_SRP(:, 1) * p_star(1), x_star_2BP_2_SRP(:, 2) * p_star(1), x_star_2BP_2_SRP(:, 3) * p_star(1), Color = "g");
title("IC2")
xlabel("X [km]")
ylabel("Y [km]")
zlabel("Z [km]")
legend("Earth","Moon","S/C CR3BP","S/C 2BP","S/C 2BP + SRP")

ax3 = nexttile;
earthy(R_E, [-d_E_barycenter; 0; 0], "Earth", 1); hold on;
earthy(R_E, [d_mn_barycenter; 0; 0], "Luna", 1); hold on;
axis equal

plot3(x_star_R_3(:, 1) * p_star(1), x_star_R_3(:, 2) * p_star(1), x_star_R_3(:, 3) * p_star(1), Color = "r"); hold on;
plot3(x_star_2BP_3(:, 1) * p_star(1), x_star_2BP_3(:, 2) * p_star(1), x_star_2BP_3(:, 3) * p_star(1), Color = "b"); hold on;
plot3(x_star_2BP_3_SRP(:, 1) * p_star(1), x_star_2BP_3_SRP(:, 2) * p_star(1), x_star_2BP_3_SRP(:, 3) * p_star(1), Color = "g");
title("IC3")
xlabel("X [km]")
ylabel("Y [km]")
zlabel("Z [km]")
legend("Earth","Moon","S/C CR3BP","S/C 2BP","S/C 2BP + SRP")

sgtitle(sprintf("CR3BP and 3rd-Body Perturbed 2BP Orbits with and without SRP in Synodic Frame \n"))

%% Helper Functions

function [f_0] = f0_cartesian(x, mu)
    rvec = x(1:3);
    vvec = x(4:6);
    r = norm(rvec);

    a_cartesian = -mu ./ r .^ 3 * rvec;

    f_0 = [vvec; a_cartesian];
end

function [B] = B_cartesian(x, mu)
    B = [zeros(3); eye(3)];
end

function [xdot] = gauss_planetary_eqn(f_0, B, a_d)
    xdot = f_0 + B * a_d;
end

function [a_3B] = third_body_perturbation(rvec_sc, rvec_3B, mu)
    a_3B = -mu * ((rvec_sc - rvec_3B) / norm(rvec_sc - rvec_3B) ^ 3 + rvec_3B / norm(rvec_3B) ^ 3);
end

function [a_SRP] = SRP_perturbation(t, rvec, A_over_m, G_0, a, n)
    
    nu = n * t;
    sun_vector = [cos(nu); sin(nu); 0];

    rvec_E = sun_vector * a;
    rvec_sc = rvec_E + rvec;

    d = norm(rvec_sc);

    a_SRP = A_over_m * G_0 / d ^ 2 * sun_vector;
end

function [a_drag] = drag_perturbation(rvec, vvec, r_o, C_D, A_over_m)
    v = norm(vvec);
    r = norm(rvec);
    alt = r - r_o;

    rho = 1.02e7 * alt ^ -7.172 * 1e9; % Formula (up to 1000 km ish) from 

    a_drag = -1 / 2 * rho * v * vvec * A_over_m * C_D;
end

function [p_star] = nondimensionalized_quantities_CR3BP(r_1, r_2, m_1, m_2, mu_barycenter)
    l_star = r_1 + r_2; % [km]

    G = 6.6743015e-11 / 1000 ^ 2; % [N km2 / kg2]
    if isempty(mu_barycenter)
        m_star = m_1 + m_2; % [kg]
        mu_barycenter = G * m_star;
    else
        m_star = mu_barycenter / G;
    end

    t_star = sqrt(l_star ^ 3 / mu_barycenter); % [s]

    p_star = [l_star, m_star, t_star];
end

function [x_star, tstar] = nondimensionalize_from_cartesian(x_cartesian, t, p_star, mu_barycenter, d)
    l_star = p_star(1);
    t_star = p_star(3);
    
    tstar = t / t_star;

    omega_star_NR = sqrt(mu_barycenter / d ^ 3) * t_star; % [rad / s] Earth-Moon mean motion about barycenter
    C_RN = make_R(omega_star_NR * tstar, 3);

    rvec_star_N = x_cartesian(1:3)' / l_star;
    vvec_star_N = x_cartesian(4:6)' / l_star * t_star;

    rvec_star_R = C_RN * rvec_star_N;
    % Transport theorem
    vvec_star_R = C_RN * (vvec_star_N + cross([0; 0; omega_star_NR], rvec_star_N));

    x_star = [rvec_star_R; vvec_star_R];
end

function [x_star_array, tstar_array] = nondimensionalize_from_cartesian_array(x_cartesian_array, t_array, p_star, mu_barycenter, d)
    x_star_array = zeros(size(x_cartesian_array, 1), 6);
    tstar_array = zeros(size(x_cartesian_array, 1), 1);

    for ind = 1:size(x_cartesian_array, 1)
        x_cartesian = x_cartesian_array(ind, :);

        [x_star_array(ind, :), tstar_array(ind)] = nondimensionalize_from_cartesian(x_cartesian, t_array(ind), p_star, mu_barycenter, d);
    end
end

function [x_cartesian_N, t] = dimensionalize_to_cartesian(x_star, tstar, p_star, mu_barycenter, d)
    l_star = p_star(1);
    t_star = p_star(3);

    t = tstar * t_star;

    omega_RN = -sqrt(mu_barycenter / d ^ 3); % [rad / s] Earth-Moon mean motion about barycenter
    C_NR = make_R(omega_RN * t, 3);

    rvec_cartesian_R = x_star(1:3)' * l_star;
    vvec_cartesian_R = x_star(4:6)' * l_star / t_star;

    rvec_cartesian_N = C_NR * rvec_cartesian_R;
    % Transport theorem
    vvec_cartesian_N = C_NR * (vvec_cartesian_R + cross([0; 0; -omega_RN], rvec_cartesian_R));

    x_cartesian_N = [rvec_cartesian_N; vvec_cartesian_N];
end

function [x_cartesian_N, t_array] = dimensionalize_to_cartesian_array(x_star_array, tstar_array, p_star, mu_barycenter, d)
    x_cartesian_N = zeros(size(x_star_array, 1), 6);
    t_array = zeros(size(x_star_array, 1), 1);

    for ind = 1:size(x_star_array, 1)
        x_star = x_star_array(ind, :);

        [x_cartesian_N(ind, :), t_array(ind)] = dimensionalize_to_cartesian(x_star, tstar_array(ind), p_star, mu_barycenter, d);
    end
end

function [xvecdot] = CR3BP_EoM(x, mu)
    % all quantities are non-dimensionalized
    % mu is mass fraction between two main bodies

    rvec_R = x(1:3);
    vvec_R = x(4:6);

    x = rvec_R(1);
    y = rvec_R(2);
    z = rvec_R(3);

    xdot = vvec_R(1);
    ydot = vvec_R(2);
    zdot = vvec_R(3);

    r_13 = sqrt((x + mu) ^ 2 +  y ^ 2 + z ^ 2); % norm(rvec - rvec_1)
    r_23 = sqrt((x - 1 + mu) ^ 2 +  y ^ 2 + z ^ 2); % norm(rvec - rvec_2)

    xddot = 2 * ydot + x - (1 - mu) * (x + mu) / r_13 ^ 3 - mu * (x - 1 + mu) / r_23 ^ 3;
    yddot = -2 * xdot + y - (1 - mu) * y / r_13 ^ 3 - mu * y / r_23 ^ 3;
    zddot = -(1 - mu) * z / r_13 ^ 3 - mu * z / r_23 ^ 3;

    avec_R = [xddot; yddot; zddot];

    xvecdot = [vvec_R; avec_R];
end

function [x_cartesian_N, x_star_R, t_CR3BP] = propagate_orbit_CR3BP(x0_star, p_star, tspan, mu_mass_ratio, mu_barycenter, d, error_tolerance)
    tolerances = odeset(RelTol=error_tolerance, AbsTol=error_tolerance);

    [t_star_CR3BP,x_star_R] = ode45(@(t,x) CR3BP_EoM(x, mu_mass_ratio), tspan, x0_star, tolerances);
    
    t_star = p_star(3);

    t_CR3BP = t_star_CR3BP * t_star;
    
    [x_cartesian_N] = dimensionalize_to_cartesian_array(x_star_R, t_star_CR3BP, p_star, mu_barycenter, d);
end