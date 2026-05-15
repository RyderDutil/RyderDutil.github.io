% ME 310
% Final Project - Rail Gun Linear Machine Analysis
% Jackie Lu, Miles Lough, Ryder Dutil
% Due: 5/12/2026

clc; clear; close all;

%% Givens
B = 2;          % Tesla
l = 1;          % m
R = 0.1;        % Ohm
L = 1e-3;       % H (1 mH)
m = 5;          % kg
c = 20;         % N*s/m
K = B * l;      % T*m
Fstall = 48000; % N (from Problem 4)

%% Problem 2

% Denominator: mLs^2 + (mR+cL)s + (cR+K^2)
den_L = [m*L, (m*R + c*L), (c*R + K^2)];

% Transfer functions (from Problem 1)
TF11 = tf([m, c],   den_L);   % i / eS
TF12 = tf([K],      den_L);   % i / Fappl
TF21 = tf([K],      den_L);   % v / eS
TF22 = tf([-L, -R], den_L);   % v / Fappl

% Pole locations
disp('Problem 2 - Pole Locations:')
disp(pole(TF11))

% Time vector
t = linspace(0, 0.5, 5000);

% Part a: eS = 1kV step, F_appl = 0
[i_a, ~] = step(1000 * TF11, t);
[v_a, ~] = step(1000 * TF21, t);

figure('Name', 'Problem 2a');
subplot(2,1,1);
plot(t, i_a, 'b', 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Current (A)');
title('2a: i(t), eS = 1 kV, Fappl = 0');
grid on;

subplot(2,1,2);
plot(t, v_a, 'r', 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Velocity (m/s)');
title('2a: v(t), eS = 1 kV, Fappl = 0');
grid on;

% Part b: eS = 1kV, F_appl = 50 N (headwind - superposition)
% i(t) = TF11*eS + TF12*Fappl
% v(t) = TF21*eS + TF22*Fappl
[i_eS,   ~] = step(1000 * TF11, t);
[i_F50,  ~] = step(50   * TF12, t);
[v_eS,   ~] = step(1000 * TF21, t);
[v_F50,  ~] = step(50   * TF22, t);

i_b = i_eS + i_F50;
v_b = v_eS + v_F50;

figure('Name', 'Problem 2b');
subplot(2,1,1);
plot(t, i_b, 'b', 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Current (A)');
title('2b: i(t), eS = 1 kV, Fappl = 50 N (headwind)');
grid on;

subplot(2,1,2);
plot(t, v_b, 'r', 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Velocity (m/s)');
title('2b: v(t), eS = 1 kV, Fappl = 50 N (headwind)');
grid on;

% Part c: eS = 1kV, F_appl = -50 N (tailwind)
[i_Fn50, ~] = step(-50 * TF12, t);
[v_Fn50, ~] = step(-50 * TF22, t);

i_c = i_eS + i_Fn50;
v_c = v_eS + v_Fn50;

figure('Name', 'Problem 2c');
subplot(2,1,1);
plot(t, i_c, 'b', 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Current (A)');
title('2c: i(t), eS = 1 kV, Fappl = -50 N (tailwind)');
grid on;

subplot(2,1,2);
plot(t, v_c, 'r', 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Velocity (m/s)');
title('2c: v(t), eS = 1 kV, Fappl = -50 N (tailwind)');
grid on;

%% Problem 5

% Denominator without L
den_noL = [m*R, (c*R + K^2)];

% Transfer functions without L
TF11_noL = tf([m, c], den_noL);
TF12_noL = tf([K],    den_noL);
TF21_noL = tf([K],    den_noL);
TF22_noL = tf([-R],   den_noL);

% Pole locations
disp('Problem 5 - Poles with L:')
disp(pole(TF11))
disp('Problem 5 - Single pole without L:')
disp(pole(TF11_noL))

% Part b - Step response without L, eS = 1kV, F_appl = 0
[i_5b, t_out] = step(1000 * TF11_noL, t);
[v_5b, ~]     = step(1000 * TF21_noL, t);

figure('Name', 'Problem 5b');
subplot(2,1,1);
plot(t_out, i_5b, 'b', 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Current (A)');
title('5b: i(t) without L, eS = 1 kV');
grid on;

subplot(2,1,2);
plot(t_out, v_5b, 'r', 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Velocity (m/s)');
title('5b: v(t) without L, eS = 1 kV');
grid on;

% Part c - i(t) with L vs without L
[i_withL, ~] = step(1000 * TF11,    t);
[i_noL,   ~] = step(1000 * TF11_noL, t);

figure('Name', 'Problem 5c');
plot(t, i_withL, 'b',   'LineWidth', 1.5); hold on;
plot(t, i_noL,   'r--', 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Current (A)');
title('5c: i(t) with L vs without L, eS = 1 kV');
legend('i(t) with L = 1 mH', 'i(t) with L = 0');
grid on;

% Part d - v(t) with L vs without L
[v_withL, ~] = step(1000 * TF21,    t);
[v_noL,   ~] = step(1000 * TF21_noL, t);

figure('Name', 'Problem 5d');
plot(t, v_withL, 'b',   'LineWidth', 1.5); hold on;
plot(t, v_noL,   'r--', 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Velocity (m/s)');
title('5d: v(t) with L vs without L, eS = 1 kV');
legend('v(t) with L = 1 mH', 'v(t) with L = 0');
grid on;

%% Problem 6

figure('Name', 'Problem 6 - TF11');
bode(TF11);
title('TF11: i(s) / eS(s)');
grid on;

figure('Name', 'Problem 6 - TF12');
bode(TF12);
title('TF12: i(s) / Fappl(s)');
grid on;

figure('Name', 'Problem 6 - TF21');
bode(TF21);
title('TF21: v(s) / eS(s)');
grid on;

figure('Name', 'Problem 6 - TF22');
bode(TF22);
title('TF22: v(s) / Fappl(s)');
grid on;

figure('Name', 'Problem 6 - All TFs');
bodeplot(TF11, TF12, TF21, TF22);
legend('TF11: i/eS', 'TF12: i/Fappl', 'TF21: v/eS', 'TF22: v/Fappl');
title('Problem 6: Bode Plots - All Transfer Functions');
grid on;

%% Problem 7

% Simulation time
Tsim = 5;

% Build Simulink model
mdl = 'RailGun_P7';
if bdIsLoaded(mdl)
    close_system(mdl, 0);
end
new_system(mdl);
open_system(mdl);

% Blocks
add_block('simulink/Sources/Step', [mdl '/eS'], ...
    'Time', '0', 'Before', '0', 'After', '1000', 'Position', [50 80 80 110]);

add_block('simulink/Sources/Step', [mdl '/Fappl'], ...
    'Time', '0', 'Before', '0', 'After', '50', 'Position', [50 230 80 260]);

add_block('simulink/Continuous/Transfer Fcn', [mdl '/TF11'], ...
    'Numerator', mat2str([m, c]), 'Denominator', mat2str(den_L), ...
    'Position', [160 70 280 120]);

add_block('simulink/Continuous/Transfer Fcn', [mdl '/TF12'], ...
    'Numerator', mat2str([K]), 'Denominator', mat2str(den_L), ...
    'Position', [160 220 280 270]);

add_block('simulink/Continuous/Transfer Fcn', [mdl '/TF21'], ...
    'Numerator', mat2str([K]), 'Denominator', mat2str(den_L), ...
    'Position', [160 330 280 380]);

add_block('simulink/Continuous/Transfer Fcn', [mdl '/TF22'], ...
    'Numerator', mat2str([-L, -R]), 'Denominator', mat2str(den_L), ...
    'Position', [160 430 280 480]);

add_block('simulink/Math Operations/Sum', [mdl '/Sum_i'], ...
    'Inputs', '++', 'Position', [340 130 370 160]);

add_block('simulink/Math Operations/Sum', [mdl '/Sum_v'], ...
    'Inputs', '++', 'Position', [340 380 370 410]);

add_block('simulink/Sinks/Scope', [mdl '/Scope_i'], ...
    'Position', [450 120 490 160]);

add_block('simulink/Sinks/Scope', [mdl '/Scope_v'], ...
    'Position', [450 370 490 410]);

% Connections
add_line(mdl, 'eS/1',    'TF11/1', 'autorouting', 'on');
add_line(mdl, 'eS/1',    'TF21/1', 'autorouting', 'on');
add_line(mdl, 'Fappl/1', 'TF12/1', 'autorouting', 'on');
add_line(mdl, 'Fappl/1', 'TF22/1', 'autorouting', 'on');
add_line(mdl, 'TF11/1',  'Sum_i/1', 'autorouting', 'on');
add_line(mdl, 'TF12/1',  'Sum_i/2', 'autorouting', 'on');
add_line(mdl, 'TF21/1',  'Sum_v/1', 'autorouting', 'on');
add_line(mdl, 'TF22/1',  'Sum_v/2', 'autorouting', 'on');
add_line(mdl, 'Sum_i/1', 'Scope_i/1', 'autorouting', 'on');
add_line(mdl, 'Sum_v/1', 'Scope_v/1', 'autorouting', 'on');

save_system(mdl);
disp('Simulink model built successfully.')

% Part a: eS = 1kV, F_appl = 50 N
disp('Running Part 7a: eS = 1kV, Fappl = 50 N')
set_param([mdl '/Fappl'], 'After', '50', 'Time', '0');
set_param(mdl, 'StopTime', num2str(Tsim));
sim(mdl);

% Part b: eS = 1kV, F_appl = 0 N
disp('Running Part 7b: eS = 1kV, Fappl = 0 N')
set_param([mdl '/Fappl'], 'After', '0', 'Time', '0');
sim(mdl);

% Part c: eS = 1kV, F_appl = F_stall at t = 1s
disp('Running Part 7c: eS = 1kV, Fappl = Fstall at t = 1s')
set_param([mdl '/Fappl'], 'After', num2str(Fstall), 'Time', '1');
sim(mdl);
