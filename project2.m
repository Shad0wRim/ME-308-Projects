%% calculate the values symbolically
syms tau_d tau_d_var delta delta_var
syms m_e I L
vars = [tau_d, delta];
uncs = [tau_d_var, delta_var];

omega_d = 2*sym(pi) / tau_d;
omega_d_std = sqrt(calc_var(omega_d, vars, uncs));

xi = delta / sqrt(4*sym(pi)^2 + delta^2);
xi_std = sqrt(calc_var(xi, vars, uncs));

omega = omega_d / sqrt(1-xi^2);
omega_std = sqrt(calc_var(omega, vars, uncs));

k_e = m_e*omega^2;
k_e_std = sqrt(calc_var(k_e, vars, uncs));

c_c = 2*m_e*omega;
c_c_std = sqrt(calc_var(c_c, vars, uncs));

c = xi*c_c;
c_std = sqrt(calc_var(c_c, vars, uncs));

young = k_e * L^3 / I / 3;
young_std = sqrt(calc_var(young, vars, uncs));


%% ruler geometry
mass = 0.025 * 20/30.48; % mass of full ruler
m_e = mass * (33 / 140); % beam continuous system
thickness=0.003175;      % thickness of ruler at thickest part
width=0.034925;          % full width of ruler
thickWidth=0.01905;      % length of thickest part
edgeWidth=...            % length of angled part
    width - thickWidth;
crossArea=...            % cross sectional area
    (thickWidth + edgeWidth/2)*thickness; 
L=0.20;                  % length of the ruler
volume=crossArea*L;      % volume of ruler
density=mass/volume;     % density of ruler
I = ...                  % second moment of area of cross section
    thickness^3 * (thickWidth + edgeWidth/2) / 12;


%% read the data and find the peaks
data = readmatrix("highframe.csv");
% plot(data(:,1), data(:,2));
selectedData = repmat(data(:,1) >= 0.2 & data(:,1) <= 0.5, 1,2);
data = reshape(data(selectedData), [], 2);
offset = 7.5e-5;
data(:, 2) = data(:, 2) + offset;
peakLocs = repmat(islocalmax(data(:,2)),1,2);
peaks = reshape(data(peakLocs), [], 2);


%% calculate various values
taus = diff(peaks(:,1)); % time differences
tau_d = mean(taus);      % average time diff
tau_d_var = var(taus);
tau_d_std = std(taus);

ratios = peaks(1:end-1, 2) ./ peaks(2:end, 2);
delta = mean(log(ratios)); % delta decrement
delta_var = var(log(ratios));
delta_std = std(log(ratios));

omega_d = subs(omega_d); % damped frequency
omega_d_std = subs(omega_d_std);
xi = subs(xi);           % damping ratio
xi_std = subs(xi_std);
omega = subs(omega);     % natural frequency
omega_std = subs(omega_std);
k_e = subs(k_e);         % equivalent spring constant
k_e_std = subs(k_e_std);
c_c = subs(c_c);         % critical damping coefficient
c_c_std = subs(c_c_std);
c = subs(c);             % damping coefficient
c_std = subs(c_std);
young = subs(young);     % modulus of elasticity
young_std = subs(young_std);
% imprecise because of Iz precision
% accepted value is about 1.3 GPa
% polypropelene


%% plot and print the data
figure(1);
plot(data(:,1), data(:,2), peaks(:,1), peaks(:,2), "r*");

fprintf("------------------------------------\n");
fprintf("tau_d       = %f s\n", tau_d);
fprintf("tau_d_std   = %f s\n\n", tau_d_std);
fprintf("delta       = %f \n", delta);
fprintf("delta_std   = %f \n\n", delta_std);
fprintf("omega_d     = %f rad/s\n", omega_d);
fprintf("omega_d_std = %f rad/s\n\n", omega_d_std);
fprintf("xi          = %f \n", xi);
fprintf("xi_std      = %f \n\n", xi_std);
fprintf("omega       = %f rad/s\n", omega);
fprintf("omega_std   = %f rad/s\n\n", omega_std);
fprintf("k_e         = %f N/m\n", k_e);
fprintf("k_e_std     = %f N/m\n\n", k_e_std);
fprintf("c           = %f Ns/m\n", c);
fprintf("c_std       = %f Ns/m\n\n", c_std);
fprintf("E           = %f Gpa\n", young / 1e9);
fprintf("E_std       = %f Gpa\n", young_std / 1e9);
fprintf("------------------------------------\n");


%% test values with ode solver
% to test that my values make sense
k_e = double(k_e);
c = double(c);

ode = @(t,x) [x(2); -c/m_e*x(2)-k_e/m_e*x(1)];
tspan = [0, 0.3];
x0 = [-4.053956E-4 + offset, 0.0413350214];
opts = odeset(RelTol=1e-6);
[ts, xs] = ode45(ode, tspan, x0, opts);

figure(2);
plot(ts, xs(:, 1));

%% local functions
% calculate the propogated variance
function var = calc_var(func, vars, uncs)
diffs = arrayfun(@(v) diff(func, v), vars);
var = sum(diffs.^2.*uncs);
end