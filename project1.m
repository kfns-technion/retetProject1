%% Question 7--------------------------------------------------------------
syms f(x,l)
f(x,l)= (4*x)- (4*x*l)/(sqrt(x^2+1))-4.905;
x = linspace(-10,10,1000);
l = [0.2,5];

figure(1)
hold on
grid on
plot(x,f(x,l(1)),'LineWidth', 2,'Color', 'r');
plot(x,f(x,l(2)),'LineWidth', 2,'Color', 'g');

title("Example for Question 7")
legend({"l=0.2,function have 1 solution","l=5,function have 3 solutions"})
ax = gca;
ax.XAxisLocation = 'origin'; 
ax.YAxisLocation = 'origin';
ylabel('f(x,l)')
xlabel('x')
grid on;
hold off

a = 1; m = 1; g = 9.81; k = 4; theta = 30*pi/180;

syms V(x,l)
V(x,l)= -m*g*x*sin(theta)+0.5*k*(sqrt(x^2+a^2)-l)^2;
lcr = [0.2,2.5,3.14296,5,7];
x = linspace(-9,11,1000);
figure(2)
hold on
grid on
plot(x,V(x,lcr(1)),'LineWidth', 2,'Color', 'r');
plot(x,V(x,lcr(2)),'LineWidth', 2,'Color', 'g');
plot(x,V(x,lcr(3)),'LineWidth', 2,'Color', 'b');
plot(x,V(x,lcr(4)),'LineWidth', 2,'Color', 'c');
plot(x,V(x,lcr(5)),'LineWidth', 2,'Color', 'y');
ylabel('V(x)')
xlabel('x')
legend('l_0 = 0.2','l_0 = 2.5','l_0 = l_(cr)','l_0 = 5','l_0 = 7')
title('Potential energy V for different l_0 values')

ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';

hold off
%% Question 8--------------------------------------------------------------

a = 1; m = 1; g = 9.81; k = 4; theta = 30*pi/180;
l0_critical = 3.14296;
l0_values = linspace(0, 7, 400);
l0_values(180) = l0_critical;
x_solutions = zeros(3,length(l0_values));

for i = 1:length(l0_values)
    l0 = l0_values(i);
    current_solutions = [];
    xvalues = fzeros(@(x) k*x - (k*x*l0_values(i)) / sqrt(x^2 + a^2) - m*g*sin(theta),-9,11) ;

    % Сохраняем решения
    x_solutions(1:length(xvalues),i) = xvalues;
end

x_solutions = sort_matrix_by_columns(x_solutions);
%% Plot for Q8 (calculation takes a lot of time)---------------------------
% Plots
x_solutions(2,180)=-1.07035;
x_solutions(1,180)=-1.07035;
x_solutions(3,180)=4.28704;
figure;
hold on;
plot(l0_values(180:end), x_solutions(1, 180:end), 'r-', 'LineWidth', 2);
plot(l0_values(180:end), x_solutions(2, 180:end), 'g-', 'LineWidth', 2, 'LineStyle','--');
plot(l0_values, x_solutions(3,:), 'b-', 'LineWidth', 2);

% plot([l0_critical l0_critical], [-5 10], 'k--', 'LineWidth', 1.5);
xline(l0_critical, 'LineWidth',1.5,'Color','black','LineStyle','--')

hold off;
xlabel('l_0 ');
ylabel('x_{eq} ');
title('X_{eq} - l_0');
grid on;
legend('First solution', 'Second soluion', 'Third Solution', 'l_{cr}', 'Location', 'best');

%% Question 9--------------------------------------------------------------

phase_portrait(0.2);
phase_portrait(5);

%% Question 10-------------------------------------------------------------

tspan = [0 100];
x0 = [6.16167 + 2; 0];
options = odeset('RelTol',1e-6,'AbsTol',1e-8);
[t,x] = ode45(@fQ10, tspan, x0, options);

u_lin = exp(-0.05*t).*(2*cos(1.978*t) + 0.05*sin(1.978*t));
x_lin = u_lin + 6.16167; % Calculate x for the linearized system

figure(1);
plot(t, x(:,1), 'r-', 'LineWidth', 1.5);
hold on;
plot(t, x_lin, 'g--', 'LineWidth', 1.5);
hold off;

legend('Nonlinear system', 'Linearized system', 'Location', 'best');
xlabel('Time (t)');
ylabel('x(t)'); % Changed y-label
title('Comparison of Nonlinear and Linearized Systems (x(t))'); % Changed title
grid on;
text(30, max(max(x(:,1)), max(x_lin))*0.9, sprintf('u(0) = %.1f, u''(0) = %.1f', 2, 0),'HorizontalAlignment','left','VerticalAlignment','top');

x0_2 = [6.16167 + 6.7; 0]; % New initial condition

[t_2,x_2] = ode45(@fQ10, tspan, x0_2, options);

% Calculate x for the linearized system
u_lin_2 = exp(-0.05*t_2).*(6.7*cos(1.978*t_2) + 0.168*sin(1.978*t_2));
x_lin_2 = u_lin_2 + 6.16167;

%Plot the results
figure(2);
plot(t_2, x_2(:,1), 'r-', 'LineWidth', 1.5);
hold on;
plot(t_2, x_lin_2, 'g--', 'LineWidth', 1.5);
hold off;

legend('Nonlinear system', 'Linearized system', 'Location', 'best');
xlabel('Time (t)');
ylabel('x(t)');
title('Comparison of Nonlinear and Linearized Systems (x(t))');
grid on;
text(30, max(max(x(:,1)), max(x_lin))*0.5, sprintf('u(0) = %.1f, u''(0) = %.1f', 6.7, 0),'HorizontalAlignment','left','VerticalAlignment','top');
%% Question 11-------------------------------------------------------------

% System parameters
m = 1;      % Mass 
k = 4;      % Stiffness coefficient 
c = 1e-3;   % Nonlinear damping coefficient 

% Initial conditions
x0_1 = 12;  % Initial displacement for the first case 
x0_2 = 50;  % Initial displacement for the second case
v0 = 0;     % Initial velocity

tspan = [0 100]; 

% Function for the right-hand side of the differential equation
fQ11 = @(t, x) [x(2); -c*abs(x(1))*x(2)^3 - k*x(1)/m]; 

% Solve for the first set of initial conditions
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8); % Set relative and absolute tolerances for ode45
[t1, x1] = ode45(fQ11, tspan, [x0_1; v0], options);

% Solve for the second set of initial conditions
[t2, x2] = ode45(fQ11, tspan, [x0_2; v0], options); 

% Analytical envelope equation
omega = sqrt(k/m); % Natural frequency
A0_1 = x0_1; % Initial amplitude for case 1
A0_2 = x0_2; % Initial amplitude for case 2

A_analytic_1 = @(t) A0_1 ./ (1 + (6*c*A0_1^3*omega^4)/(5*pi*k) * t).^(1/3);
A_analytic_2 = @(t) A0_2 ./ (1 + (6*c*A0_2^3*omega^4)/(5*pi*k) * t).^(1/3);

figure; 
% Numerical solution for the first case
subplot(2, 1, 1); 
plot(t1, x1(:, 1), 'b-', 'LineWidth', 1.5); 
hold on;
plot(t1, A_analytic_1(t1), 'r--', 'LineWidth', 1.5); 
plot(t1, -A_analytic_1(t1), 'r--', 'LineWidth', 1.5);
hold off; 
xlabel('Time (t) (s)'); 
ylabel('x(t) (m)');
title(['Numerical Solution and Analytical Envelope, x(0) = ', num2str(x0_1)]); 
legend('Numerical Solution', 'Analytical Envelope');
grid on;
% Numerical solution for the second case
subplot(2, 1, 2); 
plot(t2, x2(:, 1), 'b-', 'LineWidth', 1.5); 
hold on;
plot(t2, A_analytic_2(t2), 'r--', 'LineWidth', 1.5); 
plot(t2, -A_analytic_2(t2), 'r--', 'LineWidth', 1.5); 
hold off; 
xlabel('Time (t) (s)'); 
ylabel('x(t) (m)'); 
title(['Numerical Solution and Analytical Envelope, x(0) = ', num2str(x0_2)]);
legend('Numerical Solution', 'Analytical Envelope'); 
grid on; 

%% Question 12 b-----------------------------------------------------------

% System parameters
m = 1; c1 = 0.1; c3 = 0.001; k = 4; l0 = 5; a = 1; x_eq = 6.16167;
% M0 values (only 0.1 and 10)
M0_values = [0.1, 10];
colors = ['b', 'm']; % Colors for each M0

% Calculate normalized frequency
K_xeq = k - (k*l0)/sqrt(x_eq^2 + a^2);
omega_n = sqrt(K_xeq/m);
% Calculate dimensionless damping coefficients
delta1 = c1/(m*omega_n);
delta3 = (c3*omega_n)/m;

% Frequency range of input moment (Logarithmic spacing)
omega = logspace(-2, 2, 200);

% Preallocate matrices
A_matrix = zeros(length(M0_values), length(omega));
phi_matrix = zeros(length(M0_values), length(omega));

% Calculation of A and phi values
for m_idx = 1:length(M0_values)
    M0 = M0_values(m_idx);
    F = (a*M0)/(K_xeq*(x_eq^2+a^2));

    A_prev = F; % Initial guess

    for i = 1:length(omega)
        Omega_i = omega(i)/omega_n;

        iterations = 0;
        max_iterations = 1000;
        while iterations < max_iterations
            iterations = iterations + 1;
            right_side = 1 / ((-Omega_i^2 + 1)^2 + (delta1*Omega_i + (3/4)*delta3*Omega_i^3*abs(A_prev)^2)^2);
            A_new = sqrt(right_side)*F;
            
            if abs(A_new - A_prev) < 1e-6
                A_matrix(m_idx, i) = A_new;
                % Correct phase calculation using angle
                z = (-Omega_i^2 + 1) + 1i*(delta1*Omega_i + (3/4)*delta3*Omega_i^3*abs(A_new)^2);
                phi_matrix(m_idx, i) = rad2deg(angle(z));
                break;
            end
            A_prev = A_new;
        end
        if iterations == max_iterations
            fprintf('Did not converge for omega = %g, M0 = %g\n', omega(i), M0);
        end
    end
end

% Plotting Bode
figure(1);

% Amplitude plot (top subplot)
subplot(2, 1, 1);
lg_A = loglog(omega, A_matrix(1,:), omega, A_matrix(2,:));
lg_A(1).Color = colors(1);
lg_A(2).Color = colors(2);
lg_A(1).DisplayName = sprintf('M0 = %g', M0_values(1));
lg_A(2).DisplayName = sprintf('M0 = %g', M0_values(2));
lg_A(1).LineWidth = 1.5;
lg_A(2).LineWidth = 1.5;
xlabel('Input Moment Frequency, ω (rad/s)');
ylabel('Amplitude, |A|');
title('Amplitude vs. Input Moment Frequency');
grid on;
legend('Location', 'southwest');

% Phase plot (bottom subplot)
subplot(2, 1, 2);
lg_phi = plot(omega, phi_matrix(1,:),omega, phi_matrix(2,:)); % Use plot for linear scale on Y
set(gca, 'XScale', 'log');
lg_phi(1).Color = colors(1);
lg_phi(2).Color = colors(2);
lg_phi(1).DisplayName = sprintf('M0 = %g', M0_values(1));
lg_phi(2).DisplayName = sprintf('M0 = %g', M0_values(2));
lg_phi(1).LineWidth = 1.5;
lg_phi(2).LineWidth = 1.5;
lg_phi(2).LineStyle="--";
xlabel('Input Moment Frequency, ω (rad/s)');
ylabel('Phase, φ (degrees)');
title('Phase vs. Input Moment Frequency');
grid on;
legend('Location', 'southwest');

%Ploting linear
figure(2)
plot(omega(1:140),A_matrix(1,1:140),'DisplayName',sprintf('M0 = %g', M0_values(1)),'LineWidth',1.5,'Color','r');
xlabel('Input Moment Frequency, ω (rad/s)');
ylabel('Amplitude, |A|');
title('Amplitude vs. Input Moment Frequency');
grid on;
legend('Location',  'northeast');

figure(3)
plot(omega(1:140),A_matrix(2,1:140),'DisplayName',sprintf('M0 = %g', M0_values(2)),'LineWidth',1.5,'Color','m');
xlabel('Input Moment Frequency, ω (rad/s)');
ylabel('Amplitude, |A|');
title('Amplitude vs. Input Moment Frequency');
grid on;
legend('Location',  'northeast');

%% Question 12 c+d-----------------------------------------------------------

% System parameters
m = 1; c1 = 0.1; c3 = 0.001; k = 4; l0 = 5; a = 1; x_eq = 6.16167;
M0_values = [0.1, 10];
omega_multipliers = [0.1, 5];

% Calculate normalized frequency
K_xeq = k - (k*l0)/sqrt(x_eq^2 + a^2);
omega_n = sqrt(K_xeq/m);

% Preallocate arrays for results
Amplitudes = zeros(length(M0_values), length(omega_multipliers));
Phases = zeros(length(M0_values), length(omega_multipliers));

% Counter for figure numbers
figure_counter = 1;

for m_idx = 1:length(M0_values)
    M0 = M0_values(m_idx);
    for w_idx = 1:length(omega_multipliers)
        omega = omega_multipliers(w_idx) * omega_n;

        % Calculate time for 100 cycles
        T = 100 * (2*pi/omega);
        dt = 2*pi/(omega*500); % At least 500 points per cycle
        t = 0:dt:T;

        % System ODE
        odefun = @(t,x) [x(2); (a*M0*cos(omega*t))/(m*(a^2+x_eq^2)) - (c1/m)*x(2) - (c3/m)*x(2)^3 - (k/m)*((sqrt(x_eq^2+a^2)-l0)/sqrt(x_eq^2+a^2))*x(1)];

        % Solve ODE
        [t,x] = ode45(odefun, t, [1; 0]);

        % Analyze last 10 cycles (for amplitude and phase calculation)
        last_cycles_indices = t > t(end) - 10*(2*pi/omega);
        if any(last_cycles_indices)
            first_index = find(last_cycles_indices, 1);
            if first_index > 1
                tvec = t(last_cycles_indices) - t(first_index - 1);
            else
                tvec = t(last_cycles_indices)-t(1);
            end
            x_out = x(last_cycles_indices,1);
            x_in = M0 * cos(omega * tvec);

            % Least squares fit
            Model = [cos(omega*tvec), sin(omega*tvec), ones(size(tvec))];
            consts = Model \ [x_in,x_out];
            A_in = consts(1,1) + 1i*consts(2,1);
            A_out = consts(1,2) + 1i*consts(2,2);
            Amplitudes(m_idx, w_idx) = abs(A_out);
            Phases(m_idx, w_idx) = rad2deg(angle(A_out));

            % Plot ALL cycles in a new figure with desired style
            figure(figure_counter);
            plot(t, x(:, 1), 'b-', 'LineWidth', 1.5); % Blue line, thicker
            title(sprintf('Response: ω = %.2fωₙ, M₀ = %.2f', omega_multipliers(w_idx), M0), 'FontSize', 12); % Formatted title
            xlabel('t', 'FontSize', 12); % Larger font size for labels
            ylabel('x [m]', 'FontSize', 12);
            grid on;
            set(gca, 'FontSize', 12); % Larger font size for axis ticks
            xlim([t(1), t(end)]); % Установка пределов по оси X
            ylim([min(x(:,1))-0.2 max(x(:,1))+0.2]); % Установка пределов по оси Y с небольшим запасом

            % Increment figure counter
            figure_counter = figure_counter + 1;

        else
            fprintf('Warning: less than 10 cycles simulated for M0 = %g, omega = %g. Skipping this point.\n', M0, omega);
            continue;
        end
    end
end

% Display Results in the console
disp('Results:');
for m_idx = 1:length(M0_values)
    for w_idx = 1:length(omega_multipliers)
        fprintf('M0 = %g, ω = %gωₙ: Amplitude = %g, Phase = %g degrees\n', ...
            M0_values(m_idx), omega_multipliers(w_idx), Amplitudes(m_idx, w_idx), Phases(m_idx, w_idx));
    end
end

%% 12 e--------------------------------------------------------------------

% System parameters
m = 1; c1 = 0.1; c3 = 0.001; k = 4; l0 = 5; a = 1; x_eq = 6.16167;
M0_values = [0.1, 10];
% omega_multipliers = logspace(-1, 1, 100);
omega = linspace(0.1, 7, 100);

% Calculate normalized frequency
K_xeq = k - (k*l0)/sqrt(x_eq^2 + a^2);
omega_n = sqrt(K_xeq/m);

% Calculate dimensionless damping coefficients
delta1 = c1/(m*omega_n);
delta3 = (c3*omega_n)/m;

% Preallocate arrays for results
Amplitudes = zeros(length(M0_values), length(omega));
Phases = zeros(length(M0_values), length(omega));

% Counter for figure numbers
figure_counter = 1;

colors = ['b', 'm']; % Colors for each M0

% Preallocate matrices
A_matrix = zeros(length(M0_values), length(omega));
phi_matrix = zeros(length(M0_values), length(omega));

% Calculation of A and phi values
for m_idx = 1:length(M0_values)
    M0 = M0_values(m_idx);
    F = (a*M0)/(K_xeq*(x_eq^2+a^2));

    A_prev = F; % Initial guess

    for i = 1:length(omega)
        Omega_i = omega(i)/omega_n;

        iterations = 0;
        max_iterations = 1000;
        while iterations < max_iterations
            iterations = iterations + 1;
            right_side = 1 / ((-Omega_i^2 + 1)^2 + (delta1*Omega_i + (3/4)*delta3*Omega_i^3*abs(A_prev)^2)^2);
            A_new = sqrt(right_side)*F;
            
            if abs(A_new - A_prev) < 1e-6
                A_matrix(m_idx, i) = A_new;
                % Correct phase calculation using angle
                z = (-Omega_i^2 + 1) + 1i*(delta1*Omega_i + (3/4)*delta3*Omega_i^3*abs(A_new)^2);
                phi_matrix(m_idx, i) = rad2deg(angle(z));
                break;
            end
            A_prev = A_new;
        end
        if iterations == max_iterations
            fprintf('Did not converge for omega = %g, M0 = %g\n', omega(i), M0);
        end
    end
end

for m_idx = 1:length(M0_values)
    M0 = M0_values(m_idx);
    for w_idx = 1:length(omega)
        omega_now = omega(w_idx);

        % Calculate time for 100 cycles
        T = 100 * (2*pi/omega_now);
        dt = 2*pi/(omega_now*500); % At least 500 points per cycle
        t = 0:dt:T;

        % System ODE
        odefun = @(t,x) [x(2); (a*M0*cos(omega_now*t))/(m*(a^2+x_eq^2)) - (c1/m)*x(2) - (c3/m)*x(2)^3 - (k/m)*((sqrt(x_eq^2+a^2)-l0)/sqrt(x_eq^2+a^2))*x(1)];

        % Solve ODE
        [t,x] = ode45(odefun, t, [0; 0]);

        % Analyze last 10 cycles (for amplitude and phase calculation)
        last_cycles_indices = t > t(end) - 10*(2*pi/omega_now);
        if any(last_cycles_indices)
            first_index = find(last_cycles_indices, 1);
            if first_index > 1
                tvec = t(last_cycles_indices) - t(first_index - 1);
            else
                tvec = t(last_cycles_indices)-t(1);
            end
            x_out = x(last_cycles_indices,1);
            x_in = M0 * cos(omega_now * tvec);

            % Least squares fit
            Model = [cos(omega_now*tvec), sin(omega_now*tvec), ones(size(tvec))];
            consts = Model \ [x_in,x_out];
            A_in = consts(1,1) + 1i*consts(2,1);
            A_out = consts(1,2) + 1i*consts(2,2);
            Amplitudes(m_idx, w_idx) = abs(A_out);
            Phases(m_idx, w_idx) = rad2deg(angle(A_out));

        end
    end
end

%Ploting linear
figure(2)
hold on
plot(omega,A_matrix(1,:),'DisplayName',sprintf('analytical for M0 = %g', M0_values(1)),'LineWidth',1.5,'Color','r');
plot(omega,Amplitudes(1,:),'DisplayName',sprintf('numerical for M0 = %g', M0_values(1)),'LineWidth',1.5,'Color','b','LineStyle','--');
xlabel('Input Moment Frequency, ω (rad/s)');
ylabel('Amplitude, |A|');
title('Amplitude vs. Input Moment Frequency');
grid on;
hold off;
legend('Location', 'northeast');

figure(3)
hold on
plot(omega,A_matrix(2,:),'DisplayName',sprintf('analytical for M0 = %g', M0_values(2)),'LineWidth',1.5,'Color','r');
plot(omega,Amplitudes(2,:),'DisplayName',sprintf('numerical for M0 = %g', M0_values(2)),'LineWidth',1.5,'Color','b','LineStyle','--');
xlabel('Input Moment Frequency, ω (rad/s)');
ylabel('Amplitude, |A|');
title('Amplitude vs. Input Moment Frequency');
grid on;
hold off
legend('Location', 'northeast');

%% Functions 
function xval = fzeros(f,xmin,xmax)
% collect all zeros between interval of function
% example:
% fzeros(@(x) (-x^2 - 5*x - 3 + exp(x)),-5,5)
% fzeros(@(x) sin(pi*x) - cos(4*pi*x),0,2)
x = [xmin:0.01:xmax] ;
k = 0 ;
while true
    try
        k = k + 1;
        End = x(k) ;
        A(k) = fzero(f,[x(k) x(k+1)]) ;
    catch
    end
    if End == x(end)
        break
    end
end
if f(0) == 0
    xval = double(unique(single(A))) ;
else
    xval = double(unique(single(A))) ;
    xval(xval == 0) = [] ;
end
end 

function sorted_matrix = sort_matrix_by_columns(matrix)
% Sorts values in a matrix by columns, collecting elements from different rows.
%
% Input:
%   matrix - The input matrix.
%
% Output:
%   sorted_matrix - The matrix with sorted values by columns.

    if ~ismatrix(matrix)
        error('Input argument must be a matrix.');
    end

    num_cols = size(matrix, 2); % Get the number of columns
    num_rows = size(matrix, 1); %Get the number of rows
    sorted_matrix = NaN(size(matrix)); % Initialize the sorted matrix with NaN

    for col = 1:num_cols % Iterate through each column
        values = matrix(:, col); % Extract all values from the current column
        
        %Remove NaN values for correct sorting
        values(isnan(values)) = [];
        
        sorted_values = sort(values); % Sort the extracted values

        %Restore the original column size with NaN
        sorted_col = NaN(num_rows, 1);
        sorted_col(1:length(sorted_values)) = sorted_values;
        
        sorted_matrix(:, col) = sorted_col; % Place the sorted values back into the sorted matrix
    end
end

function phase_portrait(l0)
    % Function for plotting the phase portrait
    M = 0; % External moment
    k = 4;

    function xdot = my_system(t, x)
        xdot = [x(2); M / (x(1)^2 + 1) + 4.905 - 0.1*x(2) - 0.001*(x(2))^3 - 4*x(1) + (k*x(1)*l0) / sqrt(x(1)^2 + 1)];
    end

    figure;
    hold on;

    % Initial conditions for the phase portrait
    x0_values = [ -13, 0, 6];
    v0_values = [-4,7];
    num_initial_conditions = min(length(x0_values) * length(v0_values), 5);

    trajectory_handles = []; % Store handles for legend
    count = 0;
    for i = 1:length(x0_values)
        for j = 1:length(v0_values)
            if count < num_initial_conditions
                x0 = x0_values(i);
                v0 = v0_values(j);
                tspan = [0 200];
                [t, x] = ode45(@my_system, tspan, [x0; v0]);
                h = plot(x(:,1), x(:,2)); % Plot the trajectory and store the handle
                trajectory_handles = [trajectory_handles, h];
                count = count + 1;
            end
        end
    end

    % Finding and stability analysis of equilibrium points with higher accuracy
    options = optimoptions('fsolve','Display','off','TolFun',1e-12,'TolX',1e-12); % Increased accuracy
    f = @(x) [x(2); M / (x(1)^2 + 1) + 4.905 - 0.1*x(2) - 0.001*(x(2))^3 - 4*x(1) + (k*x(1)*l0) / sqrt(x(1)^2 + 1)];

    % Function for calculating the Jacobian matrix
    function J = jacobian_matrix(x)
        J = [0, 1;
            -M*2*x(1)/((x(1)^2+1)^2) - 4 + k*l0*(1/(sqrt(x(1)^2+1)) - x(1)^2/((x(1)^2+1)^(3/2))), -0.1 - 0.003*(x(2)^2)];
    end

    x0_guesses = [[0, 0]; [1, 1]; [-1, -1]; [6,0]; [-4,0]]; % Initial guesses for fsolve 
    equilibrium_points = []; % Store unique equilibrium points
    for i = 1:size(x0_guesses, 1)
        x_eq = fsolve(f, x0_guesses(i,:), options);
        % Check if the point has already been found 
        already_found = false;
        for j = 1:size(equilibrium_points,1)
            if norm(x_eq - equilibrium_points(j,:)) < 1e-10 % higher tolerance for comparing equilibrium points
                already_found = true;
                break;
            end
        end
        if ~already_found
            equilibrium_points = [equilibrium_points; x_eq]; % Add the new equilibrium point
            J = jacobian_matrix(x_eq); % Calculate the Jacobian matrix
            eigenvalues = eig(J); % Calculate the eigenvalues

            if all(real(eigenvalues) < 0) % All eigenvalues have negative real parts
                plot(x_eq(1), x_eq(2), 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 8, 'LineWidth', 2); % Stable point (filled green circle)
            elseif any(real(eigenvalues) > 0) % At least one eigenvalue has a positive real part
                plot(x_eq(1), x_eq(2), 'ro', 'MarkerSize', 8, 'LineWidth', 2); % Unstable point (empty red circle)
            end
        end
    end

    title(['Phase Portrait (l0 = ' num2str(l0) ')']);
    xlabel('x (Position)');
    ylabel('x'' (Velocity)');
    grid on;
    hold off;

    legend_entries = {};
    for k = 1:length(trajectory_handles)
        legend_entries{k} = sprintf('x0 = %g, v0 = %g', x0_values(ceil(k/length(v0_values))), v0_values(mod(k-1,length(v0_values))+1));
    end
    legend(trajectory_handles, legend_entries, 'Location', 'best');

    disp(['Equilibrium points for l0 = ' num2str(l0) ':']);
    disp(equilibrium_points);

end

function xdot = fQ10(t,x)
    xdot = [x(2); -0.1*x(2) - 0.001*x(2)^3 - 4*x(1) + 20*x(1)/sqrt(x(1)^2 + 1) + 4.905];
end

