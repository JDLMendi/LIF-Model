clc;

rng(2020) % Setting RNG

% Parameters
t_max = 150e-3;    % Total simulation time (seconds)
dt = 1e-3;        % Time step (seconds)
tau = 20e-3;      % Membrane time constant (seconds)
el = -60e-3;      % Resting potential (volts)
vr = -70e-3;      % Reset potential (volts)
vth = -50e-3;     % Threshold potential (volts)
r = 100e6;        % Membrane resistance (ohms)
i_mean = 25e-11;  % Mean input current (amperes)

max_step_float = t_max/dt;
max_step = floor(max_step_float);

n = 100; % Number of samples per simulation

% Begin plotting figure and hold on
figure;
title("The Leaky-Integrated-And-Fire Model")
xlabel('Time (s)');
ylabel('Membrane Potential (V)');
hold on;

% V = zeros(1, max_step);
tvec = 0:dt:(max_step-1)*dt;
simulations = zeros(n, max_step);

% Loop for multiple simulation
for simulation = 1:n
    V = repmat(el, 1, max_step); % Reset membrane potential for each simulation
    
    for step = 2:max_step
        tvalue = tvec(step);
        
        rand_num = 2 * rand() - 1;
        
        i = i_mean * (1 + 0.1 * (max_step_float) ^ 0.5 * rand_num); % Calculating input current in neuron
        V(step) = V(step-1) + dt/tau * (el - V(step-1)  + r * i); % Calculating Membrane Potential of current step
        
        if V(step) >= vth
            V(step) = vr;
        end
    end
    
    % Plot simulation
    plot(tvec, V, 'black' ,'Marker', 'o', 'LineStyle', 'none', 'MarkerSize', 1, 'HandleVisibility', 'off'); % HandleVisibility helps to prevent it from displaying on legends
    
    simulations(simulation, :) = V;
    
    % % Calculating the mean MP for each step in the simulation
    % v_mean = v_mean + (V - v_mean) / simulation;
    
    % % Plot membrane potential for this simulation
    % plot(tvec, V, 'black' ,'Marker', 'o', 'LineStyle', 'none', 'MarkerSize', 1.2);
    
end

v_mean = mean(simulations); % Calculating mean membrane for each step
v_std = std(simulations); % Calculating Standard Deviation for each step

% Plot mean membrane potential
plot(tvec, v_mean, 'r', 'LineWidth', 1, 'DisplayName', 'Mean')

% Plot standard deviation
plot(tvec, v_mean - v_std, 'black', 'LineWidth', 0.5, 'HandleVisibility', 'off')
plot(tvec, v_mean + v_std, 'black', 'LineWidth', 0.5, 'DisplayName', 'Mean ± STD')

% Plotting the mean sample
% plot(tvec, v_mean, 'b', 'LineWidth', 2);

legend('Location', 'best')

hold off;
