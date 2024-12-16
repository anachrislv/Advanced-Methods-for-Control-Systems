% Parameters
N_particles = 1000; % Number of particles
N_steps = 15;
particles = zeros(N_particles, N_steps); % Particle matrix
weights = ones(N_particles, 1) / N_particles; % Initial weights
x_est_PF = zeros(1, N_steps); % State estimation

% Initial true state (given by data_generation.m)
x = randn(1);
w_std = sqrt(10);
v_std = sqrt(1);
alp = 1;
bet = 5;

% Generate measurements (same as data_generation.m)
y = zeros(1, N_steps);
for k = 1:N_steps
    w = randn(1) * w_std;
    x(k+1) = 0.5 * x(k) + bet * x(k) / (1 + x(k)^2) + 8 * cos(1.2 * k) + w;
    v = randn(1) * v_std;
    y(k) = alp * x(k) + x(k)^2 / 20 + v;
end
x = x(1:end-1);

% Initialization: Draw particles from prior N(0,1)
particles(:, 1) = randn(N_particles, 1);

% Particle Filter Loop without Resampling
for k = 2:N_steps
    % Prediction Step: Propagate particles using the system dynamics
    for i = 1:N_particles
        w = sqrt(10) * randn; % Process noise
        particles(i, k) = 0.5 * particles(i, k-1) + ...
                          bet * particles(i, k-1) / (1 + particles(i, k-1)^2) + ...
                          8 * cos(1.2 * k) + w;
    end
    
    % Update weights using measurement likelihood
    for i = 1:N_particles
        y_pred = alp * particles(i, k) + particles(i, k)^2 / 20; % Measurement model
        weights(i) = exp(-(y(k) - y_pred)^2 / (2 * v_std^2)); % Likelihood
    end
    
    % Normalize weights (no resampling)
    weights = weights / sum(weights);
    
    % State Estimation: Weighted mean of particles
    x_est_PF(k) = sum(weights .* particles(:, k));
end

% Plot Results
figure;
plot(1:N_steps, x, 'g', 'LineWidth', 1.5, 'DisplayName', 'True State');
hold on;
plot(1:N_steps, x_est_PF, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Particle Filter Without Resampling');
legend;
title('Particle Filter State Estimation Without Resampling');
xlabel('Time Step');
ylabel('State');
grid on;
