% Parameters
N_particles = 1000; % Number of particles
N_steps = 10; % Number of time steps
particles = zeros(N_particles, N_steps); % Particle matrix (with resampling)
particles_noRes = zeros(N_particles, N_steps); % Particle matrix (without resampling)
weights = ones(N_particles, 1) / N_particles; % Initial weights
weights_noRes = ones(N_particles, 1) / N_particles; % Initial weights for no resampling
x_est_PF = zeros(1, N_steps); % State estimation with resampling
x_est_PF_noRes = zeros(1, N_steps); % State estimation without resampling

% Initialize particles (Prior: N(0,1))
particles(:, 1) = randn(N_particles, 1);
particles_noRes(:, 1) = particles(:, 1);

% Particle Filter Loop WITH Resampling
for k = 2:N_steps
    % Prediction Step
    for i = 1:N_particles
        w = sqrt(10) * randn; % Process noise
        particles(i, k) = 0.5 * particles(i, k-1) + ...
                          5 * particles(i, k-1) / (1 + particles(i, k-1)^2) + ...
                          8 * cos(1.2 * k) + w;
    end
    
    % Update weights using measurement likelihood
    for i = 1:N_particles
        y_pred = 1 * particles(i, k) + particles(i, k)^2 / 20; % Measurement model
        weights(i) = exp(-(y(k) - y_pred)^2 / (2 * 1)); % Measurement noise variance = 1
    end
    
    % Normalize weights
    weights = weights / sum(weights);
    
    % Resampling Step
    cdf = cumsum(weights); % Cumulative distribution function
    resample_indices = zeros(N_particles, 1);
    for i = 1:N_particles
        u = rand;
        resample_indices(i) = find(cdf >= u, 1, 'first');
    end
    particles(:, k) = particles(resample_indices, k);
    weights = ones(N_particles, 1) / N_particles; % Reset weights after resampling
    
    % State Estimation: Mean of particles
    x_est_PF(k) = mean(particles(:, k));
end

% Particle Filter Loop WITHOUT Resampling
for k = 2:N_steps
    % Prediction Step
    for i = 1:N_particles
        w = sqrt(10) * randn; % Process noise
        particles_noRes(i, k) = 0.5 * particles_noRes(i, k-1) + ...
                                5 * particles_noRes(i, k-1) / (1 + particles_noRes(i, k-1)^2) + ...
                                8 * cos(1.2 * k) + w;
    end
    
    % Update weights using measurement likelihood
    for i = 1:N_particles
        y_pred = 1 * particles_noRes(i, k) + particles_noRes(i, k)^2 / 20; % Measurement model
        weights_noRes(i) = exp(-(y(k) - y_pred)^2 / (2 * 1)); % Measurement noise variance = 1
    end
    
    % Normalize weights
    weights_noRes = weights_noRes / sum(weights_noRes);
    
    % State Estimation: Weighted mean of particles
    x_est_PF_noRes(k) = sum(weights_noRes .* particles_noRes(:, k));
end

% Display State Estimations
disp('State estimation with Resampling:');
disp(x_est_PF);
disp('State estimation without Resampling:');
disp(x_est_PF_noRes);

% Plot Results
figure;
plot(1:N_steps, x(1:N_steps), 'g', 'LineWidth', 1.5, 'DisplayName', 'True State');
hold on;
plot(1:N_steps, x_est_PF, 'b--', 'LineWidth', 1.5, 'DisplayName', 'Particle Filter (with Resampling)');
plot(1:N_steps, x_est_PF_noRes, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Particle Filter (without Resampling)');
legend;
title('Particle Filter State Estimation');
xlabel('Time Step');
ylabel('State');
grid on;

% Compute MSE
MSE_PF = sum((x(1:N_steps) - x_est_PF).^2) / N_steps;
MSE_PF_noRes = sum((x(1:N_steps) - x_est_PF_noRes).^2) / N_steps;

fprintf('MSE for Particle Filter (with Resampling): %.4f\n', MSE_PF);
fprintf('MSE for Particle Filter (without Resampling): %.4f\n', MSE_PF_noRes);
