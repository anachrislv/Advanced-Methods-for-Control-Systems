X_min = -35;
X_max = 35;
N = 1001;

x_values = linspace(X_min, X_max, N); % Discretize the state space
P_Prior = zeros(N, 1); % Initialize the prior probability distribution
x_0_mean = 0; % Mean of the initial state
x_0_std = 1; % Standard deviation of the initial state

% Initialize Prior distribution
for i = 1:N
    P_Prior(i) = (1 / sqrt(2 * pi * x_0_std^2)) * exp(-(x_values(i) - x_0_mean)^2 / (2 * x_0_std^2));
end

x_est_Bayesian = zeros(1, 10); % Initialize the Bayesian estimation

for k = 1:10
    y_meas = y(k); % Current measurement at step k
    
    % Correction Step
    for x_val_ind = 1:N
        x_val_i = x_values(x_val_ind); % Current state value in the discretized space
        
        % Compute Likelihood (measurement model)
        likelihood_unnormalized = (1 / sqrt(2 * pi * v_std^2)) * ...
                                  exp(-(y_meas - (alp * x_val_i + x_val_i^2 / 20))^2 / (2 * v_std^2));
        % Unnormalized posterior probability
        p_posterior_unnormalized(x_val_ind) = likelihood_unnormalized * P_Prior(x_val_ind);
    end
    
    % Normalize the posterior probability
    normalization_factor = sum(p_posterior_unnormalized) * (x_values(2) - x_values(1));
    p_posterior = p_posterior_unnormalized / normalization_factor;

    % Υπολογισμός εκτίμησης ως μέση τιμή της posterior κατανομής
    x_est_Bayesian(k) = sum(x_values(:)' .* p_posterior(:)') * (x_values(2) - x_values(1));

    % Visualization of Prior and Posterior distributions
    figure
    plot(x_values, P_Prior, 'r', 'DisplayName', 'Prior');
    hold on;
    plot(x_values, p_posterior, 'b', 'DisplayName', 'Posterior');
    quiver(y_meas, 0, 0, 0.5, 'k', 'DisplayName', 'Measurement'); % Show the measurement
    quiver(x(k), 0, 0, 0.55, 'g', 'DisplayName', 'True State'); % Show the true state
    legend;
    title(['Time Step ', num2str(k)]);
    hold off;

    % Prediction Step
    for x_val_new_ind = 1:N
        x_val_i_new = x_values(x_val_new_ind); % Next state value
        for x_val_cur_ind = 1:N
            x_val_i_cur = x_values(x_val_cur_ind); % Current state value
            
            % Transition probability (state evolution model)
            transition_mean = 0.5 * x_val_i_cur + bet * x_val_i_cur / (1 + x_val_i_cur^2) + 8 * cos(1.2 * k);
            p_Transition_x_to_x_pl = (1 / sqrt(2 * pi * w_std^2)) * ...
                                     exp(-(x_val_i_new - transition_mean)^2 / (2 * w_std^2));
                                 
            % Contribution to the new prior
            p_x_new_to_integrate(x_val_cur_ind) = p_Transition_x_to_x_pl * p_posterior(x_val_cur_ind);
        end
        
        % Compute new Prior by summing contributions
        p_x_new(x_val_new_ind) = sum(p_x_new_to_integrate) * (x_values(2) - x_values(1));
    end
    
    % Update Prior for the next step
    P_Prior = p_x_new(:);
end
