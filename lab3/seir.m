% Particle Filter για το SEIR Μοντέλο
% Αρχικοποίηση Παραμέτρων
N_particles = 1000; % Αριθμός σωματιδίων
N_steps = 200; % Χρονικά βήματα
dt = 1; % Βήμα χρόνου

% Παράμετροι του Μοντέλου
r1 = 0.1;
r2 = 0.3; 
b1 = 0.005;
b2 = 0.1;
b3 = 0.05;
b4 = 0.2;

% Τυχαίοι πολλαπλασιαστικοί θόρυβοι
sigma_w1 = 0.5;
sigma_w2 = 0.3;
sigma_w3 = 0.3;

% Αρχικές Συνθήκες
I0 = 0.001 * rand(1);
S = 1 - I0;
E = 0;
I = I0;

% Αρχικοποίηση Σωματιδίων
S_particles = ones(N_particles, 1) * S;
E_particles = zeros(N_particles, 1);
I_particles = ones(N_particles, 1) * I;

% Μετρήσεις y
y = zeros(1, N_steps);
x_true = zeros(3, N_steps); % Πραγματική κατάσταση
x_true(:, 1) = [S; E; I];

% Αποθήκευση Εκτιμήσεων Κατάστασης
S_est = zeros(1, N_steps);
E_est = zeros(1, N_steps);
I_est = zeros(1, N_steps);

% Particle Filter Loop
for k = 1:N_steps
    % Γενιά Μετρήσεων
    y(k) = (0.2 * x_true(2, k) + x_true(3, k)) * exp(randn(1) * 0.3);
    
    % Προβλεψη για κάθε σωματίδιο
    for p = 1:N_particles
        % Τυχαίοι πολλαπλασιαστικοί θόρυβοι
        mul_dist1 = exp(randn(1) * sigma_w1);
        mul_dist2 = exp(randn(1) * sigma_w2);
        mul_dist3 = exp(randn(1) * sigma_w3);

        % Εξέλιξη καταστάσεων για κάθε σωματίδιο
        S_particles(p) = S_particles(p) + dt * ((-r1 * E_particles(p) * S_particles(p) ...
                        - r2 * I_particles(p) * S_particles(p)) * mul_dist1 ...
                        + b1 * (1 - E_particles(p) - I_particles(p) - S_particles(p)) * mul_dist2);

        E_particles(p) = E_particles(p) + dt * ((r1 * S_particles(p) * E_particles(p) ...
                        + r2 * I_particles(p) * S_particles(p)) * mul_dist1 ...
                        - (b2 + b3) * E_particles(p) * mul_dist3);

        I_particles(p) = I_particles(p) + dt * (b2 * E_particles(p) - b4 * I_particles(p)) * mul_dist3;
    end

    % Υπολογισμός Βαρών
    weights = exp(-0.5 * ((y(k) - (0.2 * E_particles + I_particles)).^2) / 0.3^2);
    weights = weights / sum(weights); % Κανονικοποίηση των βαρών

    % Systematic Resampling
    cumulative_weights = cumsum(weights);
    cumulative_weights(end) = 1; % Εξασφάλιση αριθμητικής ακρίβειας
    resample_indices = zeros(N_particles, 1);
    u = (0:N_particles-1)' / N_particles + rand / N_particles; % Ομοιόμορφη κατανομή
    j = 1;

    for i = 1:N_particles
        while u(i) > cumulative_weights(j)
            j = j + 1;
        end
        resample_indices(i) = j;
    end

    % Ανάδειγμα σωματιδίων
    S_particles = S_particles(resample_indices);
    E_particles = E_particles(resample_indices);
    I_particles = I_particles(resample_indices);

    % Εκτίμηση Κατάστασης
    S_est(k) = mean(S_particles);
    E_est(k) = mean(E_particles);
    I_est(k) = mean(I_particles);

    % Εξέλιξη πραγματικής κατάστασης
    mul_dist1_true = exp(randn(1) * sigma_w1);
    mul_dist2_true = exp(randn(1) * sigma_w2);
    mul_dist3_true = exp(randn(1) * sigma_w3);

    S_true = x_true(1, k) + dt * ((-r1 * x_true(2, k) * x_true(1, k) ...
             - r2 * x_true(3, k) * x_true(1, k)) * mul_dist1_true ...
             + b1 * (1 - x_true(2, k) - x_true(3, k) - x_true(1, k)) * mul_dist2_true);

    E_true = x_true(2, k) + dt * ((r1 * x_true(1, k) * x_true(2, k) ...
             + r2 * x_true(3, k) * x_true(1, k)) * mul_dist1_true ...
             - (b2 + b3) * x_true(2, k) * mul_dist3_true);

    I_true = x_true(3, k) + dt * (b2 * x_true(2, k) - b4 * x_true(3, k)) * mul_dist3_true;

    x_true(:, k+1) = [S_true; E_true; I_true];
end

% Plot Results
figure;
subplot(3,1,1);
plot(1:N_steps, x_true(1, 1:N_steps), 'g', 'LineWidth', 1.5, 'DisplayName', 'True S');
hold on;
plot(1:N_steps, S_est, 'b--', 'LineWidth', 1.5, 'DisplayName', 'Estimated S');
legend;
title('Susceptible State (S)');

subplot(3,1,2);
plot(1:N_steps, x_true(2, 1:N_steps), 'g', 'LineWidth', 1.5, 'DisplayName', 'True E');
hold on;
plot(1:N_steps, E_est, 'b--', 'LineWidth', 1.5, 'DisplayName', 'Estimated E');
legend;
title('Exposed State (E)');

subplot(3,1,3);
plot(1:N_steps, x_true(3, 1:N_steps), 'g', 'LineWidth', 1.5, 'DisplayName', 'True I');
hold on;
plot(1:N_steps, I_est, 'b--', 'LineWidth', 1.5, 'DisplayName', 'Estimated I');
legend;
title('Infectious State (I)');

sgtitle('Particle Filter State Estimation for SEIR Model');
