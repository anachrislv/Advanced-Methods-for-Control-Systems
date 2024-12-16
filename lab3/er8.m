% Particle Filter για το SEIR Μοντέλο με Πρόβλεψη 10 Βημάτων
% Αρχικοποίηση Παραμέτρων
N_particles = 1000; % Αριθμός σωματιδίων
N_steps = 200; % Χρονικά βήματα
N_future = 10; % Βήματα πρόβλεψης
dt = 0.1; % Μικρό βήμα χρόνου για σταθερότητα

% Παράμετροι του Μοντέλου
r1 = 0.1;
r2 = 0.3; 
b1 = 0.005;
b2 = 0.1;
b3 = 0.05;
b4 = 0.2;

% Τυχαίοι πολλαπλασιαστικοί θόρυβοι
sigma_w1 = 0.2; % Μικρότεροι θόρυβοι για σταθερότητα
sigma_w2 = 0.2;
sigma_w3 = 0.2;

% Αρχικές Συνθήκες
I0 = 0.01; % Αρχική τιμή I
E0 = 0.05; % Αρχική τιμή E
S = 1 - I0 - E0;
I = I0;
E = E0;

% Αρχικοποίηση Σωματιδίων
S_particles = ones(N_particles, 1) * S;
E_particles = ones(N_particles, 1) * E;
I_particles = ones(N_particles, 1) * I;

% Μετρήσεις y
y = zeros(1, N_steps);
x_true = zeros(3, N_steps);
x_true(:, 1) = [S; E; I];

% Αποθήκευση Εκτιμήσεων Κατάστασης
S_est = zeros(1, N_steps);
E_est = zeros(1, N_steps);
I_est = zeros(1, N_steps);

% Particle Filter Loop
for k = 1:N_steps
    % Γενιά Μετρήσεων
    y(k) = (0.2 * x_true(2, k) + x_true(3, k)) * exp(randn(1) * 0.1); % Μικρότερος θόρυβος στη μέτρηση
    
    % Πρόβλεψη για κάθε σωματίδιο
    for p = 1:N_particles
        mul_dist1 = exp(randn(1) * sigma_w1);
        mul_dist2 = exp(randn(1) * sigma_w2);
        mul_dist3 = exp(randn(1) * sigma_w3);

        S_particles(p) = S_particles(p) + dt * ((-r1 * E_particles(p) * S_particles(p) ...
                        - r2 * I_particles(p) * S_particles(p)) * mul_dist1 ...
                        + b1 * (1 - E_particles(p) - I_particles(p) - S_particles(p)) * mul_dist2);

        E_particles(p) = E_particles(p) + dt * ((r1 * S_particles(p) * E_particles(p) ...
                        + r2 * I_particles(p) * S_particles(p)) * mul_dist1 ...
                        - (b2 + b3) * E_particles(p) * mul_dist3);

        I_particles(p) = I_particles(p) + dt * (b2 * E_particles(p) - b4 * I_particles(p)) * mul_dist3;

        % Περιορισμός για αριθμητική σταθερότητα
        S_particles(p) = max(min(S_particles(p), 1), 0);
        E_particles(p) = max(min(E_particles(p), 1), 0);
        I_particles(p) = max(min(I_particles(p), 1), 0);
    end

    % Υπολογισμός Βαρών
    weights = exp(-0.5 * ((y(k) - (0.2 * E_particles + I_particles)).^2) / 0.1^2);
    weights = weights / sum(weights);

    % Systematic Resampling
    cumulative_weights = cumsum(weights);
    cumulative_weights(end) = 1; % Αποφυγή αριθμητικού λάθους
    resample_indices = zeros(N_particles, 1);
    u = (0:N_particles-1)' / N_particles + rand / N_particles;
    j = 1;

    for i = 1:N_particles
        while u(i) > cumulative_weights(j)
            j = j + 1;
        end
        resample_indices(i) = j;
    end

    S_particles = S_particles(resample_indices);
    E_particles = E_particles(resample_indices);
    I_particles = I_particles(resample_indices);

    % Εκτίμηση Κατάστασης
    S_est(k) = mean(S_particles);
    E_est(k) = mean(E_particles);
    I_est(k) = mean(I_particles);

    % Εξέλιξη πραγματικής κατάστασης
    x_true(1, k+1) = x_true(1, k) + dt * (-r1 * x_true(2, k) * x_true(1, k) ...
                      - r2 * x_true(3, k) * x_true(1, k) + b1 * (1 - sum(x_true(:, k))));
    x_true(2, k+1) = x_true(2, k) + dt * (r1 * x_true(1, k) * x_true(2, k) ...
                      + r2 * x_true(3, k) * x_true(1, k) - (b2 + b3) * x_true(2, k));
    x_true(3, k+1) = x_true(3, k) + dt * (b2 * x_true(2, k) - b4 * x_true(3, k));
end

% Πρόβλεψη 10 Βημάτων στο Μέλλον
S_future = S_particles;
E_future = E_particles;
I_future = I_particles;

S_pred = zeros(1, N_future);
E_pred = zeros(1, N_future);
I_pred = zeros(1, N_future);

for t = 1:N_future
    for p = 1:N_particles
        mul_dist1 = exp(randn(1) * sigma_w1);
        mul_dist2 = exp(randn(1) * sigma_w2);
        mul_dist3 = exp(randn(1) * sigma_w3);

        S_future(p) = S_future(p) + dt * ((-r1 * E_future(p) * S_future(p) ...
                     - r2 * I_future(p) * S_future(p)) * mul_dist1 ...
                     + b1 * (1 - E_future(p) - I_future(p) - S_future(p)) * mul_dist2);

        E_future(p) = E_future(p) + dt * ((r1 * S_future(p) * E_future(p) ...
                     + r2 * I_future(p) * S_future(p)) * mul_dist1 ...
                     - (b2 + b3) * E_future(p) * mul_dist3);

        I_future(p) = I_future(p) + dt * (b2 * E_future(p) - b4 * I_future(p)) * mul_dist3;
    end

    % Περιορισμός για σταθερότητα
    S_future = max(min(S_future, 1), 0);
    E_future = max(min(E_future, 1), 0);
    I_future = max(min(I_future, 1), 0);

    S_pred(t) = mean(S_future);
    E_pred(t) = mean(E_future);
    I_pred(t) = mean(I_future);
end

% Οπτικοποίηση Πρόβλεψης
figure;
plot(1:N_future, S_pred, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Predicted S');
hold on;
plot(1:N_future, E_pred, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Predicted E');
plot(1:N_future, I_pred, 'g-', 'LineWidth', 1.5, 'DisplayName', 'Predicted I');
legend;
title('Prediction of SEIR States 10 Steps Ahead');
xlabel('Future Time Steps');
ylabel('State Values');
grid on;
