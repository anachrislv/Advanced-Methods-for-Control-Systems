% Ορισμός Παραμέτρων Συστήματος
A = [0, 1, 0, 0, 0, 0;
     0, 0, -9.81, 0, 0, 0;
     0, 0, 0, 1, 0, 0;
     0, 0, 9.81, 0, 0, 0;
     0, 0, 0, 0, 0, 1;
     0, 0, 0, 0, 0, 0]; % Παράδειγμα δυναμικού συστήματος
B = [0, 0, 0;
     1, 0, 0;
     0, 0, 0;
     0, 1, 0;
     0, 0, 0;
     0, 0, 1]; % Είσοδοι

C = eye(6); % Έξοδοι
D = zeros(size(C, 1), size(B, 2)); % Χωρίς απευθείας όρους

% Έλεγχος Ελέγξιμου Συστήματος
Co = ctrb(A, B);
rank_Co = rank(Co);

if rank_Co < size(A, 1)
    disp('Το σύστημα δεν είναι πλήρως ελέγξιμο.');
    return;
else
    disp('Το σύστημα είναι πλήρως ελέγξιμο.');
end

% Ορισμός Τροχιάς
trajectory = [0, 0, 0;
              1, 0, 0;
              1, 1, 0;
              1, 1, 1;
              0, 1, 1;
              0, 0, 1;
              0, 0, 0]; % Σημεία της τροχιάς

time_step = 0.1; % Βήμα χρόνου
T = 10; % Συνολικός χρόνος
time = 0:time_step:T;

% Προετοιμασία Αποτελεσμάτων
x = zeros(6, length(time));
u = zeros(3, length(time));

% Ρύθμιση LQR Ελεγκτή
Q = diag([100, 10, 100, 10, 100, 10]); % Βαρύτητες καταστάσεων
R = diag([0.1, 0.1, 0.1]); % Βαρύτητες εισόδων

[K, ~, ~] = lqr(A, B, Q, R);

% Προσομοίωση
current_state = zeros(6, 1); % Αρχική κατάσταση
trajectory_index = 1;

for k = 1:length(time)
    % Ενημέρωση επιθυμητής θέσης
    if trajectory_index <= size(trajectory, 1)
        desired_position = trajectory(trajectory_index, :)';
    else
        desired_position = trajectory(end, :)';
    end

    % Υπολογισμός Σφάλματος
    error = desired_position - current_state(1:3);

    % Υπολογισμός Εισόδου
    u(:, k) = -K * [error; current_state(4:end)];

    % Ενημέρωση κατάστασης
    current_state = A * current_state + B * u(:, k);

    % Μετάβαση στο επόμενο σημείο της τροχιάς
    if norm(error) < 0.1 && trajectory_index < size(trajectory, 1)
        trajectory_index = trajectory_index + 1;
    end

    % Αποθήκευση κατάστασης
    x(:, k) = current_state;
end

% Γραφική Απεικόνιση
figure;
plot3(trajectory(:, 1), trajectory(:, 2), trajectory(:, 3), 'ro-', 'LineWidth', 2);
hold on;
plot3(x(1, :), x(2, :), x(3, :), 'b-', 'LineWidth', 2);
grid on;
xlabel('x');
ylabel('y');
zlabel('z');
legend('Desired Trajectory', 'Actual Trajectory');
title('Tracking Performance');
