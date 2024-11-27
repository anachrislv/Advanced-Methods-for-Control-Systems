% Έλεγχος Δυνατότητας Ελέγχου
Co = ctrb(A, B);
rank_Co = rank(Co);

disp('Βαθμίδα Πίνακα Ελέγχου:');
disp(rank_Co);

if rank_Co < size(A, 1)
    disp('Το σύστημα δεν είναι πλήρως ελέγξιμο. Χρήση ελέγξιμου υποσυστήματος...');

    % Εργασία με ελέγξιμο υποσύστημα
    A_ctrl = A(1:rank_Co, 1:rank_Co);
    B_ctrl = B(1:rank_Co, :);
    Q_ctrl = eye(rank_Co); % Μοναδιαία μήτρα για καταστάσεις
    R_ctrl = eye(size(B_ctrl, 2)); % Μοναδιαία μήτρα για εισόδους

    try
        [K_ctrl, S_ctrl, e_ctrl] = lqr(A_ctrl, B_ctrl, Q_ctrl, R_ctrl);
        disp('Πίνακας Κέρδους για το Ελέγξιμο Υποσύστημα:');
        disp(K_ctrl);
        disp('Ιδιοτιμές Κλειστού Βρόχου:');
        disp(e_ctrl);

        % Μετατροπή του Κέρδους Πίσω στο Αρχικό Σύστημα
        K = zeros(size(B, 2), size(A, 1));
        K(:, 1:rank_Co) = K_ctrl;
        disp('Πίνακας Κέρδους για το Αρχικό Σύστημα:');
        disp(K);
    catch ME
        disp('Αδυναμία υπολογισμού του LQR ελεγκτή.');
        disp(ME.message);
    end
else
    disp('Το σύστημα είναι πλήρως ελέγξιμο. Υπολογισμός LQR...');

    % Απλός Υπολογισμός LQR
    Q = eye(size(A, 1)); % Μοναδιαία μήτρα για καταστάσεις
    R = eye(size(B, 2)); % Μοναδιαία μήτρα για εισόδους
    try
        [K, S, e] = lqr(A, B, Q, R);
        disp('Πίνακας Κέρδους K:');
        disp(K);
        disp('Ιδιοτιμές Κλειστού Βρόχου:');
        disp(e);
    catch ME
        disp('Αδυναμία υπολογισμού του LQR ελεγκτή.');
        disp(ME.message);
    end
end
