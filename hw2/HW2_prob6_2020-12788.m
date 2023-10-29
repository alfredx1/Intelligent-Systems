% Load data from 'KF_data_a.mat'
load('KF_data_a.mat');

xt = x;  % Ground-truth data
zt = y;  % Noisy measurements

%% 6a
x1_true = xt(1,:);
x2_true = xt(3,:);

figure;
plot(x1_true, x2_true, 'b-', 'LineWidth', 2, 'DisplayName', 'True Position');
xlabel('x1t');
ylabel('x2t');
title('Particle True Position');
legend('Location', 'Best');
grid on;
hold on;


%% 6b
figure;
plot(x1_true, x2_true, 'b-', 'LineWidth', 2, 'DisplayName', 'True Position');
hold on;
plot(zt(1,:), zt(2,:), 'r-o', 'MarkerSize', 2, 'DisplayName', 'Noisy Observations');
xlabel('x1t');
ylabel('x2t');
title('Particle Noisy Observations');
legend('Location', 'Best');
grid on;
hold on;


%% 6c
F = [1, 1, 0, 0; 0, 0.98, 0, 0; 0, 0, 1, 1; 0, 0, 0, 0.98];
H = [1, 0, 0, 0; 0, 0, 1, 0];
Sigma_w = eye(4);
Sigma_z = [2500, 0; 0, 2500];
x_hat = zeros(4, 1);
P = eye(4);

estimated_positions = [];

for t = 1:size(zt, 2)  % Change size(zt, 1) to size(zt, 2)
    % Generate process noise wt with N(0, 1) for dimensions 1 and 2
    wt = randn(4, 1);
    wt(3:4) = 0;  % Set dimensions 3 and 4 to zero
    
    % Generate measurement noise vt with N(0, Σz)
    vt = mvnrnd(zeros(2, 1), Sigma_z)';
    
    % Prediction step including process noise
    x_hat = F * x_hat + wt;
    P = F * P * F' + Sigma_w;
    
    % Update step
    z_hat = H * x_hat + vt;  % Include measurement noise vt
    y = zt(:, t) - z_hat;  
    S = H * P * H' + Sigma_z;
    K = P * H' / S;
    x_hat = x_hat + K * y;
    P = (eye(4) - K * H) * P;
    
    % Store estimated position
    estimated_positions = [estimated_positions; x_hat(1), x_hat(3)];
end
% Plot the true and estimated positions
% Plot the true and estimated positions
figure;
plot(x1_true, x2_true, 'b-', 'LineWidth', 2, 'DisplayName', 'True Position');
hold on;
plot(estimated_positions(:, 1), estimated_positions(:, 2), 'r-o', 'MarkerSize', 5, 'DisplayName', 'Estimated Position');
xlabel('x1t');
ylabel('x2t');
legend('Location', 'Best');
title('Kalman filter');
grid on;
hold off;

