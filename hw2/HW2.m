
%%1b
h = 0.7;  % Probability of rain transitioning to rain
u = 0.3;  
initial_prob = 0.7;  % Initial probability
time_steps = 30;  

time = zeros(1, time_steps + 1);
probability = zeros(1, time_steps + 1);

time(1) = 0;
probability(1) = initial_prob;

for t = 1:time_steps
    new_prob = h * probability(t) + u * (1 - probability(t));
    time(t + 1) = t;
    probability(t + 1) = new_prob;
end

figure;
plot(time, probability, 'bo-', 'LineWidth', 2);
title('Rain probability Over Time');
xlabel('Time Steps');
ylabel('Probability');
ylim([0, 1]);
%% 4a

sigma_x_values = linspace(0, 0.01, 4);
sigma_z_values = linspace(0, 0.01, 4);
t_values = 0:10;
sigma_t_values = zeros(length(t_values), length(sigma_x_values), length(sigma_z_values));

for i = 1:length(sigma_x_values)
    for j = 1:length(sigma_z_values)
        sigma_x = sigma_x_values(i);
        sigma_z = sigma_z_values(j);
        
        % Initialize σ²t at t=0
        sigma_t = 0; 
        
        for t = 1:length(t_values)
            sigma_t = (sigma_t + sigma_x) * sigma_z / (sigma_t + sigma_x + sigma_z);
            sigma_t_values(t, i, j) = sigma_t;
        end
    end
end

figure;
hold on;

for i = 1:length(sigma_x_values)
    for j = 1:length(sigma_z_values)
        plot(t_values, sigma_t_values(:, i, j), 'DisplayName', ['σ²x = ', num2str(sigma_x_values(i)), ', σ²z = ', num2str(sigma_z_values(j))]);
    end
end

xlabel('Time Steps (t)');
ylabel('σ²t');
title('σ²t behaviour');
legend('Location', 'best');
grid on;
hold off;
%% 6

%%
