clc;
clear all;
close all;

%% HMM consists of 5 elements.
hmm.N = 3;      % number of states. (weather - Sunny(1), Cloudy(2), Rainy(3))
hmm.M = 4;      % number of observations. (activity - Walking(1), Shopping(2), Cleaning(3), Movies(4)
hmm.Pi = [0.4; 0.2; 0.4];        % initial probability.
hmm.A = [0.5 0.4 0.1; 0.3 0.3 0.4; 0.2 0.4 0.4]; % transition model.   N*N matrix.
                            % A_ij = probability of transition from i-th
                            % state to j-th state.
hmm.B = [0.3 0.2 0.35 0.15; 0.5 0.2 0.1 0.2; 0.15 0.25 0.15 0.45];     % observation model. N*M matrix.
                                        % B_ij = probability of getting j-th observation from i-th state. 
O = [1 1 4 2 3 3 2 4 2 1];  %observation seq   
% TODO: implement the functions 'hmm_filtering', 'hmm_smoothing',
% 'hmm_likelihood', and'hmm_viterbi'.
%% Problem (a)
[filtP, filtering_result] = hmm_filtering(hmm, O);
state_labels = {'Sunny', 'Cloudy', 'Rainy'};
time_step_labels = 1:size(filtering_result, 1);
fprintf('State (Weather): ');
disp(state_labels);
disp('Filtering Result (Posterior Probabilities):');
disp(filtering_result);

% %% Problem (b)
[smooP, smoothing_result] = hmm_smoothing(hmm, O, length(O));
% smoothing_result
activity_labels = {'Walking', 'Shopping', 'Cleaning', 'Movies'};

fprintf('Emission Probabilities (B matrix):\n');
for i = 1:hmm.N
    fprintf('State %d (Weather %s):\n', i, state_labels{i});
    fprintf('Observation  Activity\n');
    for j = 1:hmm.M
        fprintf('  %d:         %s: %.4f\n', j, activity_labels{j}, hmm.B(i, j));
    end
    fprintf('\n');
end
% %% Problem (c)
hmm_visualize_bar(filtering_result, smoothing_result);

% 
% %% Problem (d)
[seq_likelihood] = hmm_likelihood(hmm, O);
fprintf('Likelihood of the observation sequence: %.6f\n', seq_likelihood);

% %% Problem (e)
[optimal_seq] = hmm_viterbi(hmm, O);
most_likely_weather = state_labels(optimal_seq);

fprintf('Most likely sequence of weather based on observation seq:\n');
disp(most_likely_weather);

function [filtP, filtering_result] = hmm_filtering(hmm, O)
    N = hmm.N;
    T = length(O);

    filtering_result = zeros(N, T);
    filtP = hmm.Pi .* hmm.B(:, O(1)); % Initial filtering probabilities
    filtering_result(:, 1) = filtP;

    % Forward pass from 2 to T
    for k = 2:T
        predicted_filtP = hmm.A * filtP;
        filtP = predicted_filtP .* hmm.B(:, O(k));
        filtering_result(:, k) = filtP;
    end
end

function [smooP, smoothing_result] = hmm_smoothing(hmm, O, T)
    N = hmm.N;

    smoothing_result = zeros(N, T);
    beta = ones(N, 1);

    % Backward pass from T to 1
    for k = T:-1:1
        smoothing_result(:, k) = beta;
        smooP = hmm.Pi' * beta;
        if k > 1
            beta = hmm.A' * (beta .* hmm.B(:, O(k))) / sum(beta);
        end
    end
end


function state_label = getStateLabel(state)
    state_labels = {'Sunny', 'Cloudy', 'Rainy'};
    state_label = state_labels{state};
end


function [seq_likelihood] = hmm_likelihood(hmm, O)
    N = hmm.N;
    T = length(O);

    alpha = zeros(T, N);
    c = zeros(T, 1);
    seq_likelihood = 0;  % Initialize seq_likelihood as 0.

    % Initialization
    alpha(1, :) = hmm.Pi .* hmm.B(:, O(1));
    c(1) = 1 / sum(alpha(1, :));
    alpha(1, :) = alpha(1, :) * c(1);

    for t = 2:T
        alpha(t, :) = (alpha(t - 1, :) * hmm.A) .* hmm.B(:, O(t))';
        c(t) = 1 / sum(alpha(t, :));
        alpha(t, :) = alpha(t, :) * c(t);
    end

    % product of scaling factors c
    seq_likelihood = prod(c);
end



function optimal_seq = hmm_viterbi(hmm, observed_seq)
    N = hmm.N;
    T = length(observed_seq);

    delta = zeros(T, N);
    psi = zeros(T, N);
    
    % first delta values
    delta(1, :) = hmm.Pi .* hmm.B(:, observed_seq(1));
    
    % Forward pass
    for t = 2:T
        for j = 1:N
            % find prob and find the maximum
            temp = delta(t-1, :) .* hmm.A(:, j)';
            [delta(t, j), psi(t, j)] = max(temp);
            delta(t, j) = delta(t, j) * hmm.B(j, observed_seq(t));
        end
    end

    optimal_seq = zeros(1, T);
    [~, optimal_seq(T)] = max(delta(T, :));
    
    for t = T-1:-1:1
        optimal_seq(t) = psi(t+1, optimal_seq(t+1));
    end
end

