function [] = hmm_visualize_bar(filtering_result, smoothing_result)
% Visualize the sequential variations of filtering and smoothing results
% using 2D bar graph.
% filtering_result : N*T matrix.
% smoothing_result : N*T matrix.
% No output.

if size(filtering_result) ~= size(smoothing_result)
    error('dimension mismatch!');
end
N = size(filtering_result,1);
T = size(filtering_result,2);

x_filt = zeros(1, T);
y_filt = zeros(T, N);
x_smoo = zeros(1, T);
y_smoo = zeros(T, N);


for i=1:T
    x_filt(i) = i;
    x_smoo(i) = i;
    for j=1:N
        y_filt(i, j) = filtering_result(j,i);
        y_smoo(i, j) = smoothing_result(j,i);
    end
    % Normalize for who didn't normalize the probabilities by mistake.
    norm_y_filt(i,:) = hmm_normalize(y_filt(i,:));
    norm_y_smoo(i,:) = hmm_normalize(y_smoo(i,:));
end


subplot(2,1,1);
bar(x_filt, norm_y_filt);
xlabel('Sequence');
ylabel('Probability(filtering)');
axis([0.5 T+0.5  0 1]);
legend('Sunny', 'Cloudy', 'Rainy');

subplot(2,1,2);
bar(x_smoo, norm_y_smoo);
xlabel('Sequence');
ylabel('Probability(smoothing)');
axis([0.5 T+0.5 0 1]);
legend('Sunny', 'Cloudy', 'Rainy');

        
