function [vec_normalize] = hmm_normalize(given_vector)
% make sum of elements in given_vector be one.
% given_vector must have size of n X 1 ot 1 X n

n = length(given_vector);
summ = sum(given_vector);
vec_normalize = zeros(size(given_vector));
if summ==0
    fprintf('divide by zero');
end
vec_normalize = given_vector / summ;
end