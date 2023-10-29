clear all; 
clc;
% express factor as a matrix
n1 = struct('parent', [], 'cpt', []);
n2 = struct('parent', [], 'cpt', []);
n3 = struct('parent', [1,2], 'cpt', []);
n4 = struct('parent', [3], 'cpt', []);
n5 = struct('parent', [3], 'cpt', []);


pb = 0.01;
pe = 0.009;
pa_be = 0.98;
pa_bNOTe = 0.89;
pa_NOTbe = 0.14;
pa_NOTbNOTe = 0.01;
pj_a = 0.65;
pj_NOTa =0.08;
pm_a = 0.94;
pm_NOTa = 0.03;

n1.cpt = [1, pb;
        0, 1-pb];
n2.cpt = [1, pe;
        0, 1-pe];
n3.cpt = [1, 1, 1, pa_be;
        1, 1, 0, pa_bNOTe;
        1, 0, 1, pa_NOTbe;
        1, 0, 0, pa_NOTbNOTe;
        0, 1, 1, 1-pa_be;
        0, 1, 0, 1-pa_bNOTe;
        0, 0, 1, 1-pa_NOTbe;
        0, 0, 0, 1-pa_NOTbNOTe];
n4.cpt = [1, 1, pj_a;
        1, 0, pj_NOTa;
        0, 1, 1-pj_a;
        0, 0, 1-pj_NOTa];
n5.cpt = [1, 1, pm_a;
        1, 0, pm_NOTa;
        0, 1, 1-pm_a;
        0, 0, 1-pm_NOTa];
bn = [n1, n2, n3, n4, n5];


f1 = to_factor(4, [4,5], [1,1], bn);
f2 = to_factor(3, [4,5], [1,1], bn);
ans1 = product_factor(f1, f2);
ans2 = sum_out(ans1, 3);

pb_jm = [ [1;0],[1;1], [1;1], eliminate(1, [4,5], [1,1], bn);
    [1;0], [1;1], [0;0],  eliminate(1, [4,5], [1,0], bn);
    [1;0], [0;0], [1;1], eliminate(1, [4,5], [0,1], bn);
    [1;0], [0;0], [0;0], eliminate(1,[4,5], [0, 0], bn)];

f_true_where = find(pb_jm(:, 1) == 1);
f_false_where = find(pb_jm(:, 1) == 0);
pb_jm = [ones(length(f_true_where), 1), pb_jm(f_true_where, 1:0), pb_jm(f_true_where, 1+1:end); ...
                        zeros(length(f_false_where), 1), pb_jm(f_false_where, 1:1-1), pb_jm(f_false_where, 1+1:end)];

pe_jm = [[1;1], [1;1], [1;0], eliminate(2, [4,5], [1,1], bn);
    [1;1], [0;0], [1;0], eliminate(2, [4,5], [1,0], bn);
    [0;0], [1;1], [1;0], eliminate(2, [4,5], [0,1], bn);
    [0;0], [0;0], [1;0], eliminate(2,[4,5], [0, 0], bn)];

f_true_where = find(pe_jm(:, 1) == 1);
f_false_where = find(pe_jm(:, 1) == 0);
pe_jm = [ones(length(f_true_where), 1), pe_jm(f_true_where, 1:0), pe_jm(f_true_where, 1+1:end); ...
                        zeros(length(f_false_where), 1), pe_jm(f_false_where, 1:1-1), pe_jm(f_false_where, 1+1:end)];


function dist = eliminate(x, e, obsv, bn)
    order = get_ordering(x, bn);
    factors = [];
    for i = length(order):-1:1
        var = order(i);
        factors = [factors, to_factor(var, e, obsv, bn)];
        if ~ismember(var, [e, x])
            tmp = []; % where relevant factors are stord
            tmp_ind = [];
            index = 1;
            
            for j = 1:length(factors)
                if find(factors(j).args == var)
                    tmp = [tmp, factors(j)]; % relevant factors
                    tmp_ind = [tmp_ind, j];
                end
            end
            
            % delete from factors array
            factors(tmp_ind) = [];
            
            product = struct('factor', [], 'args', []); % product of relevant factors
            for j = 1:length(tmp)
                product = product_factor(product, tmp(j));
            end
            
            factors = [factors, sum_out(product, var)]; % summed out
        end
    end

    % pointwise product at last
    for i = 1:length(factors)
        if i==1
            product = factors(i);
        else
            product = product_factor(product, factors(i));
        end
    end
    alpha = sum(product.factor, 1);
    alpha = 1/alpha(2);
    % Normalize
    for i = 1:length(product.factor)
        product.factor(i,end) = product.factor(i,end)*alpha;
    end
    dist = product.factor(:, end);
end

function o = get_ordering(var, bn)
    o = 1:length(bn);
    o(1) = var;
    i = 2;
    while i-1 ~= var
        o(i) = i-1;
        i = i+1;
    end
end

function f = to_factor(var, e, obsv, bn)
    n_parent = length(bn(var).parent);
    if ~n_parent
        f = struct('factor', bn(var).cpt, 'args', [var]);
    elseif ismember(var, e)
        is_true = obsv(find(e==var));
        cpt = bn(var).cpt(find(bn(var).cpt(:, 1) == is_true), 2:end);
        f = struct('factor', cpt, 'args', [bn(var).parent]);
    else
        f = struct('factor', bn(var).cpt, 'args', [var, bn(var).parent]);
    end
end

function f = product_factor(f1, f2)
    if isempty(f1.factor) && isempty(f1.args)
        f = f2;
        return;
    end
    % find common arg
    [common_arg, arg_index_f1, arg_index_f2] = intersect(f1.args, f2.args);
    % find which rows have common arg == true
    f1_true_where = find(f1.factor(:, arg_index_f1) == 1);
    f2_true_where = find(f2.factor(:, arg_index_f2) == 1);
    f1_false_where = find(f1.factor(:, arg_index_f1) == 0);
    f2_false_where = find(f2.factor(:, arg_index_f2) == 0);

    % make augmented f1 and f2 so taht args of common are in the far left
    f1_aug = [ones(length(f1_true_where), 1), f1.factor(f1_true_where, 1:arg_index_f1-1), f1.factor(f1_true_where, arg_index_f1+1:end); ...
                            zeros(length(f1_false_where), 1), f1.factor(f1_false_where, 1:arg_index_f1-1), f1.factor(f1_false_where, arg_index_f1+1:end)];
    f2_aug = [ones(length(f2_true_where), 1), f2.factor(f2_true_where, 1:arg_index_f2-1), f2.factor(f2_true_where, arg_index_f2+1:end); ...
                            zeros(length(f2_false_where), 1), f2.factor(f2_false_where, 1:arg_index_f2-1), f2.factor(f2_false_where, arg_index_f2+1:end)];

    fact = zeros(2^(length(f1.args)+length(f2.args)-1) , length(f1.args)+length(f2.args));
    
    for i = 1:2^length(f1.args)/2
        for j = 1:2^(length(f2.args)-1)
            fact((i-1)*2^(length(f2.args)-1) + j, 1:length(f1.args)) = f1_aug(i, 1:end-1); % table part of f1
            fact((i-1)*2^(length(f2.args)-1) + j, end) = f1_aug(i, end); % prob part of f1
        end
        fact((i-1)*2^(length(f2.args)-1)+1 : i*2^(length(f2.args)-1), length(f1.args) + 1 : end-1) = f2_aug(1:end/2, 2:end-1);
        fact((i-1)*2^(length(f2.args)-1)+1 : i*2^(length(f2.args)-1), end) = fact((i-1)*2^(length(f2.args)-1)+1 : i*2^(length(f2.args)-1), end) .* f2_aug(1:end/2, end);
    end
    for i = 2^length(f1.args)/2+1:2^length(f1.args)
        for j = 1:2^(length(f2.args)-1)
            fact((i-1)*2^(length(f2.args)-1) + j, 1:length(f1.args)) = f1_aug(i, 1:end-1);
            fact((i-1)*2^(length(f2.args)-1) + j, end) = f1_aug(i, end);
        end
        fact((i-1)*2^(length(f2.args)-1)+1 : i*2^(length(f2.args)-1), length(f1.args) + 1 : end-1) = f2_aug(end/2+1:end, 2:end-1);
        fact((i-1)*2^(length(f2.args)-1)+1 : i*2^(length(f2.args)-1), end) = fact((i-1)*2^(length(f2.args)-1)+1 : i*2^(length(f2.args)-1), end) .* f2_aug(end/2+1:end, end);
    end
    
    f = struct('factor', fact, 'args', [common_arg, f1.args(f1.args ~= common_arg), f2.args(f2.args ~= common_arg)]);
end

function fsum = sum_out(f, arg) % sum about arg
    arg_index = find(f.args == arg);
    if isempty(arg_index)
        fsum = f;
        return;
    end
    f_true_where = find(f.factor(:, arg_index) == 1);
    f_false_where = find(f.factor(:, arg_index) == 0);
    f_aug = struct('factor', [ones(length(f_true_where), 1), f.factor(f_true_where, 1:arg_index-1), f.factor(f_true_where, arg_index+1:end); ...
                            zeros(length(f_false_where), 1), f.factor(f_false_where, 1:arg_index-1), f.factor(f_false_where, arg_index+1:end)], 'args', [arg, f.args(f.args ~= arg)]);
    fm = zeros(size(f_aug.factor,1)/2, size(f_aug.factor,2) - 1);
    for i = 1:size(f_aug.factor, 1)/2
        fm(i, 1:end-1) = f_aug.factor(i, 2:end-1);
        if f_aug.factor(i, 1:end-1) ~= f_aug.factor(i + size(f_aug.factor, 1)/2, 1:end-1)
            disp("Error");
            return
        end
        fm(i, end) = f_aug.factor(i, end) + f_aug.factor(i + size(f_aug.factor, 1)/2, end);
    end
    fsum = struct('factor', fm, 'args', f_aug.args(2:end));
end