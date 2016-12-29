% Use l1 minimization to compute the coefficients in a gpC surrogate model.
% Use subspace to compute the rotation.
%
%
% Author :  Xiu Yang
% Date   :  3/30/2015
% Update : 12/20/2016

% Check dim, num_sample and data set before use

clear;
close all;

load('data.mat');

% Dimension of the system, i.e., number of i.i.d. Gaussian random variables
dim = 12;
% Polynomial order of the Hermite expansion
poly_order = 3;
% Number of samples
num_sample = 120;
% Number of validation samples
num_valid_sample = 40;
valid_input = input(num_sample+1:num_sample+num_valid_sample,:);
valid_output = output(num_sample+1:num_sample+num_valid_sample,:);
input(num_sample+1:end,:)= [];
output(num_sample+1:end,:)= [];
% Total number of the basis functions
num_basis = nchoosek(dim+poly_order, poly_order);
% Load the pre-computed stiff matrix
load(['kernel_dim' num2str(dim) '_P' num2str(poly_order) '.mat']);
indx_mat = full_tensor(@tensor, dim, poly_order);
display('finish reading data');
% Number of iterations in the rotation procedure
num_iteration = 6;
% Coefficient in each iteration
recon_coef = zeros(num_basis, num_iteration+1);
% Rotation matrix
rotate_mat = eye(dim);

% Line 39 to 51  are preparing for the cross-validation. You can see "recon"
% means the part for reconstruction, and "valid" means the part for validation.
num_sets = 4;
% reconstruction index
recon_indx = cell(num_sets,1);
% validation index
valid_indx = cell(num_sets,1);
% reconstruction measure matrix
recon_matrix = cell(num_sets,1);
% validation measure matrix
valid_matrix = cell(num_sets,1);
% reconstruction observation
recon_obser = cell(num_sets,1);
% validation measure matrix
valid_obser = cell(num_sets,1);

% The measurement matrix
measure_mat = zeros(num_sample, num_basis);
for k = 1:num_sample
  measure_mat(k,:) = eval_tensor_poly(@hermite_norm, input(k,:)', dim, poly_order, indx_mat)';
end

for i = 1:num_sets
  valid_indx{i} = i:num_sets:num_sample;
  recon_indx{i} = setxor(1:num_sample, valid_indx{i});
end

% Cross-validation step
for i = 1:num_sets
  recon_matrix{i} = measure_mat(recon_indx{i}, :);
  valid_matrix{i} = measure_mat(valid_indx{i}, :);
  recon_obser{i} = output(recon_indx{i}, :);
  valid_obser{i} = output(valid_indx{i}, :);
end

delta0 = 0.05;
num_delta = 5;
d_delta = 10^(1/(num_delta-1));
delta_r = zeros(num_delta, 1);
delta_v = zeros(num_delta, 1);
opts = spgSetParms('verbosity', 0);
for i = 1:num_delta
  delta_r(i) = delta0;
  tmp = 0;
  for j = 1:num_sets
    c = spg_bpdn(recon_matrix{j}, recon_obser{j}, delta_r(i), opts);
    % check the error
    tmp = tmp + norm(valid_matrix{j}*c-valid_obser{j}, 2.0);
  end
  delta_v(i) = tmp/num_sets;
  delta0 = delta0/d_delta;
end
[tmp, min_ind] = min(delta_v);
delta = delta_r(min_ind)*sqrt(num_sets/(num_sets-1));

% This is the standard l1 minimization (or OMP, you can switch)
c = spg_bpdn(measure_mat, output, delta, opts);
recon_coef(:,1) = c;

% Examine the error by using the validation data
valid_measure_mat = zeros(num_valid_sample, num_basis);
for k = 1:num_valid_sample
  valid_measure_mat(k,:) = eval_tensor_poly(@hermite_norm, valid_input(k,:)', dim, poly_order, indx_mat)';
end
my_output = valid_measure_mat*recon_coef(:,1);
display('RMSE of standard l1');
norm(my_output-valid_output)/norm(valid_output)

% Rotate
display('Start rotation');
rotate_pnt = input;

tmp_delta = delta/(sqrt(num_sets/(num_sets-1)));

for iter = 1:num_iteration
  % Genrate gradient matrix
  A = zeros(dim);
  for i = 1:dim
    A(i,i) = c'*kernel{i,i}*c/2.0;
    for j = i+1:dim
      A(i,j) = c'*kernel{i,j}*c;
    end 
  end 
  A = A+A';

  % Compute the rotation matrix
  [U,S,V]=svd(A);

  % Rotate the sampling points
  rotate_pnt = rotate_pnt*U;
  rotate_mat = rotate_mat*U;

  % Generate new measurement matrix since we have "new" sampling ponits.
  for k = 1:num_sample
    measure_mat(k,:) = eval_tensor_poly(@hermite_norm, rotate_pnt(k,:)', dim, poly_order, indx_mat)';
  end

  % Cross-validation to estimate the error bound.
  for i = 1:num_sets
    recon_matrix{i} = measure_mat(recon_indx{i}, :);
    valid_matrix{i} = measure_mat(valid_indx{i}, :);
    recon_obser{i} = output(recon_indx{i}, :);
    valid_obser{i} = output(valid_indx{i}, :);
  end
  delta0 = tmp_delta;
  num_delta2 = 3;
  d_delta2 = 10^(1/(num_delta2-1));
  delta2_r = zeros(num_delta2, 1);
  delta2_v = zeros(num_delta2, 1);
  for i = 1:num_delta2
    delta2_r(i) = delta0;
    tmp = 0;
    for j = 1:num_sets
      c = spg_bpdn(recon_matrix{j}, recon_obser{j}, delta2_r(i), opts);
      % check the error
      tmp = tmp + norm(valid_matrix{j}*c-valid_obser{j}, 2.0);
    end
    delta2_v(i) = tmp/num_sets;
    delta0 = delta0/d_delta2;
  end
  [tmp, min_ind] = min(delta2_v);
  delta2 = delta2_r(min_ind)*sqrt(num_sets/(num_sets-1));
  
  c = spg_bpdn(measure_mat, output, delta2, opts);
  recon_coef(:,iter+1) = c;
end

% Examine the error by using the validation data
valid_input = valid_input*rotate_mat;
for k = 1:num_valid_sample
  valid_measure_mat(k,:) = eval_tensor_poly(@hermite_norm, valid_input(k,:)', dim, poly_order, indx_mat)';
end
my_output = valid_measure_mat*recon_coef(:,end);
display('RMSE after rotations');
norm(my_output-valid_output)/norm(valid_output)
