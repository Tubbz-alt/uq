% Generate input and output samples by using a function as the black box. In
% practice, the black box is your simulation code, experiment, etc.
%
% Author : Xiu Yang
% Date   : 12/20/2016

clear;
close all;

% Dimension of the system, i.e., number of i.i.d. Gaussian random variables
dim = 12;
% Number of samples
num_sample = 160;
% "Black box", here we use a third-order polynomial
f = @(x) sum(x) + 0.25*sum(x)^2 + 0.025*sum(x)^3;
% Generate input samples
randn('seed', 1);
% The input is a MxP matrix, where M is the number of samples and P is the
% dimension.
input = randn(dim, num_sample)';
% Compute the output samples from the black box
output = zeros(num_sample, 1);
for k = 1:num_sample
  output(k) = f(input(k,:));
end

save('data.mat', 'input', 'output');

% Write out inputs and outputs for testing  purposes
save('test-output/input.txt', 'input', '-ascii', '-tabs')
save('test-output/output.txt', 'output', '-ascii', '-tabs')
