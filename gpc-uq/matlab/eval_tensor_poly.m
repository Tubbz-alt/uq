% Evaluate the polynomial of tensor product form.
%
% Input :
%     poly  : A function of polynomial returns a vector with the values of the
%             polynomials at x according to the index matrix ind_mat.
%     x     : Point to evaluate the polynomial. Can be a scalar or column vector.
%     dim   : dimension of x.
%     order : Order of the polynomial.
%     ind_mat: Index matrix.
%
% Output :
%     val   : value of the polynomial.
%
% Author : Xiu Yang
% Update : July 5th, 2012

function [val] = eval_tensor_poly(polynomial, x, dim, order, ind_mat)

[row, col] = size(ind_mat);
val = ones(row, 1);

val_mat = polynomial(x, order);

for i = 1:row
  indices = find(ind_mat(i,:));
  for j = 1:length(indices)
    val(i) = val(i)*val_mat(indices(j), ind_mat(i, indices(j))+1);
  end
end
