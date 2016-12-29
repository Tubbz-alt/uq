% Prepare following matrix E(f_i f_j)
% 
% Assume that f = \sum c_k \phi_k, we have f_{x_i} = \sum c_k \phi_k_{x_i},
% therefor, let c=(c_1, c_2, ...)^T, \phi = (\phi_1, \phi_2, ...)^T we have
%
% E(f_i f_j) = c^T E(\phi_{x_i} \phi_{x_j}^T) c. Here c is changing for each
% iteration while E(\phi_{x_i} \phi_{x_j}^T) can be computed offline. This code
% compute this matrix from exact integral.

% Author : Xiu Yang
% Date   : 12/31/2014
%
% Usage : assign dimension on line 19 and the order of gPC exansion on line 20.

clear;
clc;

dim = 12;
poly_order = 3;
num_basis = nchoosek(dim+poly_order, poly_order);
indx_mat = full_tensor(@tensor, dim, poly_order);

kernel = cell(dim, dim);

tstart = tic;

parfor i = 1:dim
  for j = 1:dim
    N = zeros(num_basis);
    indx_k = find(indx_mat(:,i));
    indx_l = find(indx_mat(:,j));
    indx_mat_k = indx_mat(indx_k,:);
    indx_mat_l = indx_mat(indx_l,:);
    c_i = sqrt(indx_mat_k(:,i));
    c_j = sqrt(indx_mat_l(:,j));
    indx_mat_k(:,i) = indx_mat_k(:,i) - 1;
    indx_mat_l(:,j) = indx_mat_l(:,j) - 1;

    for k = 1:size(indx_k,1)
      for l = 1:size(indx_l,1)
        if (indx_mat_k(k,:)==indx_mat_l(l,:))
          N(indx_k(k),indx_l(l)) = c_i(k)*c_j(l);
          continue;
        end
      end
    end
    kernel{i,j} = sparse(N); 
  end
end

tellapsed = toc(tstart)

save(['kernel_dim' num2str(dim) '_P' num2str(poly_order) '.mat'], 'kernel');
