""" Generate a gPC kernel matrix.

Prepare following matrix E(f_i f_j)

Assume that f = \sum c_k \phi_k, we have f_{x_i} = \sum c_k \phi_k_{x_i}, Therefore, let
c=(c_1, c_2, ...)^T, \phi = (\phi_1, \phi_2, ...)^T we have
E(f_i f_j) = c^T E(\phi_{x_i} \phi_{x_j}^T) c. Here c is changing for each iteration while
E(\phi_{x_i} \phi_{x_j}^T) can be computed offline. This code compute this matrix from exact
integral.

Author : Xiu Yang (translated by Nathan Baker)
Date   : 12/31/2014 (translated in 2016)

Usage : assign dimension and order of gPC expansion. """
from scipy.special import binom

def num_basis(dim, poly_order):
    """ Calculate the number of basis functions given the dimension and gPC order. """
    return binom(dim+poly_order, poly_order)

def kernel(dim, poly_order):
    """ Generate gPC kernel matrix with dimension dim and order poly_order. """
    nbasis = num_basis(dim, poly_order)
    # indx_mat = full_tensor(@tensor, dim, poly_order);
    # kernel = cell(dim, dim);
