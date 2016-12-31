""" Generate test data for gPC UQ. """
import numpy

def test_function(x):
    """ The test function evaluated over an array x.  If x is a matrix, this function sums over the
    columns. """
    y = x.sum(axis=1)
    return (y + 0.25*y*y + 0.025*y*y*y)


def load_data_from_file(input_path):
    """ Load input data from a file and return in an array """
    data = []
    with open(input_path, "rt") as input_file:
        for line in input_file:
            line.strip()
            data = data + [float(w) for w in line.split()]
    return numpy.array(data)

def reshape_data(data, num_dim, num_samp):
    """ Reshape data returns from load_data_from_file based on the given number of dimensions and
    samples. """
    return data.reshape((num_samp, num_dim))