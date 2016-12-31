""" Generate test data for gPC UQ """
import numpy
from numpy.linalg import norm
import unittest
import os
from .gen_data import test_function, load_data_from_file, reshape_data

TEST_PATH = os.path.dirname(__file__)
TEST_INPUT_PATH = os.path.join(TEST_PATH, "matlab-output/input.txt")
TEST_OUTPUT_PATH = os.path.join(TEST_PATH, "matlab-output/output.txt")
TEST_INPUT_DIM = 12
TEST_INPUT_SAMPLES = 160
TEST_EPSILON = 0.0001


class TestFunction(unittest.TestCase):
    """ Test the input function """
    def setUp(self):
        """ Set up test case """
        array = load_data_from_file(TEST_INPUT_PATH)
        self.input = array.reshape((TEST_INPUT_SAMPLES, TEST_INPUT_DIM))
        self.output = load_data_from_file(TEST_OUTPUT_PATH)
    def test_sum(self):
        """ Test the test function sum """
        fun = test_function(self.input)
        test_sum = fun.sum()
        real_sum = 496.5897
        self.assertTrue(numpy.abs(test_sum - real_sum) < TEST_EPSILON)
    def test_norm(self):
        """ Test the test function norm """
        fun = test_function(self.input)
        test_norm = norm(fun, 2)
        real_norm = 103.2338
        self.assertTrue(numpy.abs(test_norm - real_norm) < TEST_EPSILON)

class TestLoad(unittest.TestCase):
    """ Test data loading """
    def setUp(self):
        """ Set up test case """
        self.input = load_data_from_file(TEST_INPUT_PATH)
        self.output = load_data_from_file(TEST_OUTPUT_PATH)
    def test_input_sum(self):
        """ Test data loading correctly """
        test_sum = self.input.sum()
        real_sum = 60.8333
        self.assertTrue(numpy.abs(test_sum - real_sum) < TEST_EPSILON)
    def test_output_sum(self):
        """ Test data loading correctly """
        test_sum = self.output.sum()
        real_sum = 496.5897
        self.assertTrue(numpy.abs(test_sum - real_sum) < TEST_EPSILON)
    def test_shape(self):
        """ Test data shaping """
        shaped = reshape_data(data=self.input, num_dim=TEST_INPUT_DIM, num_samp=TEST_INPUT_SAMPLES)
        diag = numpy.diag(shaped)
        test_sum = diag.sum()
        real_sum = 4.1207
        self.assertTrue(numpy.abs(test_sum - real_sum) < TEST_EPSILON)

if __name__ == "__main__":
    unittest.main()