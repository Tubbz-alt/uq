""" Test the generation of the kernel matrix """
import unittest
from gpc import kernel

TEST_INPUT_DIM = 12
TEST_POLY_ORDER = 3

class TestFunction(unittest.TestCase):
    """ Test the kernel matrix """
    def setUp(self):
        """ Setup test cases """
        pass
    def test_num_basis(self):
        """ Test the number of basis functions """
        test_num_basis = kernel.num_basis(dim=TEST_INPUT_DIM, poly_order=TEST_POLY_ORDER)
        real_num_basis = 455
        self.assertEqual(real_num_basis, test_num_basis)

if __name__ == "__main__":
    unittest.main()