import os
import sys
import unittest

# This is added to load the modules from parent directory
sys.path.insert(1, os.path.dirname(os.getcwd()))
sys.path.insert(1, os.path.dirname(os.getcwd().rsplit("/", 1)[0]))


if __name__ == '__main__':
    suite = unittest.TestLoader().discover(start_dir=os.getcwd() + '/syntheticDataGeneration', pattern="test_*")
    unittest.TextTestRunner(verbosity=2).run(suite)

