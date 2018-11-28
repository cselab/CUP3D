"""
Python wrapper for the pybind11-based .so.
"""

import os
import sys
sys.path.append(os.path.dirname(__file__))
from _cubismup3d import *  # For now just import everything.
sys.path.pop()


# Add any Python code here, such as various utility functions.
