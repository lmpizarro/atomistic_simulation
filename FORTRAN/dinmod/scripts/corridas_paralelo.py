# -*- coding: utf-8 -*-

#!/usr/bin/python

###############################################################################
#
# Script para correr el programa de dinámica molecular en paralelo

import os
import os.path
import shutil
import subprocess
import errno
import numpy as np

# Paralelización
from mpi4py import MPI

