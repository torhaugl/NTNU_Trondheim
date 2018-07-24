from input_file import *
from print_file import print_fozard, print_output
import classes
from numpy.random import random
from time import time
import datetime

now = datetime.datetime.now()

# Runs tests on Biofilm


# INITIALIZE BIOFILM
# 1-2 cells in 10 particles per vortex
bf = classes.Biofilm()
for vortex in bf.vortex_arr:
    if vortex.z == 0:
        for _ in range(10):
            vortex.add_cell(400 + 400 * random(), 2, 0)

# print before
print_fozard("data/before_small.dat", 0)
print("Start-time:", now.strftime("%Y-%m-%d %H:%M"))

# CALCULATE N STEPS
start_time = time()
classes.time_step(N, bf)
delta_time = time() - start_time

# print
print("Model time (min): %.2f" % (bf.time_step * dt))
print("Used time(min): %.2f" % (delta_time / 60))
print_fozard("data/small.out", delta_time)
print_output("data/small.dat", bf)
