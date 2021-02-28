import numpy as np
import matplotlib.pyplot as plt

import ZLFFI12

def AddData0(data0):
  data0[62,0] = 0.252525
  data0[62,1] = 0.01
  data0[62,2] = 0
  data0[62,3] = 1.25253
  data0[62,4] = 0.286388
  data0[63,0] = 0.272727
  data0[63,1] = 0.01
  data0[63,2] = 0
  data0[63,3] = 1.27273
  data0[63,4] = 0.28408
  data0[64,0] = 0.292929
  data0[64,1] = 0.01
  data0[64,2] = 0
  data0[64,3] = 1.29293
  data0[64,4] = 0.281771
  data0[65,0] = 0.313131
  data0[65,1] = 0.01
  data0[65,2] = 0
  data0[65,3] = 1.31313
  data0[65,4] = 0.278477
  data0[72,0] = 0.454545
  data0[72,1] = 0.01
  data0[72,2] = 0
  data0[72,3] = 1.45455
  data0[72,4] = 0.249494
  data0[73,0] = 0.474747
  data0[73,1] = 0.01
  data0[73,2] = 0
  data0[73,3] = 1.47475
  data0[73,4] = 0.244833
  data0[74,0] = 0.494949
  data0[74,1] = 0.01
  data0[74,2] = 0
  data0[74,3] = 1.49495
  data0[74,4] = 0.239579
  data0[78,0] = 0.575758
  data0[78,1] = 0.01
  data0[78,2] = 0
  data0[78,3] = 1.57576
  data0[78,4] = 0.215697
  data0[79,0] = 0.59596
  data0[79,1] = 0.01
  data0[79,2] = 0
  data0[79,3] = 1.59596
  data0[79,4] = 0.208667
  data0[80,0] = 0.616162
  data0[80,1] = 0.01
  data0[80,2] = 0
  data0[80,3] = 1.61616
  data0[80,4] = 0.201637
  data0[87,0] = 0.757576
  data0[87,1] = 0.01
  data0[87,2] = 0
  data0[87,3] = 1.75758
  data0[87,4] = 0.145038
  data0[88,0] = 0.777778
  data0[88,1] = 0.01
  data0[88,2] = 0
  data0[88,3] = 1.77778
  data0[88,4] = 0.135542
  data0[93,0] = 0.878788
  data0[93,1] = 0.01
  data0[93,2] = 0
  data0[93,3] = 1.87879
  data0[93,4] = 0.0838003
  data0[94,0] = 0.89899
  data0[94,1] = 0.01
  data0[94,2] = 0
  data0[94,3] = 1.89899
  data0[94,4] = 0.072452
  data0[95,0] = 0.919192
  data0[95,1] = 0.01
  data0[95,2] = 0
  data0[95,3] = 1.91919
  data0[95,4] = 0.0611037
  data0[96,0] = 0.939394
  data0[96,1] = 0.01
  data0[96,2] = 0
  data0[96,3] = 1.93939
  data0[96,4] = 0.0497554
  data0[97,0] = 0.959596
  data0[97,1] = 0.01
  data0[97,2] = 0
  data0[97,3] = 1.9596
  data0[97,4] = 0.038407
  done=True


  ZLFFI12.AddData0(data0)
