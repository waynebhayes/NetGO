#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Copyright 2015, Institute for Systems Biology.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

Author: William Poole
Email: william.poole@systemsbiology.org / tknijnen@systemsbiology.org
Created: June 2015
"""

print("If this code fails (as Python is wont to do), try ebm.sh--which is also about 10-50x faster...")
print("... and if ebm.sh isn't fast enough, try the compiled C version (another 10x-100x faster) in my libwayne repo")

import sys
import numpy as np
from EmpiricalBrownsMethod import *
from scipy.stats import pearsonr

if len(sys.argv) < 2:
    raw_data = sys.stdin
else:
    raw_data = open(sys.argv[1])

data = []
pvals = []
lineNum = 0
for line in raw_data:
    if lineNum > 0:
        L = line.replace("\n", "").split("\t")
        pvals.append([float(L[1])])
        data.append([float(l) for l in L[2:]])
    lineNum = lineNum + 1
raw_data.close()

data = np.array(data)
pvals = np.array(pvals)
print("ebm", EmpiricalBrownsMethod(data, pvals))
