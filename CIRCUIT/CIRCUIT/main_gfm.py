# -*- coding: utf-8 -*-
"""
Created on Wed Aug 20 15:27:18 2025

@author: Greg
"""

import numpy as np
import sympy as sp
from DatasetGeneration import DatasetGenerator
from CircuitTraining import ApparatusCircuit
from CircuitGenerator import CircuitDrawer
from ResultSaver import save_to_csv


# GFM
input_dir = r"E:\DataSynchrnous\BaiduSyncdisk\Matlab\论文\CircuitGenerator\VOC\new_training_data"
label_dir = r"E:\DataSynchrnous\BaiduSyncdisk\Matlab\论文\CircuitGenerator\VOC\training_label"
output_dir = r"E:\DataSynchrnous\BaiduSyncdisk\Python\tool\STNET\Test_GFM"

filename = "test_gfm.csv"
source = {"path": [input_dir, label_dir],
          "total_point":[5,5,5,20]
          }
case1 = DatasetGenerator(source)

save_to_csv(case1.Dataset, output_dir, filename)

# Network Preparation and Training
print("\n生成电路图中...")
drawer = CircuitDrawer(case1.analyzer)
drawer.draw_combined_circuits()
print("电路图已生成！")

system_gfm = ApparatusCircuit("test1", output_dir, filename)
system_gfm.train_circuit("pysr")