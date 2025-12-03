# -*- coding: utf-8 -*-
"""
Created on Sat May  3 19:57:37 2025

@author: Duange Guo
"""

"""
Train the circuit model for the apparatus
"""

import numpy as np
from ParameterFunction import ParaFunction
from ResultSaver import load_training_dataset


class ApparatusCircuit(object):
    """
    The circuit class for apparatus
    """

    def __init__(self, name, para_path, filename):
        """
        name is the apparatus's name
        
        para_sheet is the csv data file containning O and C

        """
        
        self.name = name
        # dataset in panda form
        self.dataset = load_training_dataset(para_path, filename)
        # form the operation dataset
        col_names = ['p', 'v', 'q', 'x']
        self.O = self.dataset[col_names].to_numpy()
        
        rc_cols = [col for col in self.dataset.columns if col.startswith('rc_params_')]
        # 统计rl_params组数（按a/b/c...分组）
        rl_groups = set(col.split('_')[2] for col in self.dataset.columns 
                       if col.startswith('rl_params_') and len(col.split('_')) >= 4)
        # 统计rlc_params组数（按1/2/3...分组）
        rlc_groups = set(col.split('_')[2] for col in self.dataset.columns 
                        if col.startswith('rlc_params_') and len(col.split('_')) >= 4)

        self.op_num =  len(self.O)           # the number of operation points
        self.rc_num =  int(len(rc_cols)/4)   # the number of rc branch for a circuit
        self.rl_num =  len(rl_groups)        # the number of rl branch for a circuit
        self.rlc_num = len(rlc_groups)       # the number of rlc branch for a circuit
        
        # the function for each branch in one circuit
        self.branch_rc = {"R": None,
                          "C": None}
        self.branch_rl = [
            {
                "id": i,  # begin with 1
                "R": None,
                "L": None
            }
            for i in range(1, self.rl_num + 1)  # generate number self.rl_num
        ]
        self.branch_rlc = [
            {
                "id": i,  # begin with 1
                "R": None,
                "L": None,
                "C": None,
                "G": None
            }
            for i in range(1, self.rlc_num + 1)  # generate number self.rlc_num
        ]

    
    def train_circuit(self, network):
        """
        Train the circuit with chosen network: "pysr", "kan", "mlp"
        """
        
        # check firstly
        

        
        
        # for RC circuit
        Y_rc_r = self.dataset["rc_params_R"].values
        Y_rc_c = self.dataset["rc_params_C"].values

        RC_R = ParaFunction("RC_R", self.O, Y_rc_r, network)
        RC_C = ParaFunction("RC_C", self.O, Y_rc_c, network)
        self.branch_rc["R"] = RC_R
        self.branch_rc["C"] = RC_C
        
        # for RL circuit
        for i in range(self.rl_num):
            itemR = f"rl_params_{chr(97 + i)}_R"
            itemL = f"rl_params_{chr(97 + i)}_L"
            Y_rl_r = self.dataset[itemR].values
            Y_rl_l = self.dataset[itemL].values
            
            RL_R = ParaFunction(f"RL_R{i+1}", self.O, Y_rl_r, network)
            RL_L = ParaFunction(f"RL_L{i+1}", self.O, Y_rl_l, network)

            self.branch_rl[i]["R"]  = RL_R  # i is the list index
            self.branch_rl[i]["L"]  = RL_L
        
        # for RLC circuit
        for i in range(self.rlc_num):
            itemR = f"rlc_params_{i+1}_R"
            itemL = f"rlc_params_{i+1}_L"
            itemC = f"rlc_params_{i+1}_C"
            itemG = f"rlc_params_{i+1}_g_m"
            Y_rlc_r = self.dataset[itemR].values
            Y_rlc_l = self.dataset[itemL].values
            Y_rlc_c = self.dataset[itemC].values
            Y_rlc_g = self.dataset[itemG].values

            RLC_R = ParaFunction(f"RLC_R{i+1}", self.O, Y_rlc_r, network)
            RLC_L = ParaFunction(f"RLC_L{i+1}", self.O, Y_rlc_l, network)
            RLC_C = ParaFunction(f"RLC_C{i+1}", self.O, Y_rlc_c, network)
            RLC_G = ParaFunction(f"RLC_G{i+1}", self.O, Y_rlc_g, network)

            self.branch_rlc[i]["R"]  = RLC_R  # i is the list index
            self.branch_rlc[i]["L"]  = RLC_L
            self.branch_rlc[i]["C"]  = RLC_C
            self.branch_rlc[i]["G"]  = RLC_G
        
        
if __name__ == "__main__":
    pass

    
# system = ApparatusCircuit("test1", output_dir, filename)

# system.train_circuit("kan")

# system.train_circuit("pysr")

# system.train_circuit("mlp")
 