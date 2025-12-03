# -*- coding: utf-8 -*-
"""
Created on Sat May 10 20:41:04 2025

@author: Duange Guo
"""

import random
import numpy as np
import sympy as sp
import control as ct
import pandas as pd
import scipy.io as sio
from CircuitGenerator import vectfit_auto_rescale, model, SystemAnalyzer


class DatasetGenerator(object):
    """
    To generate the dataset for training
    """

    def __init__(self, source):
        """
        
        """
        
        self.source = source
        self.omega = np.logspace(-1, 3, 1000)  # frequency range: 0.1 ~ 1000 rad/s
        self.n_points = 3  # sampling point (0,0.5,1)
        
        self.Dataset = dict()
        self.analyzer = None
        
        # response from matlab dataset
        if isinstance(self.source, dict):
            pn = 8
            vn = 8
            qn = 8
            xn = 1
            
            data_path = self.source['path'][0]
            label_path = self.source['path'][1]
            path = [data_path, label_path]
            total_point_list = self.source['total_point']
            p_number = total_point_list[0]
            v_number = total_point_list[1]
            q_number = total_point_list[2]
            x_number = total_point_list[3]
            
            random_p = random.sample(range(1, p_number+1), pn)
            random_v = random.sample(range(1, v_number+1), vn)
            random_q = random.sample(range(1, q_number+1), qn)
            random_x = random.sample(range(1, x_number+1), xn)
            
            self.generate_matlab_dataset(path, random_p, random_v, random_q, random_x)
            
        # response from Minegrid model
        elif isinstance(self.source, list):
            test_f, test_s = self.generate_analytical_dataset()
        
        # response from dataset
        elif isinstance(self.source, str):
            test_f, test_s = self.generate_excel_dataset()
        
        else:
            print("source should be list containing analytical state space or path string")
            test_f, test_s = None, None
    

    def generate_matlab_dataset(self, path, random_p, random_v, random_q, random_x):
        data_path = path[0]
        label_path = path[1]
        iteration = 1
        for p in random_p:
            for v in random_v:
                for q in random_q:
                    for x in random_x:
                        file_name = f"P{p}_V{v}_Q{q}_X{x}.mat"
                        label_data = sio.loadmat(f'{label_path}\{file_name}')
                        operation_variable = ['p','v','q','x']
                        for var in operation_variable:
                            data = eval(f"{label_data[var]}")
                            if var not in self.Dataset:
                                self.Dataset[var] = []
                            self.Dataset[var].append(data[0][0]) 

                        mat_data = sio.loadmat(f'{data_path}\{file_name}')
                        A, B, C, D = mat_data['A'], mat_data['B'], mat_data['C'], mat_data['D']
                        sys_linearized = ct.ss(A, B, C, D)
                        # compute \mathcal V
                        mag, phase, w = ct.bode(sys_linearized, omega=self.omega, dB=True, Hz=False, plot=False)
                        mag = mag[0,0,:]
                        phase = phase[0,0,:]
                        test_f = mag * np.exp(1j * phase * np.pi / 180)
                        test_s = 1j * w
                        self._generate(test_f, test_s)
                        print(f"{file_name} is processed with total number:{iteration}")
                        iteration = iteration+1


    def generate_excel_dataset(self):
        # 使用pandas读取CSV文件到DataFrame
        Mdata = pd.read_csv(self.source)
        # 将DataFrame转换为numpy数组以便后续处理
        Mdata = Mdata.to_numpy()
        # 计算f(s)的样本：第0列为频率(omega)，第1列为幅度，第2列为角度（转换为弧度）
        test_f = Mdata[:, 1] * np.exp(1j * Mdata[:, 2] * np.pi / 180)
        w = Mdata[:, 0]
        # 将频率转换为复数域s=jω
        test_s = 1j * w
    
    
    def generate_analytical_dataset(self):
        """
        Generate the dataset for training
        """
        # p_samples = np.linspace(0, 1, self.n_points)
        # v_samples = np.linspace(0, 1, self.n_points)
        # q_samples = np.linspace(0, 1, self.n_points)
        # x_samples = np.linspace(0, 1, self.n_points)
        # # loop for dataset generation \mathcal O
        # for i, p_val in enumerate(p_samples):
        #     for j, v_val in enumerate(v_samples):
        #         for k, q_val in enumerate(q_samples):
        #             for l, x_val in enumerate(x_samples):
        #                 print(f"\nProcessing point (p={p_val:.1f}, v={v_val:.1f}, q={q_val:.1f}, x={x_val:.1f})")
        #                 # restore \mathcal O
        #                 operation_variable = ['p_val','v_val','q_val','x_val']
        #                 for var in operation_variable:
        #                     data = eval(var)
        #                     if var not in self.Dataset:
        #                         self.Dataset[var] = []
        #                     self.Dataset[var].append(data)        
        #                     A = self.source[0]
        #                     B = self.source[1]
        #                     C = self.source[2]
        #                     D = self.source[3]
        #                     x1, x2, x3, x4, x5, x6, u = sp.symbols('x1 x2 x3 x4 x5 x6 u')
        #                     A_lin = A.subs({x1: x1_val, x2: x2_val, u: u_val})
        #                     B_lin = B.subs({x1: x1_val, x2: x2_val, u: u_val})
        #                     A_np = np.array(A_lin, dtype=float)
        #                     B_np = np.array(B_lin, dtype=float)
        #                     C = np.array([[1, 0, 0, 1, 0, 0]])
        #                     D = np.array([[0]])
        #                     sys_linearized = ct.ss(A_np, B_np, C, D)
        #                     # compute \mathcal V
        #                     mag, phase, w = ct.bode(sys_linearized, omega=self.omega, dB=True, Hz=False, plot=False)
        #                     test_f = mag * np.exp(1j * phase * np.pi / 180)
        #                     test_s = 1j * w
                        

    def _generate(self, test_f, test_s):

        poles, residues, d, h = vectfit_auto_rescale(test_f, test_s)
        fitted = model(test_s, poles, residues, d, h)
        input_lines = [
            f'poles: {", ".join(f"{p.real:g}{p.imag:+g}j" for p in poles)}',
            f'residues: {", ".join(f"{r.real:g}{r.imag:+g}j" for r in residues)}',
            f'offset: {d:g}',
            f'slope: {h:g}'
        ]
        # compute \mathcal C
        self.analyzer = SystemAnalyzer()
        self.analyzer.analyze(input_lines=input_lines)
        branch_type = ['rc_params', 'rl_params', 'rlc_params']
        element_type = ['R','L','C','g_m']
        
        for item in branch_type:
            # the data type of rc is different from rl and rlc
            if isinstance(self.analyzer.output_data[item], list):
                branch_num = len(self.analyzer.output_data[item])
            else:
                branch_num = 1
            for element in element_type:
                for num in range(branch_num):  # index beginning with 0
                    # for RL and RLC circuit
                    if isinstance(self.analyzer.output_data[item], list):
                        i_d = self.analyzer.output_data[item][num]['id']
                        items = f"{item}_{i_d}_{element}"
                        data = self.analyzer.output_data[item][num][element]
                        ## special processing for gfm and gfl
                        rlc_num = 1
                        rl_item = ['e','f']
                        ##
                        if (num==rlc_num-1) and (item == 'rlc_params'):
                            for t_item in rl_item:
                                items_t = f"rl_params_{t_item}_{element}"
                                if items_t not in self.Dataset:
                                    self.Dataset[items_t] = []
                                self.Dataset[items_t].append(data)
                            continue
                        if items not in self.Dataset:
                            self.Dataset[items] = []
                        self.Dataset[items].append(data)
                    # for RC circuit
                    else:
                        items = f"{item}_{element}"
                        data = self.analyzer.output_data[item][element]
                        if items not in self.Dataset:
                            self.Dataset[items] = []
                        self.Dataset[items].append(data)
                            

if __name__ == '__main__':

    # 定义符号变量
    x1, x2, x3, x4, x5, x6, u = sp.symbols('x1 x2 x3 x4 x5 x6 u')

    # 定义非线性方程
    f1 = x2
    f2 = (-0.3*x1 + 0.2*(x3 - x1)) / 0.5
    f3 = x4
    f4 = (-0.2*(x3 - x1) + 0.1*(x5 - x3)) / 0.6
    f5 = x6
    f6 = (-0.1*(x5 - x3) +u) / 0.7

    # 计算雅可比矩阵（A, B）
    A = sp.Matrix([[sp.diff(f1, x1), sp.diff(f1, x2), sp.diff(f1, x3), sp.diff(f1, x4), sp.diff(f1, x5), sp.diff(f1, x6)],
                   [sp.diff(f2, x1), sp.diff(f2, x2), sp.diff(f2, x3), sp.diff(f2, x4), sp.diff(f2, x5), sp.diff(f2, x6)],
                   [sp.diff(f3, x1), sp.diff(f3, x2), sp.diff(f3, x3), sp.diff(f3, x4), sp.diff(f3, x5), sp.diff(f3, x6)],
                   [sp.diff(f4, x1), sp.diff(f4, x2), sp.diff(f4, x3), sp.diff(f4, x4), sp.diff(f4, x5), sp.diff(f4, x6)],
                   [sp.diff(f5, x1), sp.diff(f5, x2), sp.diff(f5, x3), sp.diff(f5, x4), sp.diff(f5, x5), sp.diff(f5, x6)],
                   [sp.diff(f6, x1), sp.diff(f6, x2), sp.diff(f6, x3), sp.diff(f6, x4), sp.diff(f6, x5), sp.diff(f6, x6)]])
    B = sp.Matrix([[sp.diff(f1, u)],
                   [sp.diff(f2, u)],
                   [sp.diff(f3, u)],
                   [sp.diff(f4, u)],
                   [sp.diff(f5, u)],
                   [sp.diff(f6, u)]])

    C = np.array([[1, 0, 0, 1, 0, 0]])
    D = np.array([[0]])

    source = [A,B,C,D]

    output_dir = r"E:\DataSynchrnous\BaiduSyncdisk\Python\tool\STNET\Test"
    
    filename = "test.csv"
    
    case1 = DatasetGenerator(source)
    
    # save_to_csv(case1.Dataset, output_dir, filename)
    
    # print("\n生成电路图中...")
    # drawer = CircuitDrawer(case1.analyzer)
    # drawer.draw_combined_circuits()
    # print("电路图已生成！")
    
    # system = ApparatusCircuit("test1", output_dir, filename)
    # system.train_circuit("mlp")




