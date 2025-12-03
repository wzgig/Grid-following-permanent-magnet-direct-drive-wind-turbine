# -*- coding: utf-8 -*-
"""
Created on Sat May 10 20:41:13 2025

@author: Duange Guo
"""

"""
Result Saver
"""

import numpy as np
import pandas as pd
import os


def save_to_csv(data, output_dir, filename, system_name=None, overwrite=False):
    """
    保存数据到CSV（含所有扩展功能）
    
    参数:
        data: 包含omega/mag/phase的字典
        output_dir: 目标目录（自动创建）
        filename: 文件名（不含路径）
        system_name: 系统标注信息
        overwrite: 是否允许覆盖已有文件
    """

    # 处理路径（支持中文/空格）
    output_dir = os.path.abspath(output_dir)
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, filename)
    
    # 避免意外覆盖
    if os.path.exists(output_path) and not overwrite:
        raise FileExistsError(f"文件已存在：{output_path}（设置 overwrite=True 强制覆盖）")
    
    # 创建DataFrame
    df = pd.DataFrame(data)
    if system_name:
        df["System"] = system_name  # 添加系统标注
    
    # 保存（控制浮点数格式）
    df.to_csv(output_path, index=False, float_format="%.6f")
    print(f"数据已保存到：{output_path}")
    return output_path


def load_frequency_csv(output_dir, filename):
    """
    Load csv dataset of frequency response
    """
    
    csv_path = os.path.join(output_dir, filename)
    # 使用pandas读取CSV文件到DataFrame
    Mdata = pd.read_csv(csv_path)
    # 将DataFrame转换为numpy数组以便后续处理
    Mdata = Mdata.to_numpy()
    # 计算f(s)的样本：取前160行数据，第0列为幅度，第1列为角度（转换为弧度）
    test_f = Mdata[:, 1] * np.exp(1j * Mdata[:, 2] * np.pi / 180)
    w = Mdata[:, 0]
    # 将频率转换为复数域s=jω
    test_s = 1j * w
    
    return test_s, test_f


def load_training_dataset(output_dir, filename):
    """
    Load csv dataset of training_dataset pair
    """
    
    csv_path = os.path.join(output_dir, filename)
    dataset = pd.read_csv(csv_path)

    return dataset



