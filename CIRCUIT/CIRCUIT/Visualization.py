# -*- coding: utf-8 -*-
"""
Created on Sat May  3 20:08:20 2025

@author: Duange Guo
"""

"""
Visualization Tool for Sensitivity Analysis of STNET
"""

import plotly.graph_objects as go
from plotly.subplots import make_subplots
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np


def plot_training_data(X, y):
    fig = plt.figure(figsize=(16, 16))
    
    # 3D图1: Feature1, Feature4, Output
    ax1 = fig.add_subplot(111, projection='3d')
    scatter1 = ax1.scatter(X[:, 0], X[:, 1], X[:, 2], c=y, cmap='viridis', s=100)
    ax1.set_xlabel('Feature 1 (constant)')
    ax1.set_ylabel('Feature 4')
    ax1.set_zlabel('Output')
    ax1.set_title('3D: Feature1, Feature4, Output')
    


def plot_single_feature(X, y):
    # 创建图形
    plt.figure(figsize=(12, 8))
    
    # 第四特征 vs 输出（主要关系）
    fea = 4
    plt.subplot(2, 2, 1)
    plt.scatter(X[:, fea-1], y, alpha=0.7, color='red', s=80)
    plt.xlabel(f'Feature {fea}')
    plt.ylabel('Output')
    plt.title(f'Main Relationship: Feature {fea} vs Output')
    plt.grid(True, alpha=0.3)
    
    # 添加趋势线
    z = np.polyfit(X[:, fea-1], y, 1)
    p = np.poly1d(z)
    plt.plot(X[:, fea-1], p(X[:, fea-1]), "r--", alpha=0.8, 
             label=f'Trend: y = {z[0]:.2f}x + {z[1]:.2f}')
    plt.legend()
    
    # 输出值的分布（查看异常值）
    plt.subplot(2, 2, 2)
    plt.hist(y, bins=10, alpha=0.7, color='green', edgecolor='black')
    plt.xlabel('Output Value')
    plt.ylabel('Frequency')
    plt.title('Distribution of Output Values')
    plt.axvline(np.mean(y), color='red', linestyle='--', 
                label=f'Mean: {np.mean(y):.2f}')
    plt.legend()
    
    # 第四特征的分布
    plt.subplot(2, 2, 3)
    plt.hist(X[:, fea-1], bins=10, alpha=0.7, color='blue', edgecolor='black')
    plt.xlabel(f'Feature {fea} Value')
    plt.ylabel('Frequency')
    plt.title(f'Distribution of Feature {fea}')
    
    # 箱线图查看输出值的统计
    plt.subplot(2, 2, 4)
    plt.boxplot(y, vert=True)
    plt.ylabel('Output Value')
    plt.title('Boxplot of Output Values')
    
    plt.tight_layout()
    plt.show()

