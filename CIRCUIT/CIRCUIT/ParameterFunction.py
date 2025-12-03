# -*- coding: utf-8 -*-
"""
Created on Sat May  3 20:07:31 2025

@author: Duange Guo
"""

"""
The class of trained function for circuit parameter
"""


from sklearn.neural_network import MLPRegressor
from PINN.PINN_Net import PhysicsInformedNN
import torch
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split


class ParaFunction(object):
    """
    Nonlinear Function of Parameters
    """
    
    def __init__(self, name, trainx, trainy, network, config = []):
        """
        name is the type of parameters including: R, L, G, C
        network is pysr, kan or mlp
        
        trainx, and trainy is the input and output dataset pair
        """
        
        self.name = name
        self.network = network
        self.net_model = None
        
        self.trainx = trainx
        self.trainy = trainy
        
        self.config = config
        
        
        self.result = self.train_network()


    def _mlp_model(self):
        """return mlp"""
        
        model = MLPRegressor(hidden_layer_sizes=(10, 10), max_iter=1000, random_state=42)
        
        return model
    
    
    def _pinn_model(self, X_u_train, u_train):
        """return pinn"""
        
        layers = [2, 20, 20, 20, 20, 20, 20, 20, 20, 1]
        model = PhysicsInformedNN(X_u_train, u_train, layers, 0, 1)
        
        return model


    def train_network(self):
        """
        Prepare and train the network according to the user
        """
        
        if self.network == "mlp":
            # the dataset keeps the same with kan
            dataset = self._kan_dataset(self.trainx, self.trainy)
            self.net_model = self._mlp_model()
            self.net_model.fit(dataset["train_input"], dataset["train_label"])
            
            return None
        
        if self.network == "pinn":
            self.net_model = self._pinn_model(self.trainx, self.trainy)
            self.net_model.train(0)
            
    
    def calculate_error(self, train_data, test_data):
        """
        Calculate the error for comparison
        """
        
        pass
    
    
    def output_equation(self, goal):
        """
        Output the analytical equation
        """
        
        if goal == "matlab":
            pass
        
        if goal == "latex":
            pass
        
        if goal == "symbolic":
            pass
        
