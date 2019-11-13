#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 13 15:24:05 2019

@author: himanshu
"""

import numpy as np
from europeanOption import EuropeanOption

class MonteCarloPricing:
    
    def __init__(self,
                 option,
                 stepsize, 
                 numOfPaths, 
                 volatility, 
                 underlyingModel = 'GBM', 
                 discretizationMethod = 'Euler'):
        
        self.__option = option
        self.__stepsize = stepsize
        self.__paths = numOfPaths
        self.__sigma = volatility
        self.__model = underlyingModel
        self.__method = discretizationMethod
    
    def samplePaths(self, antithetic=False):
        # Set random seed
        if self.__method == 'Euler':
            # Generate Stock prices using Euler discretization
            St = self.__euler(antithetic)
        elif self.__method == 'Exact':
            #TODO: implement
            pass
        elif self.__method == 'Milstein':
            #TODO: implement
            pass
    
        return St
    
    #TODO: Implement antithetic paths
    def __euler(self, antithetic):
        S0 = self.__option.underlyingPrice
        r = self.__option.riskFreeRate
        T = self.__option.timeToExpiry
        
        np.random.seed(seed = 1231)
        num_steps = int(T/self.__stepsize)
        
        St = np.zeros((self.__paths, num_steps+1))
        St[:,0] = S0
        
        if self.__model == 'GBM':
            for i in range(1, num_steps+1):
                dZ = np.random.standard_normal(self.__paths)
                dW = np.sqrt(self.__stepsize)*dZ
                St[:,i] = St[:,i-1] + r*St[:,i-1]*self.__stepsize + self.__sigma*St[:,i-1]*dW
            return St
        elif self.__model == 'JD':
            #TODO: implement
            pass
        elif self.__model == 'Heston':
            #TODO: implement
            pass
        elif self.__model == 'OU':
            #TODO: implement
            pass
    
    #TODO: implement
    def __exact(self, St):
        pass
    
    #TODO: implement
    def __milstein(self, St):
        pass
    
    def optionPrice(self):
        K = self.__option.strikePrice
        r = self.__option.riskFreeRate
        T = self.__option.timeToExpiry
        
        St = self.samplePaths()
        if self.__option.optionType == 'Call':
            payoff = St[:,-1] - K
        elif self.__option.optionType == 'Put':
            payoff = K - St[:,-1]
            
        payoff[payoff < 0.0] = 0.0
        mc_val = np.exp(-r*(T))*np.mean(payoff)
        
        return "{0:.2f}".format(mc_val)
    
    def __repr__(self):
        
        return "MonteCarloPricing({}, {}, {}, {}, {}, {})"\
                .format(self.__option,
                        self.__stepsize,
                        self.__paths,
                        self.__sigma,
                        self.__model,
                        self.__method)
    
if __name__ == "__main__":
    S0 = 100
    K = 110
    r = 0.10
    T = 1
    volatility = 0.25
    
    print('------------------------------------------------------------------'
          +'----------------------------')
    option = EuropeanOption(S0, K, T, 'Call', r)
    mcPricing = MonteCarloPricing(option, 0.001, 1000, volatility)
    print('Monte-Carlo Estimate for Call:', mcPricing.optionPrice())
    
    print('------------------------------------------------------------------'
          +'----------------------------')
    
    option = EuropeanOption(S0, K, T, 'Put', r)
    mcPricing = MonteCarloPricing(option, 0.001, 1000, volatility)
    print('Monte-Carlo Estimate for Put:', mcPricing.optionPrice())
    print(mcPricing)
    