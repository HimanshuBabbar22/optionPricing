#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 13 15:24:05 2019

@author: himanshu
"""

import numpy as np

class MonteCarloPricing:
    
    def __init__(self,
                 option,
                 riskFreeRate,
                 volatility,
                 stepsize, 
                 numOfPaths, 
                 underlyingModel = 'GBM', 
                 discretizationMethod = 'Euler'):
        
        self.__option = option
        self.__r = riskFreeRate
        self.__sigma = volatility
        self.__stepsize = stepsize
        self.__paths = numOfPaths
        self.__model = underlyingModel
        self.__method = discretizationMethod
    
    def samplePaths(self):
        # Set random seed
        np.random.seed(seed = 1231)
        if self.__method == 'Euler':
            # Generate Stock prices using Euler discretization
            StPos, StNeg = self.__euler()
        elif self.__method == 'Exact':
            StPos, StNeg = self.__exact()
        return StPos, StNeg
    
    #TODO: Implement antithetic paths
    def __euler(self):
        S0 = self.__option.underlyingPrice
        r = self.__r
        T = self.__option.timeToExpiry
        sigma = self.__sigma
        
        numSteps = int(T/self.__stepsize)
        
        # Building antithetic paths
        StPos = np.zeros((self.__paths, numSteps+1))
        StNeg = np.zeros((self.__paths, numSteps+1))
        StPos[:,0] = S0
        StNeg[:,0] = S0
        
        if self.__model == 'GBM':
            for i in range(1, numSteps+1):
                dZ = np.random.standard_normal(self.__paths)
                dW = np.sqrt(self.__stepsize)*dZ
                StPos[:,i] = StPos[:,i-1] + r*StPos[:,i-1]*self.__stepsize + sigma*StPos[:,i-1]*dW
                StNeg[:,i] = StNeg[:,i-1] + r*StNeg[:,i-1]*self.__stepsize + sigma*StNeg[:,i-1]*(-dW)
        elif self.__model == 'JD':
            #TODO: implement
            pass
        elif self.__model == 'Heston':
            #TODO: implement
            pass
        elif self.__model == 'OU':
            #TODO: implement
            pass
        return StPos, StNeg
    
    def __exact(self):
        S0 = self.__option.underlyingPrice
        r = self.__r
        T = self.__option.timeToExpiry
        sigma = self.__sigma
        
        numSteps = int(T/self.__stepsize)
        # Building antithetic paths
        StPos = np.zeros((self.__paths, numSteps+1))
        StNeg = np.zeros((self.__paths, numSteps+1))
        if self.__model == 'GBM':
            dZ = np.arange(0, T+self.__stepsize, self.__stepsize)*np.random.standard_normal((self.__paths, numSteps+1))
            dt = np.arange(0, T+self.__stepsize, self.__stepsize)
            StPos = S0*np.exp((r - 0.5*(sigma**2))*dt + sigma*np.sqrt(dt)*dZ)
            StNeg = S0*np.exp((r - 0.5*(sigma**2))*dt + sigma*np.sqrt(dt)*(-dZ))
        return StPos, StNeg
    
    
    def optionPrice(self):
        K = self.__option.strikePrice
        r = self.__r
        T = self.__option.timeToExpiry
        
        StPos, StNeg = self.samplePaths()
        if self.__option.optionType == 'Call':
            payoffPos = StPos[:,-1] - K
            payoffNeg = StNeg[:,-1] - K
        elif self.__option.optionType == 'Put':
            payoffPos = K - StPos[:,-1]
            payoffNeg = K - StNeg[:,-1]
            
        payoffPos[payoffPos < 0.0] = 0.0
        payoffNeg[payoffNeg < 0.0] = 0.0
        mc_val = np.exp(-r*(T))*np.mean((payoffPos+payoffNeg)/2)
        
        return "{0:.2f}".format(mc_val)
    
    def __repr__(self):
        
        return "MonteCarloPricing({}, {}, {}, {}, {}, {}, {})"\
                .format(self.__option,
                        self.__r,
                        self.__sigma,
                        self.__stepsize,
                        self.__paths,
                        self.__model,
                        self.__method)
    
if __name__ == "__main__":
    
    from option import Option
    S0 = 100
    K = 110
    r = 0.10
    T = 1
    volatility = 0.25
    
    print('------------------------------------------------------------------'
          +'----------------------------')
    option = Option(S0, K, T, 'Call', 'European')
    mcPricing = MonteCarloPricing(option, r, volatility, 0.001, 1000, discretizationMethod='Euler')
    print(mcPricing)
    print('Monte-Carlo Estimate for Call using Euler discretization:', mcPricing.optionPrice())
    
    print('------------------------------------------------------------------'
          +'----------------------------')
    
    option = Option(S0, K, T,  'Put', 'European')
    mcPricing = MonteCarloPricing(option, r, volatility, 0.001, 1000, discretizationMethod='Euler')
    print(mcPricing)
    print('Monte-Carlo Estimate for Put using Euler discretization:', mcPricing.optionPrice())
    
    print('------------------------------------------------------------------'
          +'----------------------------')
    option = Option(S0, K, T, 'Call', 'European')
    mcPricing = MonteCarloPricing(option, r, volatility, 0.001, 1000, discretizationMethod='Exact')
    print(mcPricing)
    print('Monte-Carlo Estimate for Call using Exact discretization:', mcPricing.optionPrice())
    
    print('------------------------------------------------------------------'
          +'----------------------------')
    
    option = Option(S0, K, T,  'Put', 'European')
    mcPricing = MonteCarloPricing(option, r, volatility, 0.001, 1000, discretizationMethod='Exact')
    print(mcPricing)
    print('Monte-Carlo Estimate for Put using Exact discretization:', mcPricing.optionPrice())
    
    