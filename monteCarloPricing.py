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
    
    def samplePaths(self, antithetic=False):
        # Set random seed
        np.random.seed(seed = 1231)
        if self.__method == 'Euler':
            # Generate Stock prices using Euler discretization
            St = self.__euler(antithetic)
        elif self.__method == 'Exact':
            St = self.__exact(antithetic)
        return St
    
    #TODO: Implement antithetic paths
    def __euler(self, antithetic):
        S0 = self.__option.underlyingPrice
        r = self.__r
        T = self.__option.timeToExpiry
        sigma = self.__sigma
        
        numSteps = int(T/self.__stepsize)
        
        St = np.zeros((self.__paths, numSteps+1))
        St[:,0] = S0
        
        if self.__model == 'GBM':
            for i in range(1, numSteps+1):
                dZ = np.random.standard_normal(self.__paths)
                dW = np.sqrt(self.__stepsize)*dZ
                St[:,i] = St[:,i-1] + r*St[:,i-1]*self.__stepsize + sigma*St[:,i-1]*dW
        elif self.__model == 'JD':
            #TODO: implement
            pass
        elif self.__model == 'Heston':
            #TODO: implement
            pass
        elif self.__model == 'OU':
            #TODO: implement
            pass
        return St
    
    def __exact(self, antithetic):
        S0 = self.__option.underlyingPrice
        r = self.__r
        T = self.__option.timeToExpiry
        sigma = self.__sigma
        
        numSteps = int(T/self.__stepsize)
        St = np.zeros((self.__paths, numSteps+1))
        if self.__model == 'GBM':
            dZ = np.arange(0, T+self.__stepsize, self.__stepsize)*np.random.standard_normal((self.__paths, numSteps+1))
            dt = np.arange(0, T+self.__stepsize, self.__stepsize)
            St = S0*np.exp((r - 0.5*(sigma**2))*dt + sigma*dZ)
        return St
    
    
    def optionPrice(self):
        K = self.__option.strikePrice
        r = self.__r
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
    
    