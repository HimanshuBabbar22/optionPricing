#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 15 21:10:48 2019

@author: himanshu
"""

import numpy as np
from monteCarloPricing import MonteCarloPricing
from option import Option

class BarrierOptionPricing:
    
    def __init__(self, 
                 option : Option,
                 riskFreeRate,
                 volatility,
                 stepsize=0.001, 
                 numOfPaths=1000, 
                 discretizationMethod = 'Euler'
                 ):
        
        self.__option = option
        self.__r = riskFreeRate
        self.__sigma = volatility
        self.__dt = stepsize
        self.__paths = numOfPaths
        self.__mcObj = MonteCarloPricing(option,
                                         riskFreeRate,
                                         volatility,
                                         stepsize,
                                         numOfPaths,
                                         discretizationMethod)
    
    def __upAndOutPayoff(self, St):
        K = self.__option.strikePrice
        T = self.__option.timeToExpiry
        numSteps = int(T/self.__dt)
        payoff = np.ones(self.__paths)
        
        level = self.__option.barrierLevel
        for i in range(1, numSteps+1):
            worthLessPaths = np.where((St[:,i] >= level)==True)[0]
            if len(worthLessPaths) > 0:
                payoff[worthLessPaths] = 0.0
        
        valuablePaths = np.where((payoff > 0.0)==True)[0]
        if self.__option.optionType == 'Call':
            payoff[valuablePaths] = St[valuablePaths,-1] - K
        else:
            payoff[valuablePaths] = K - St[valuablePaths,-1]
        
        payoff[payoff <= 0.0] = 0.0
        return payoff
    
    def __upAndInPayoff(self, St):
        K = self.__option.strikePrice
        T = self.__option.timeToExpiry
        numSteps = int(T/self.__dt)
        payoff = np.zeros(self.__paths)
        
        level = self.__option.barrierLevel
        paths = set()
        for i in range(1, numSteps+1):
            valuablePaths = np.where((St[:,i] >= level)==True)[0]
            paths.union(valuablePaths)
        
        if self.__option.optionType == 'Call':
            payoff[valuablePaths] = St[valuablePaths,-1] - K
        else:
            payoff[valuablePaths] = K - St[valuablePaths,-1]
        
        payoff[payoff <= 0.0] = 0.0
        return payoff
    
    def optionPrice(self):
        StPos, StNeg = self.__mcObj.samplePaths()
        
        if self.__option.barrierType == 'UO':
            payoffPos = self.__upAndOutPayoff(StPos)
            payoffNeg = self.__upAndOutPayoff(StNeg)
        elif self.__option.barrierType == 'UI':
            payoffPos = self.__upAndInPayoff(StPos)
            payoffNeg = self.__upAndInPayoff(StNeg)
        
        
        optVal = np.exp(-r*(T))*np.mean((payoffPos+payoffNeg)/2)
        return float("{0:.2f}".format(optVal))
 

if __name__ == '__main__':

    from option import Barrier
    S0 = 100
    K = 110
    r = 0.10
    T = 1
    volatility = 0.25
    
    print('------------------------------------------------------------------'
          +'----------------------------')
    option = Barrier(S0, K, T, 'Put', 'UO', 135.0)
    baPricing = BarrierOptionPricing(option, r, volatility, 0.001, 1000, discretizationMethod='Euler')
    #print(baPricing)
    print('Monte-Carlo Estimate for Up and Out Call:', baPricing.optionPrice())
    
# =============================================================================
#     print('------------------------------------------------------------------'
#           +'----------------------------')   
# 
#     option = Barrier(S0, K, T, 'CAll', 'UI', 135.0)
#     baPricing = BarrierOptionPricing(option, r, volatility, 0.001, 1000, discretizationMethod='Euler')
#     #print(baPricing)
#     print('Monte-Carlo Estimate for Up and In Call:', baPricing.optionPrice())    
# =============================================================================
        
        