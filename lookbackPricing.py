#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 16 02:26:53 2019

@author: himanshu
"""

import numpy as np
from monteCarloPricing import MonteCarloPricing
from option import Option

class LookbackOptionPricing:
    
    def __init__(self,
                 option:Option,
                 riskFreeRate,
                 volatility,
                 stepsize=0.001, 
                 numOfPaths=1000, 
                 discretizationMethod = 'Euler'):
        
        self.__option = option
        self.__r = riskFreeRate
        self.__sigma = volatility
        self.__mcObj = MonteCarloPricing(option,
                                         riskFreeRate,
                                         volatility,
                                         stepsize,
                                         numOfPaths,
                                         discretizationMethod)

    def __payoffFunc(self, St):
        K = self.__option.strikePrice
        M = np.max(St, axis=1)
        if self.__option.optionType == 'Call':
            payoff = M - K
        else:
            payoff = K - M
            
        payoff[payoff < 0.0] = 0.0
        return payoff
    
    def optionPrice(self):
        T = self.__option.timeToExpiry
        StPos, StNeg = self.__mcObj.samplePaths()
        
        if self.__option.optionType == 'Call':
            payoffPos = self.__payoffFunc(StPos)
            payoffNeg = self.__payoffFunc(StNeg)
        else:
            payoffPos = self.__payoffFunc(StPos)
            payoffNeg = self.__payoffFunc(StNeg)
            
        optVal = np.exp(-self.__r*(T))*np.mean((payoffPos+payoffNeg)/2)
        return float("{0:.2f}".format(optVal))
    
    def __repr__(self):
        
        return "LookbackOptionPricing({}, {}, {})"\
                .format(self.__option,
                        self.__r,
                        self.__sigma)
        
if __name__ == '__main__':

    from option import Lookback
    S0 = 100
    K = 110
    r = 0.10
    T = 1
    volatility = 0.25
    
    print('------------------------------------------------------------------'
          +'----------------------------')
    option = Lookback(S0, K, T, 'Call')
    lbPricing = LookbackOptionPricing(option, r, volatility, 0.001, 1000, discretizationMethod='Euler')
    print(lbPricing)
    print('Monte-Carlo Estimate for Lookback Call:', lbPricing.optionPrice())
    
    print('------------------------------------------------------------------'
          +'----------------------------')
    option = Lookback(S0, K, T, 'Put')
    lbPricing = LookbackOptionPricing(option, r, volatility, 0.001, 1000, discretizationMethod='Euler')
    print(lbPricing)
    print('Monte-Carlo Estimate for Lookback Put:', lbPricing.optionPrice())