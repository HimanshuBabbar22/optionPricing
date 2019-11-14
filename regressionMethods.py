#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 13 18:49:43 2019

@author: himanshu
"""

import numpy as np
from scipy.optimize import curve_fit as curveFit
from monteCarloPricing import MonteCarloPricing

class RegressionMethods:
    
    def __init__(self, 
                 option, 
                 riskFreeRate,
                 volatility,
                 stepsize=0.001, 
                 numOfPaths=1000, 
                 underlyingModel = 'GBM', 
                 discretizationMethod = 'Euler'
                 ):
        
        self.__option = option
        self.__r = riskFreeRate
        self.__sigma = volatility
        self.__stepsize = stepsize
        self.__mcObj = MonteCarloPricing(option,
                                         riskFreeRate,
                                         volatility,
                                         stepsize,
                                         numOfPaths,
                                         underlyingModel,
                                         discretizationMethod)
    

    @staticmethod
    def __cHat(x, a0, a1, a2, a3):
        return a0 + a1*x + a2*(x**2) + a3*(x**3)

    def __payoffFunc(self, S):
        K = self.__option.strikePrice
        if self.__option.optionType == 'Call':
            payoff = S-K
        else:
            payoff = K-S
        payoff[payoff < 0] = 0.0
        return payoff
    
    def optionPrice(self, method=1):
        if method == 1:
            return self.__regMethod1()
        return self.__regMethod2()
        
    def __regMethod1(self):
        if not self.__option.expiryType == 'American':
            print('Not an American Option')
            return 0.0
        
        stepsize = self.__stepsize
        
        St = self.__mcObj.samplePaths()
        M = St.shape[1]-1
        V = np.zeros(St.shape)
        V[:,M] = self.__payoffFunc(St[:,M])
        
        for i in range(M-1, 0, -1):
            [opt_params, cov] = curveFit(RegressionMethods.__cHat, St[:,i], np.exp(-self.__r*stepsize)*V[:,i+1])
            actual_payoff = self.__payoffFunc(St[:,i])
            fitted_payoff = RegressionMethods.__cHat(St[:,i], *opt_params)
            V[:,i] = np.maximum(actual_payoff, fitted_payoff)
        
        optPrice = np.exp(-self.__r*stepsize)*np.mean(V[:,1])
        return float('{0:.2f}'.format(optPrice))
    
    def __regMethod2(self):
        if not self.__option.expiryType == 'American':
            print('Not an American Option')
            return 0.0
        
        stepsize = self.__stepsize
        St = self.__mcObj.samplePaths()
        M = St.shape[1]-1
        N = St.shape[0]
        g = self.__payoffFunc(St[:,M])
        tao = np.ones(N)*M
        
        for i in range(M-1, 0, -1):
            if self.__option.optionType == 'Call':
                k = np.where(St[:,i] > K)[0]
            else:
                k = np.where(K > St[:,i])[0]
            if len(k) > 3:
                [opt_params, cov] = curveFit(RegressionMethods.__cHat, 
                                                St[k,i], 
                                                np.exp(-self.__r*(tao[k]-i)*stepsize)*g[k]
                                            )
            actual_payoff = self.__payoffFunc(St[k,i])
            fitted_payoff = RegressionMethods.__cHat(St[k,i], *opt_params)
            update_ixs = np.where((actual_payoff >= fitted_payoff) == True)[0]
            g[k[update_ixs]] = actual_payoff[update_ixs]
            tao[k[update_ixs]] = i
            
        cHat0 = np.mean(np.exp(-self.__r*tao*stepsize)*g)
        
        v0 = np.mean(self.__payoffFunc(St[:,0]))
        
        optPrice = max(v0, cHat0)
        return float('{0:.2f}'.format(optPrice))
    
    def __repr__(self):
        
        return "RegressionMethods({}, {}, {})"\
                .format(self.__option,
                        self.__r,
                        self.__sigma)
        
if __name__ == "__main__":
    from option import Option
    S0 = 100
    K = 110
    r = 0.10
    T = 1
    volatility = 0.25
    
    print('------------------------------------------------------------------'
          +'----------------------------')
    option = Option(S0, K, T, 'Call', 'American')
    regPricing = RegressionMethods(option, r, volatility)
    print(regPricing)
    print('Regression Method 1 price for Call:', regPricing.optionPrice())
    
    print('------------------------------------------------------------------'
          +'----------------------------')
    
    option = Option(S0, K, T,  'Put', 'American')
    regPricing = RegressionMethods(option, r, volatility)
    print(regPricing)
    print('Regression Method 1 price for Put:', regPricing.optionPrice())
    
    print('------------------------------------------------------------------'
          +'----------------------------')
    option = Option(S0, K, T, 'Call', 'American')
    regPricing = RegressionMethods(option, r, volatility)
    print(regPricing)
    print('Regression Method 2 price for Call:', regPricing.optionPrice(2))
    
    print('------------------------------------------------------------------'
          +'----------------------------')
    
    option = Option(S0, K, T,  'Put', 'American')
    regPricing = RegressionMethods(option, r, volatility)
    print(regPricing)
    print('Regression Method 2 price for Put:', regPricing.optionPrice(2))
    