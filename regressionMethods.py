#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 13 18:49:43 2019

@author: himanshu
"""

import numpy as np
from scipy.optimize import curve_fit as curveFit
from monteCarloPricing import MonteCarloPricing
from option import Option

class RegressionMethods:
    
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
        self.__stepsize = stepsize
        self.__mcObj = MonteCarloPricing(option,
                                         riskFreeRate,
                                         volatility,
                                         stepsize,
                                         numOfPaths,
                                         discretizationMethod)
    

    @staticmethod
    def __cHat(x, a0, a1, a2, a3):
        return a0 + a1*x + a2*(x**2) + a3*(x**3)

    def __payoffFunc(self, SPos, SNeg):
        K = self.__option.strikePrice
        if self.__option.optionType == 'Call':
            payoffPos = SPos-K
            payoffNeg = SNeg-K
        else:
            payoffPos = K-SPos
            payoffNeg = K-SNeg
        payoffPos[payoffPos < 0] = 0.0
        payoffNeg[payoffNeg < 0] = 0.0
        return payoffPos, payoffNeg
    
    def optionPrice(self, method=1):
        if method == 1:
            return self.__regMethod1()
        return self.__regMethod2()
        
    def __regMethod1(self):
        if not self.__option.expiryType == 'American':
            print('Not an American Option')
            return 0.0
        
        stepsize = self.__stepsize
        
        StPos, StNeg = self.__mcObj.samplePaths()
        M = StPos.shape[1]-1
        VPos = np.zeros(StPos.shape)
        VNeg = np.zeros(StPos.shape)
        
        VPos[:,M], VNeg[:,M] = self.__payoffFunc(StPos[:,M], StNeg[:,M])
        
        for i in range(M-1, 0, -1):
            [optParamsPos, cov] = curveFit(RegressionMethods.__cHat, StPos[:,i],\
                                            np.exp(-self.__r*stepsize)*VPos[:,i+1])
            
            [optParamsNeg, cov] = curveFit(RegressionMethods.__cHat, StNeg[:,i],\
                                            np.exp(-self.__r*stepsize)*VNeg[:,i+1])
            
            actPayoffPos, actPayoffNeg = self.__payoffFunc(StPos[:,i], StNeg[:,i])
            fitPayoffPos = RegressionMethods.__cHat(StPos[:,i], *optParamsPos)
            fitPayoffNeg = RegressionMethods.__cHat(StNeg[:,i], *optParamsNeg)
            
            VPos[:,i] = np.maximum(actPayoffPos, fitPayoffPos)
            VNeg[:,i] = np.maximum(actPayoffNeg, fitPayoffNeg)
        
        optPrice = np.exp(-self.__r*stepsize)*np.mean((VPos[:,1]+VNeg[:,1])/2)
        return float('{0:.2f}'.format(optPrice))
    
    def __regMethod2(self):
        if not self.__option.expiryType == 'American':
            print('Not an American Option')
            return 0.0
        
        stepsize = self.__stepsize
        StPos, StNeg = self.__mcObj.samplePaths()
        M = StPos.shape[1]-1
        N = StPos.shape[0]
        gPos, gNeg = self.__payoffFunc(StPos[:,M], StNeg[:,M])
        taoPos, taoNeg = np.ones(N)*M, np.ones(N)*M
        
        for i in range(M-1, 0, -1):
            if self.__option.optionType == 'Call':
                kPos = np.where(StPos[:,i] > K)[0]
                kNeg = np.where(StNeg[:,i] > K)[0]
            else:
                kPos = np.where(K > StPos[:,i])[0]
                kNeg = np.where(K > StNeg[:,i])[0]
            if (len(kPos) > 3) and (len(kNeg) > 3):
                [optParamsPos, cov] = curveFit(RegressionMethods.__cHat, 
                                                StPos[kPos,i], 
                                                np.exp(-self.__r*(taoPos[kPos]-i)*stepsize)*gPos[kPos]
                                            )
                [optParamsNeg, cov] = curveFit(RegressionMethods.__cHat, 
                                                StNeg[kNeg,i], 
                                                np.exp(-self.__r*(taoNeg[kNeg]-i)*stepsize)*gNeg[kNeg]
                                            )
                
            actPayoffPos, actPayoffNeg = self.__payoffFunc(StPos[kPos,i], StNeg[kNeg,i])
            fitPayoffPos = RegressionMethods.__cHat(StPos[kPos,i], *optParamsPos)
            fitPayoffNeg = RegressionMethods.__cHat(StNeg[kNeg,i], *optParamsNeg)
            
            updateIdxsPos = np.where((actPayoffPos >= fitPayoffPos) == True)[0]
            updateIdxsNeg = np.where((actPayoffNeg >= fitPayoffNeg) == True)[0]
            gPos[kPos[updateIdxsPos]] = actPayoffPos[updateIdxsPos]
            gNeg[kNeg[updateIdxsNeg]] = actPayoffNeg[updateIdxsNeg]
            taoPos[kPos[updateIdxsPos]] = i
            taoNeg[kNeg[updateIdxsNeg]] = i
            
        cHatPos0 = np.exp(-self.__r*taoPos*stepsize)*gPos
        cHatNeg0 = np.exp(-self.__r*taoNeg*stepsize)*gNeg
        cHat0 = np.mean((cHatPos0+cHatNeg0)/2)
        
        vPos0, vNeg0 = self.__payoffFunc(StPos[:,0], StNeg[:,0])
        v0 = np.mean((vPos0+vNeg0)/2)
        
        optPrice = max(v0, cHat0)
        return float('{0:.2f}'.format(optPrice))
    
    def __repr__(self):
        
        return "RegressionMethods({}, {}, {})"\
                .format(self.__option,
                        self.__r,
                        self.__sigma)
        
if __name__ == "__main__":
    from option import American
    S0 = 100
    K = 110
    r = 0.10
    T = 1
    volatility = 0.25
    
    print('------------------------------------------------------------------'
          +'----------------------------')
    option = American(S0, K, T, 'Call')
    regPricing = RegressionMethods(option, r, volatility)
    print(regPricing)
    print('Regression Method 1 price for Call:', regPricing.optionPrice())
    
    print('------------------------------------------------------------------'
          +'----------------------------')
    
    option = American(S0, K, T, 'Put')
    regPricing = RegressionMethods(option, r, volatility)
    print(regPricing)
    print('Regression Method 1 price for Put:', regPricing.optionPrice())
    
    print('------------------------------------------------------------------'
          +'----------------------------')
    option = American(S0, K, T, 'Call')
    regPricing = RegressionMethods(option, r, volatility)
    print(regPricing)
    print('Regression Method 2 price for Call:', regPricing.optionPrice(2))
    
    print('------------------------------------------------------------------'
          +'----------------------------')
    
    option = American(S0, K, T, 'Put')
    regPricing = RegressionMethods(option, r, volatility)
    print(regPricing)
    print('Regression Method 2 price for Put:', regPricing.optionPrice(2))
    