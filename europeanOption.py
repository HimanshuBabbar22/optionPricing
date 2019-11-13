#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 12 17:33:39 2019

@author: himanshu
"""

import numpy as np
import scipy.optimize as sop
from scipy.stats import norm


class EuropeanOption:
    
    def __init__(self,
                 underlyingPrice, 
                 strikePrice, 
                 timeToExpiry,
                 optionType,
                 riskFreeRate):
        
        self.underlyingPrice = int(underlyingPrice)
        self.strikePrice = int(strikePrice)
        self.timeToExpiry = int(timeToExpiry)
        self.optionType = optionType
        self.riskFreeRate = int(riskFreeRate)
        
        
    def __repr__(self):
        optType = self.optionType
        if optType is not None:
            return "EuropeanOption({}, {}, {}, {}, {})"\
                  .format(self.underlyingPrice,
                          self.strikePrice,
                          self.timeToExpiry,
                          optType,
                          self.riskFreeRate)
        
        
if __name__ == "__main__":
    S0 = 100
    K = 110
    r = 0.10
    T = 1
    volatility = 0.25
    
    option = EuropeanOption(S0, K, T, 'Call', r)
    print(option)
    