#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 12 17:33:39 2019

@author: himanshu
"""

class Option:
    
    def __init__(self,
                 underlyingPrice, 
                 strikePrice, 
                 timeToExpiry,
                 optionType,
                 expiryType):
        
        self.underlyingPrice = float(underlyingPrice)
        self.strikePrice = float(strikePrice)
        self.timeToExpiry = float(timeToExpiry)
        self.optionType = optionType
        self.expiryType = expiryType
    
    def __repr__(self):
        
        return "Option({}, {}, {}, {}, {})"\
              .format(self.underlyingPrice,
                      self.strikePrice,
                      self.timeToExpiry,
                      self.optionType,
                      self.expiryType)
        
        
if __name__ == "__main__":
    S0 = 100
    K = 110
    T = 1
    volatility = 0.25
    
    option = Option(S0, K, T, 'Call', 'European')
    print(option)
    