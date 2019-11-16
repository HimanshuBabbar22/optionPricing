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
                 timeToExpiry):
        
        self.underlyingPrice = float(underlyingPrice)
        self.strikePrice = float(strikePrice)
        self.timeToExpiry = float(timeToExpiry)
        
    def __repr__(self):
        
        return "Option({}, {}, {}, {}, {})"\
              .format(self.underlyingPrice,
                      self.strikePrice,
                      self.timeToExpiry)
        
class European(Option):
    
    def __init__(self,
                 underlyingPrice, 
                 strikePrice, 
                 timeToExpiry,
                 optionType):
        
        super().__init__(underlyingPrice, strikePrice, timeToExpiry)
        self.optionType = optionType
        self.expiryType = 'European'
        
    def __repr__(self):
        
        return "European({}, {}, {}, {})"\
              .format(self.underlyingPrice,
                      self.strikePrice,
                      self.timeToExpiry,
                      self.optionType)
              
class American(Option):
    
    def __init__(self,
                 underlyingPrice, 
                 strikePrice, 
                 timeToExpiry,
                 optionType):
        
        super().__init__(underlyingPrice, strikePrice, timeToExpiry)
        self.optionType = optionType
        self.expiryType = 'American'
        
    def __repr__(self):
        
        return "American({}, {}, {}, {})"\
              .format(self.underlyingPrice,
                      self.strikePrice,
                      self.timeToExpiry,
                      self.optionType)

class Barrier(Option):
    
    def __init__(self,
                 underlyingPrice, 
                 strikePrice, 
                 timeToExpiry,
                 optionType,
                 barrierType,
                 barrierLevel):
        
        super().__init__(underlyingPrice, strikePrice, timeToExpiry)
        self.optionType = optionType
        self.barrierType = barrierType
        self.barrierLevel = barrierLevel
        self.expiryType = 'Barrier'
        
    def __repr__(self):
        
        return "Barrier({}, {}, {}, {}, {}, {})"\
              .format(self.underlyingPrice,
                      self.strikePrice,
                      self.timeToExpiry,
                      self.optionType,
                      self.barrierType,
                      self.barrierLevel)
    
    
if __name__ == "__main__":
    S0 = 100
    K = 110
    T = 1
    volatility = 0.25
    
    euOption = European(S0, K, T, 'Call')
    print(euOption)
    
    amOption = American(S0, K, T, 'Call')
    print(amOption)
    
    baOption = Barrier(S0, K, T, 'Call', 'UO', 125)
    print(baOption)