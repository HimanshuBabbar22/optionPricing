#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 13 11:37:44 2019

@author: himanshu
"""

import numpy as np
from scipy.optimize import fmin
from scipy.stats import norm
from europeanOption import EuropeanOption

class BlackScholesPricing:
    
    def __init__(self, option, volatility, dividendYield=0.0):
        self.__option = option
        self.__sigma = volatility
        self.__div = dividendYield
        pass
    
    def __d1(self):
        S0 = self.__option.underlyingPrice
        K = self.__option.strikePrice
        r = self.__option.riskFreeRate
        T = self.__option.timeToExpiry
        sigma = self.__sigma
        q = self.__div
        
        d1 = (np.log(S0/K) + (r-q+((sigma**2)/2))*(T))/(sigma*np.sqrt(T))
        
        return d1

    def optionPrice(self, sigma=None):
        S0 = self.__option.underlyingPrice
        K = self.__option.strikePrice
        r = self.__option.riskFreeRate
        T = self.__option.timeToExpiry
        if sigma is None: sigma = self.__sigma
        q = self.__div
        
        d1 = self.__d1()
        d2 = d1 - (sigma*np.sqrt(T))
    
        bsPrice = None
        if self.__option.optionType == 'Call':
            bsPrice = S0*np.exp(-q*T)*norm.cdf(d1) - K*np.exp(-r*T)*norm.cdf(d2)
        elif self.__option.optionType == 'Put':
            bsPrice = -S0*np.exp(-q*T)*norm.cdf(-d1) + K*np.exp(-r*T)*norm.cdf(-d2)
        if type(bsPrice) == np.ndarray:
            bsPrice = bsPrice[0]
        
        return "{0:.2f}".format(bsPrice)
    
    def delta(self):
        T = self.__option.timeToExpiry
        q = self.__div
        
        d1 = self.__d1()
        if self.__option.optionType == 'Call':
            return "{0:.6f}".format(np.exp(-q*T)*norm.cdf(d1))
        elif self.__option.optionType == 'Put':
            return "{0:.6f}".format(np.exp(-q*T)*(norm.cdf(d1)-1))
    
    def gamma(self):
        S0 = self.__option.underlyingPrice
        T = self.__option.timeToExpiry
        sigma = self.__sigma
        q = self.__div
        
        d1 = self.__d1()
        gamma = norm.pdf(d1)*np.exp(-q*T)/(S0*sigma*np.sqrt(T))
        
        return "{0:.6f}".format(gamma)
        
    def vega(self):
        S0 = self.__option.underlyingPrice
        T = self.__option.timeToExpiry
        q = self.__div
        
        d1 = self.__d1()
        vega = S0*np.sqrt(T)*norm.pdf(d1)*np.exp(-q*T)
        
        return "{0:.6f}".format(vega)
    
    def theta(self):
        S0 = self.__option.underlyingPrice
        K = self.__option.strikePrice
        r = self.__option.riskFreeRate
        T = self.__option.timeToExpiry
        sigma = self.__sigma
        q = self.__div
        
        d1 = self.__d1()
        d2 = d1 - (sigma*np.sqrt(T))
        
        s1 = -S0*norm.pdf(d1)*sigma*np.exp(-q*T)/(2*np.sqrt(T))
        if self.__option.optionType == 'Call':
            s2 = q*S0*norm.cdf(d1)*np.exp(-q*T)\
                 - r*K*np.exp(-r*T)*norm.cdf(d2)
        elif self.__option.optionType == 'Put':
            s2 = -q*S0*norm.cdf(-d1)*np.exp(-q*T)\
                 + r*K*np.exp(-r*T)*norm.cdf(-d2)
        
        theta = s1+s2
        
        return "{0:.6f}".format(theta)
    
    def rho(self):
        K = self.__option.strikePrice
        r = self.__option.riskFreeRate
        T = self.__option.timeToExpiry
        sigma = self.__sigma
        
        d1 = self.__d1()
        d2 = d1 - (sigma*np.sqrt(T))
        
        if self.__option.optionType == 'Call':
            return "{0:.6f}".format(K*T*np.exp(-r*T)*norm.cdf(d2))
        elif self.__option.optionType == 'Put':
            return "{0:.6f}".format(-K*T*np.exp(-r*T)*norm.cdf(-d2))
    
    @classmethod
    def __errorFunc(cls, implVol, *args):
        obj, price, dividendYield = args
        err = abs(float(cls.optionPrice(obj, implVol)) - price)
        return err
    
    def impliedVolatility(self, price, dividendYield=0.0):
        args = (self, price, dividendYield)
        impVol = fmin(BlackScholesPricing.__errorFunc, 1.0, args=args, disp=0)
        return "{0:.6f}".format(impVol[0])
    
    def __repr__(self):
        
        return "BlackScholesPricing({}, {}, {})"\
                .format(self.__option,
                        self.__sigma,
                        self.__div)
    
if __name__ == "__main__":
    S0 = 100
    K = 110
    r = 0.10
    T = 1
    volatility = 0.25
    
    print('------------------------------------------------------------------'
          +'----------------------------')
    option = EuropeanOption(S0, K, T, 'Call', r)
    bsPricing = BlackScholesPricing(option, volatility)
    print(bsPricing)
    print('Black Scholes Price:', bsPricing.optionPrice())
    print('Black Scholes Delta:', bsPricing.delta())
    print('Black Scholes Gamma:', bsPricing.gamma())
    print('Black Scholes Vega:', bsPricing.vega())
    print('Black Scholes Theta:', bsPricing.theta())
    print('Black Scholes Rho:', bsPricing.rho())
    print('Black Scholes Implied-Vol:', bsPricing.impliedVolatility(12.0))
    
    print('------------------------------------------------------------------'
          +'----------------------------')
    
    
    option = EuropeanOption(S0, K, T, 'Call', r)
    bsPricing = BlackScholesPricing(option, volatility)
    print(bsPricing)
    print('Black Scholes Price:', bsPricing.optionPrice())
    print('Black Scholes Delta:', bsPricing.delta())
    print('Black Scholes Gamma:', bsPricing.gamma())
    print('Black Scholes Vega:', bsPricing.vega())
    print('Black Scholes Theta:', bsPricing.theta())
    print('Black Scholes Rho:', bsPricing.rho())
    print('Black Scholes Implied-Vol:', bsPricing.impliedVolatility(13.0))
    print('------------------------------------------------------------------'
          +'----------------------------')