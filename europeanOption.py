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
        
        self.__underlyingPrice = underlyingPrice
        self.__strikePrice = strikePrice
        self.__timeToExpiry = timeToExpiry
        self.__optionType = optionType
        self.__riskFreeRate = riskFreeRate
    
    def __d1(self, volatility, dividendYield=0.0):
        S0 = self.__underlyingPrice
        K = self.__strikePrice
        r = self.__riskFreeRate
        T = self.__timeToExpiry
        
        d1 = (np.log(S0/K) + (r-dividendYield+((volatility**2)/2))*(T))/(volatility*np.sqrt(T))
        return d1

    def blackScholesPrice(self, volatility, dividendYield=0.0):
        S0 = self.__underlyingPrice
        K = self.__strikePrice
        r = self.__riskFreeRate
        T = self.__timeToExpiry
        
        d1 = self.__d1(volatility, dividendYield)
        d2 = d1 - (volatility*np.sqrt(T))
    
        bsPrice = None
        if self.__optionType == 'Call':
            bsPrice = S0*np.exp(-dividendYield*T)*norm.cdf(d1) - K*np.exp(-r*T)*norm.cdf(d2)
        elif self.__optionType == 'Put':
            bsPrice = -S0*np.exp(-dividendYield*T)*norm.cdf(-d1) + K*np.exp(-r*T)*norm.cdf(-d2)
        if type(bsPrice) == np.ndarray:
            bsPrice = bsPrice[0]
        return "{0:.2f}".format(bsPrice)
    
    def blackScholesDelta(self, volatility, dividendYield=0.0):
        d1 = self.__d1(volatility, dividendYield)
        if self.__optionType == 'Call':
            return "{0:.6f}".format(np.exp(-dividendYield*T)*norm.cdf(d1))
        elif self.__optionType == 'Put':
            return "{0:.6f}".format(np.exp(-dividendYield*T)*(norm.cdf(d1)-1))
    
    def blackScholesGamma(self, volatility, dividendYield=0.0):
        S0 = self.__underlyingPrice
        T = self.__timeToExpiry
        
        d1 = self.__d1(volatility, dividendYield)
        gamma = norm.pdf(d1)*np.exp(-dividendYield*T)/(S0*volatility*np.sqrt(T))
        return "{0:.6f}".format(gamma)
        
    def blackScholesVega(self, volatility, dividendYield=0.0):
        S0 = self.__underlyingPrice
        T = self.__timeToExpiry
        
        d1 = self.__d1(volatility, dividendYield)
        vega = S0*np.sqrt(T)*norm.pdf(d1)*np.exp(-dividendYield*T)
        return "{0:.6f}".format(vega)
    
    def blackScholesTheta(self, volatility, dividendYield=0.0):
        S0 = self.__underlyingPrice
        K = self.__strikePrice
        r = self.__riskFreeRate
        T = self.__timeToExpiry
        
        d1 = self.__d1(volatility, dividendYield)
        d2 = d1 - (volatility*np.sqrt(T))
        
        s1 = -S0*norm.pdf(d1)*volatility*np.exp(-dividendYield*T)/(2*np.sqrt(T))
        if self.__optionType == 'Call':
            s2 = dividendYield*S0*norm.cdf(d1)*np.exp(-dividendYield*T)\
                 - r*K*np.exp(-r*T)*norm.cdf(d2)
        elif self.__optionType == 'Put':
            s2 = -dividendYield*S0*norm.cdf(-d1)*np.exp(-dividendYield*T)\
                 + r*K*np.exp(-r*T)*norm.cdf(-d2)
        
        theta = s1+s2
        return "{0:.6f}".format(theta)
    
    def blackScholesRho(self, volatility, dividendYield=0.0):
        K = self.__strikePrice
        r = self.__riskFreeRate
        T = self.__timeToExpiry
        
        d1 = self.__d1(volatility, dividendYield)
        d2 = d1 - (volatility*np.sqrt(T))
        
        if self.__optionType == 'Call':
            return "{0:.6f}".format(K*T*np.exp(-r*T)*norm.cdf(d2))
        elif self.__optionType == 'Put':
            return "{0:.6f}".format(-K*T*np.exp(-r*T)*norm.cdf(-d2))
    
    @staticmethod
    def __errorFunc(implVol, *args):
        obj, price, dividendYield = args
        err = abs(float(EuropeanOption.blackScholesPrice(obj, implVol, dividendYield)) - price)
        return err
    
    def blackScholesImpliedVol(self, price, dividendYield=0.0):
        args = (self, price, dividendYield)
        impVol = sop.fmin(EuropeanOption.__errorFunc, 1.0, args=args, disp=0)
        return "{0:.6f}".format(impVol[0])
    
    def __samplePaths(self, volatility, stepsize, paths, method='Euler'):
        S0 = self.__underlyingPrice
        r = self.__riskFreeRate
        T = self.__timeToExpiry
        
        # Set random seed
        np.random.seed(seed = 1231)
        num_steps = int(T/stepsize)
        
        St = np.zeros((paths, num_steps+1))
        St[:,0] = S0
        
        if method == 'Euler':
            # Generate Stock prices using Euler discretization
            for i in range(1, num_steps+1):
                dZ = np.random.standard_normal(paths)
                dW = np.sqrt(stepsize)*dZ
                St[:,i] = St[:,i-1] + r*St[:,i-1]*stepsize + volatility*St[:,i-1]*dW
        else:
            #TODO: implement Exact discretization
            pass
    
        return St
    
    def monteCarloEstimate(self, volatility, stepsize, paths, method='Euler'):
        K = self.__strikePrice
        r = self.__riskFreeRate
        T = self.__timeToExpiry
        
        St = self.__samplePaths(volatility, stepsize, paths)
        if self.__optionType == 'Call':
            payoff = St[:,-1] - K
        elif self.__optionType == 'Put':
            payoff = K - St[:,-1]
            
        payoff[payoff < 0.0] = 0.0
        mc_val = np.exp(-r*(T))*np.mean(payoff)
        
        return "{0:.2f}".format(mc_val)
    
    def __repr__(self):
        return "Option('{}', '{}', '{}', '{}', '{}')"\
              .format(self.__underlyingPrice,
                      self.__strikePrice,
                      self.__timeToExpiry,
                      self.__optionType,
                      self.__riskFreeRate)
        
    def __str__(self):
        return "European {} Option; Underlying Price = {}; Strike Price = {}; Time to Expiry = {} year(s)"\
              .format(self.__optionType,
                      self.__underlyingPrice, 
                      self.__strikePrice,
                      self.__timeToExpiry)
        
if __name__ == "__main__":
    S0 = 100
    K = 110
    r = 0.10
    T = 1
    volatility = 0.25
    
    print('------------------------------------------------------------------'
          +'----------------------------')
    option = EuropeanOption(S0, K, T, 'Call', r)
    print(option)
    print('Monte-Carlo Estimate:', option.monteCarloEstimate(volatility, 0.001, 1000))
    print('Black Scholes Price:', option.blackScholesPrice(volatility, 0.0))
    print('Black Scholes Delta:', option.blackScholesDelta(volatility, 0.0))
    print('Black Scholes Gamma:', option.blackScholesGamma(volatility, 0.0))
    print('Black Scholes Vega:', option.blackScholesVega(volatility, 0.0))
    print('Black Scholes Theta:', option.blackScholesTheta(volatility, 0.0))
    print('Black Scholes Rho:', option.blackScholesRho(volatility, 0.0))
    print('Black Scholes Implied-Vol:', option.blackScholesImpliedVol(12.0, 0.0))
    
    print('------------------------------------------------------------------'
          +'----------------------------')
    
    
    option = EuropeanOption(S0, K, T, 'Put', r)
    print(option)
    print('Monte-Carlo Estimate:', option.monteCarloEstimate(volatility, 0.001, 1000))
    print('Black Scholes Price:', option.blackScholesPrice(volatility, 0.0))
    print('Black Scholes Delta:', option.blackScholesDelta(volatility, 0.0))
    print('Black Scholes Gamma:', option.blackScholesGamma(volatility, 0.0))
    print('Black Scholes Vega:', option.blackScholesVega(volatility, 0.0))
    print('Black Scholes Theta:', option.blackScholesTheta(volatility, 0.0))
    print('Black Scholes Rho:', option.blackScholesRho(volatility, 0.0))
    print('Black Scholes Implied-Vol:', option.blackScholesImpliedVol(12.0, 0.0))
    print('------------------------------------------------------------------'
          +'----------------------------')