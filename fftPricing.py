#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 14 01:32:26 2019

@author: himanshu
"""

import numpy as np
from scipy.interpolate import interp1d
from option import Option

class FFTPricing:
    
    def __init__(self,
                 option : Option,
                 riskFreeRate,
                 volatility,
                 samplePoints,
                 bandwidth,
                 dampingFactor,
                 underlyingModel = 'GBM'):
        
        self.__option = option
        self.__r = riskFreeRate
        self.__sigma = volatility
        self.__N = samplePoints
        self.__B = bandwidth
        self.__alpha = dampingFactor
        self.__model = underlyingModel
        
    
    # Computes the characterstic function of a GBM.
    def __charactersticFunc(self, omega):
        S0 = self.__option.underlyingPrice
        r = self.__r
        T = self.__option.timeToExpiry
        sigma = self.__sigma
        alpha = self.__alpha
        
        if self.__model == 'GBM':
            x0 = np.log(S0)
            mu = x0 + ((r - (sigma**2)/2)*(T))
            sig = (sigma**2)*(T)/2
            omega_prime = omega + 1j*(alpha+1)
            return np.exp(-1j*mu*omega_prime - sig*(omega_prime**2))
        elif self.__model == 'VG':
            pass
    
    # Computes the Fourier Transform of a GBM.
    def __fourierTransform(self, omega):
        alpha = self.__alpha
        r = self.__r
        T = self.__option.timeToExpiry
        
        q_hat = self.__charactersticFunc(omega)
        num = np.exp(-r*(T))*q_hat
        den = (alpha - 1j*omega)*(alpha - (1j*omega) + 1)
        return num/den
    
    def optionPrice(self):
        if not self.__option.expiryType == 'European':
            print('Not a European Option')
            return 0.0
        
        K = self.__option.strikePrice
        
        N = self.__N
        B = self.__B
        alpha = self.__alpha
        
        h = B/(N-1)
        omega = np.arange(0,N)*h
        
        dk = 2*np.pi/(h*N)
        k = np.log(20) + np.arange(0,N)*dk
        
        dw = np.zeros(N)
        dw[0] = h/2
        dw[1:] = h
        
        # FFT Algorithm
        V = np.zeros(N)
        for n in range(N):
            nu_hat = self.__fourierTransform(omega)
            inner_sum = np.sum(np.exp(1j*omega*k[n])*nu_hat*dw)
            V[n] = ((np.exp(-alpha*k[n])/np.pi)*inner_sum).real
        
        val = interp1d(k, V)
        return float('{0:.2f}'.format(val(np.log(K))))
    
    def __repr__(self):
        
        return "FFTPricing({}, {}, {}, {}, {}, {})"\
                .format(self.__option,
                        self.__r,
                        self.__sigma,
                        self.__N,
                        self.__B,
                        self.__alpha)
    
if __name__ == "__main__":
    from option import European
    S0 = 100
    K = 110
    r = 0.10
    T = 1
    volatility = 0.25
    
    N = 2**10
    B = 50
    alpha = 10.0
    
    print('------------------------------------------------------------------'
          +'----------------------------')
    option = European(S0, K, T, 'Call')
    fftPricing = FFTPricing(option, r, volatility, N, B, alpha)
    print(fftPricing)
    print('FFT price for Call:', fftPricing.optionPrice())
    
    print('------------------------------------------------------------------'
          +'----------------------------')
    option = European(S0, K, T, 'Put')
    fftPricing = FFTPricing(option, r, volatility, N, B, -alpha)
    print(fftPricing)
    print('FFT price for Put:', fftPricing.optionPrice())
    