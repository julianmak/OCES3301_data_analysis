# subroutines

import numpy as np

# dirty hack to get exponents and precision of a number
def get_exponent_precision(num, precision):
    out = []
    out.append(int(np.floor(np.log10(abs(num))))) # exponent
    out.append(round(num / float(10**int(np.floor(np.log10(abs(num))))), precision))
    
    return out

def sigmai_dep(CT, SA, p):
    """
    Compute the in-situ (or potential) density from CONSERVATIVE Temperature, 
    ABSOLUTE Salinity and pressure fields using the TEOS-10 EOS

    p can be a field (computed with subroutine p_from_z say) to get in-situ 
    density or can be a number p = p_ref, then it is a potential density 
    referenced to p_ref

    Adapted from the MATLAB GSW toolbox (http://www.TEOS-10.org)

    Inputs:
    CT             = Conservative Temperature      t        deg celsius
    SA             = Absolute Salinity             s        g / kg
    p              = (reference) pressure          p        dbar

    Returns:
    sigmai_dep_out = (potential density            rho      kg / m^3

    """
    # ensures that SA is non-negative.
    SA = abs(SA)

    # deltaS = 24
    sfac = 0.0248826675584615                 # sfac   = 1/(40*(35.16504/35)).
    offset = 5.971840214030754e-1             # offset = deltaS*sfac.

    x2 = sfac * SA
    xs = np.sqrt(x2 + offset)
    ys = CT * 0.025
    z  = p * 1e-4

    v000 =  1.0769995862e-3
    v001 = -6.0799143809e-5
    v002 =  9.9856169219e-6
    v003 = -1.1309361437e-6
    v004 =  1.0531153080e-7
    v005 = -1.2647261286e-8
    v006 =  1.9613503930e-9
    v010 = -1.5649734675e-5
    v011 =  1.8505765429e-5
    v012 = -1.1736386731e-6
    v013 = -3.6527006553e-7
    v014 =  3.1454099902e-7
    v020 =  2.7762106484e-5
    v021 = -1.1716606853e-5
    v022 =  2.1305028740e-6
    v023 =  2.8695905159e-7
    v030 = -1.6521159259e-5
    v031 =  7.9279656173e-6
    v032 = -4.6132540037e-7
    v040 =  6.9111322702e-6
    v041 = -3.4102187482e-6
    v042 = -6.3352916514e-8
    v050 = -8.0539615540e-7
    v051 =  5.0736766814e-7
    v060 =  2.0543094268e-7
    v100 = -3.1038981976e-4
    v101 =  2.4262468747e-5
    v102 = -5.8484432984e-7
    v103 =  3.6310188515e-7
    v104 = -1.1147125423e-7
    v110 =  3.5009599764e-5
    v111 = -9.5677088156e-6
    v112 = -5.5699154557e-6
    v113 = -2.7295696237e-7
    v120 = -3.7435842344e-5
    v121 = -2.3678308361e-7
    v122 =  3.9137387080e-7
    v130 =  2.4141479483e-5
    v131 = -3.4558773655e-6
    v132 =  7.7618888092e-9
    v140 = -8.7595873154e-6
    v141 =  1.2956717783e-6
    v150 = -3.3052758900e-7
    v200 =  6.6928067038e-4
    v201 = -3.4792460974e-5
    v202 = -4.8122251597e-6
    v203 =  1.6746303780e-8
    v210 = -4.3592678561e-5
    v211 =  1.1100834765e-5
    v212 =  5.4620748834e-6
    v220 =  3.5907822760e-5
    v221 =  2.9283346295e-6
    v222 = -6.5731104067e-7
    v230 = -1.4353633048e-5
    v231 =  3.1655306078e-7
    v240 =  4.3703680598e-6
    v300 = -8.5047933937e-4
    v301 =  3.7470777305e-5
    v302 =  4.9263106998e-6
    v310 =  3.4532461828e-5
    v311 = -9.8447117844e-6
    v312 = -1.3544185627e-6
    v320 = -1.8698584187e-5
    v321 = -4.8826139200e-7
    v330 =  2.2863324556e-6
    v400 =  5.8086069943e-4
    v401 = -1.7322218612e-5
    v402 = -1.7811974727e-6
    v410 = -1.1959409788e-5
    v411 =  2.5909225260e-6
    v420 =  3.8595339244e-6
    v500 = -2.1092370507e-4
    v501 =  3.0927427253e-6
    v510 =  1.3864594581e-6
    v600 =  3.1932457305e-5

    v = v000 + ( 
        xs * (v100 + xs * (v200 + xs * (v300 + xs * (v400 + xs * (v500 
      + v600 * xs))))) + ys * (v010 + xs * (v110 + xs * (v210 + xs * (v310 
      + xs * (v410 + v510 * xs)))) + ys * (v020 + xs * (v120 + xs * (v220 
      + xs * (v320 + v420 * xs))) + ys * (v030 + xs * (v130 + xs * (v230 
      + v330 * xs)) + ys * (v040 + xs * (v140 + v240*xs) + ys * (v050 
      + v150 * xs + v060 * ys))))) + z * (v001 + xs * (v101 + xs * (v201 
      + xs * (v301 + xs * (v401 + v501 * xs)))) + ys * (v011 + xs * (v111
      + xs * (v211 + xs * (v311 + v411 * xs))) + ys * (v021 + xs * (v121 
      + xs * (v221 + v321 * xs)) + ys * (v031 + xs * (v131 + v231 * xs) 
      + ys * (v041 + v141 * xs + v051 * ys)))) + z * (v002 + xs * (v102 
      + xs * (v202 + xs * (v302 + v402 * xs))) + ys * (v012 + xs * (v112 
      + xs * (v212 + v312 * xs)) + ys * (v022 + xs * (v122 + v222 * xs) 
      + ys * (v032 + v132 * xs + v042 * ys))) + z * (v003 + xs * (v103 
      + v203 * xs) + ys * (v013 + v113 * xs + v023 * ys) + z * (v004 
      + v104 * xs + v014 * ys + z * (v005 + v006 * z)))))
              )

    sigmai_dep_out = (1 / v) - 1000

    return sigmai_dep_out

def linearEOS(T_vec, S_vec, rho0=1027.0, T0=10.0, S0=35.0,
                           alp0=1.6550e-1, bet0=7.6554e-1):
    
    """
    linear equation of state, numbers taken from the NEMO model, but note that NEMO uses
          
           d_a = rho / rho0 - 1

    whereas I am going to plot

           rho = rho0 * (1 + d_a)

    so all my coefficients here are multiplied by a factor of rho0 and with the 
    more conventional units of e.g. [alpha] = K-1 instead of [alpha_before] = K-1 kg m-3
    """
    alp0 /= rho0
    bet0 /= rho0
    
    dens  = np.zeros((len(T_vec), len(S_vec)))
    for j in range(len(T_vec)):
        dens[j, :] = rho0 * (1 - alp0 * (T_vec[j] - T0) + bet0 * (S_vec[:] - S0))
        
    # output parameters in a dictionary
    params = {"alp0" : alp0, "bet0" : bet0}
        
    return params, dens

def toy_nonlinearEOS(T_vec, S_vec, rho0=1027.0, T0=10.0, S0=35.0,
                                   alp0=1.6550e-1, bet0=7.6554e-1,
                                   lam1 = 5.9520e-2,       # cabbeling coeff T^2
                                   lam2 = 5.4914e-4,       # cabbeling coeff S^2
                                   nu   = 2.4341e-3,       # cabbeling coeff T S
                                   mu1  = 1.4970e-4,       # thermobaric in T
                                   mu2  = 1.1090e-5,       # thermobaric in S
                                   z    = 0):
    """
    toy nonlinear equation of state from Vallis (2006), numbers taken from the NEMO model, 
    but note that NEMO uses
    
           d_a = rho / rho0 - 1

    whereas I am going to plot

           rho = rho0 * (1 + d_a)

    so some of my coefficients here are multiplied by a factor of rho0 and with the 
    more conventional units of e.g. [alpha] = K-1 instead of [alpha_before] = K-1 kg m-3
    """
    alp0 /= rho0
    bet0 /= rho0
    nu   /= rho0
    
    dens  = np.zeros((len(T_vec), len(S_vec)))
    for j in range(len(T_vec)):
        Ta = T_vec[j] - T0
        Sa = S_vec[:] - S0
        dens[j, :] =rho0 * (1 - alp0 * (1.0 + 0.5 * lam1 * Ta + mu1 * z) * Ta 
                              + bet0 * (1.0 - 0.5 * lam2 * Sa - mu2 * z) * Sa
                              - nu * Ta * Sa)
        
    # output parameters in a dictionary
    params = {"alp0" : alp0, "bet0" : bet0, 
              "lam1" : lam1, "lam2" : lam2,
              "nu"   : nu,
              "mu1"  : mu1,  "mu2"  : mu2,
              "z"    : 0}
    return params, dens
