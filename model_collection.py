# JM (19 Jan 2022): collection of models

import numpy as np

def mod_nicolson_bailey(k, r, n=200, a=1.0):
    """
    (Modified) Nicolson-Bailey model, with non-dimensional inputs

    k    the carrying capacity
    r    some positive constant
    n    number of generations to run (optional)
    a    measure of some sort of death/survival rate (optional)
    
    Outputs the arrays "host" and "para" as a time-series in generation

    """

    host_init, para_init = k / (k-1) * np.log(k) / a, np.log(k) / a + 1e-1

    host, para = np.zeros(n), np.zeros(n)
    host[0], para[0] = host_init, para_init

    for i in range(n-1):
        host[i+1] = host[i] * np.exp(r * (1.0 - host[i] / k)) * np.exp(-a * para[i])
        para[i+1] = host[i] * (1.0 - np.exp(-a * para[i])) 

    return host, para
    
def noisy_lotka_volterra(n=200, noise=1e-16, a=2/3, b=4/3, c=1.0, d=1.0, prey_ini=2.0, pred_ini=2.0):
    """
    (Stochastic) Lotka-Volterra equation, with (optional) non-dimensional inputs

    n         number of generations to run
    noise     stochastic noise level (default is basically off)
    a         growth rate of prey
    b         death rate of prey by predation
    c         death rate of predator
    d         growth rate of predator by predation
    prey_ini  initial value for prey
    pred_ini  initial value for predator
    
    Outputs the arrays "prey" and "predator" as a time-series in non-dimensional time

    """

    prey, predator = np.zeros(n), np.zeros(n)

    prey[0], predator[0] = prey_ini, pred_ini
    dt = 0.2

    i = 0 # Euler step
    prey[i+1]     = prey[i]     + dt * (a * prey[i]               - b * prey[i] * predator[i])
    predator[i+1] = predator[i] + dt * (d * prey[i] * predator[i] - c           * predator[i])

    # Adams-Bashforth 2 scheme
    for i in range(0, n-2):
        prey[i+2]     = prey[i+1]     + dt * ( 3/2 * (a * prey[i+1]                 - b * prey[i+1] * predator[i+1])
                                              -1/2 * (a * prey[i  ]                 - b * prey[i  ] * predator[i  ])
                                             ) + noise * np.random.rand()
        predator[i+2] = predator[i+1] + dt * ( 3/2 * (d * prey[i+1] * predator[i+1] - c             * predator[i+1])
                                              -1/2 * (d * prey[i  ] * predator[i  ] - c             * predator[i  ])
                                             )
        
    return prey, predator
