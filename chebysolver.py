import numpy as np
import matplotlib.pyplot as plt



#Function to calculate chebyshev polynomials


def chebshev(N, x):

    order = N

    length = len(x)

    T_matrix = np.zeros((N, length))

    for i in range(order):

        if i == 0:

            T_matrix[0, :] = 1

        elif i ==1:

            T_matrix[i, :] = x

        else:

            T_matrix[i, :] = 2* x * T_matrix[i-1, :] - T_matrix[i-2, :]

    return T_matrix


def chebderivatives(N, x):

    theta = np.cos(x)

    n = np.linspace(0, N, N+1)

    #Calculate derivatives of cos(n*theta)

    T_cos = np.zeros((len(n), len(theta)))
    T_cos_t = np.zeros((len(n), len(theta)))
    T_cos_tt = np.zeros((len(n), len(theta)))
    T_cos_ttt = np.zeros((len(n), len(theta)))
    T_cos_tttt = np.zeros((len(n), len(theta)))

    T_x = np.zeros((len(n), len(theta)))
    T_xx = np.zeros((len(n), len(theta)))
    T_xxx = np.zeros((len(n), len(theta)))
    T_xxxx = np.zeros((len(n), len(theta)))

    for i in range(len(n)):

        T_cos[i, :] = np.cos(n[i]*theta)
        T_cos_t[i, :] = -n[i]*np.sin(n[i]*theta)
        T_cos_tt[i, :] = -n[i]**2*np.cos(n[i]*theta)
        T_cos_ttt[i, :] = n[i]**3*np.sin(n[i]*theta)
        T_cos_tttt[i, :] = n[i]**4*np.cos(n[i]*theta)

    sin = np.sin(theta)
    cos = np.cos(theta)

    for i in range(len(n)):

        T_x[i, :] = -T_cos_t[i, :]/sin
        T_xx[i, :] = (sin*T_cos_tt[i, :]-cos*T_cos_t[i, :])/(sin**3)
        T_xxx[i, :] = (-sin*sin*T_cos_ttt[i, :]+3*cos*sin*T_cos_tt[i, :]-(3*cos*cos+sin*sin)*T_cos_t[i, :])/(sin**5)
        T_xxxx[i, :] = (sin*sin*sin*T_cos_tttt[i, :]-6*cos*sin*sin*T_cos_ttt[i, :]+(15*cos*cos*sin+4*sin*sin*sin)*T_cos_tt[i, :]-(9*cos*sin*sin+15*cos*cos*cos)*T_cos_t[i, :])/(sin**7)


    return T_x, T_xx, T_xxx, T_xxxx

#define collaction points

N = 4
x = np.linspace(-1, 1, 10)
Tn = chebshev(N+1, x)
Tx, Txx, Txxx, Txxxx = chebderivatives(N, x)

G = np.zeros((len(x), N))
f = np.zeros(len(x))

#boundary conditions
f[0] = alpha0
f[1] = alpha1
f[2] = beta0
f[3] = beta1

#u(-1)
G[0, :] = Tn[:, 0]
G[1, :] = Tx[:, 0]

#u(1)
G[2, :] = Tn[:, len(x)-1]
G[3, :] = Tx[:, len(x)-1]
