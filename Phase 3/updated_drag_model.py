import numpy as np
import matplotlib.pyplot as plt
from subfunctions_EDL import *
import scipy.interpolate as sp
import subfunctions_EDL as subEDL
import functools as ft

# Drag Coefficient Modifier Function

def __get_coefficients_for_drag_model():
    mach_exp = np.array([0.25,0.5,0.65,0.7,0.8,0.9,0.95,1,1.1,1.2,1.3,1.4,1.5,1.6,1.8,1.9,2,2.2,2.5,2.6])
    MEF_exp = np.array([1,1,1,0.98,0.9,0.72,0.66,0.76,0.9,0.96,0.99,0.999,0.992,0.98,0.9,0.85,0.82,0.75,0.65,0.62])

    nsegs = np.size(mach_exp, 0) - 1

    # convert mach_exp and MEF_exp to column vectors
    mach_exp = mach_exp.reshape((len(mach_exp),1))
    MEF_exp = MEF_exp.reshape((len(MEF_exp),1))

    # initialize matrix A and vector b
    a = np.zeros((2*nsegs, 2*nsegs))
    b = np.zeros((2*nsegs, 1))

    # fill in matrix A and vector b
    for i in range(nsegs):
        a[2*i, 2*i] = mach_exp[i]
        a[2*i, 2*i+1] = 1
        a[2*i+1, 2*i] = mach_exp[i+1]
        a[2*i+1, 2*i+1] = 1
        b[2*i] = MEF_exp[i]
        b[2*i+1] = MEF_exp[i+1]

    # solve for the coefficients of the piecewise linear function
    c = np.linalg.solve(a, b)

    return c

def __full_get_MEF(v, alt, c, plot = False):
    """
    This function calculates the Mach efficiency factor for a given velocity and altitude using interpolated data from the documentation.
    It returns the Mach efficiency factor.
    
    Inputs:
        v: velocity in m/s
        alt: altitude in m
        plot: boolean, if true, plots the drag coefficient modifier
    
    Outputs:
        MEF: Mach efficiency factor
    
    """
    
    mach_exp = np.array([0.25,0.5,0.65,0.7,0.8,0.9,0.95,1,1.1,1.2,1.3,1.4,1.5,1.6,1.8,1.9,2,2.2,2.5,2.6])
    MEF_exp = np.array([1,1,1,0.98,0.9,0.72,0.66,0.76,0.9,0.96,0.99,0.999,0.992,0.98,0.9,0.85,0.82,0.75,0.65,0.62])

    nsegs = np.size(mach_exp, 0) - 1

    # convert mach_exp and MEF_exp to column vectors
    mach_exp = mach_exp.reshape((len(mach_exp),1))
    MEF_exp = MEF_exp.reshape((len(MEF_exp),1))

    # solve for the MEF at the given mach number
    MEF = 0
    mach = subEDL.v2M_Mars(v, alt)
    i = np.where(mach_exp > mach)[0][0] - 1
    MEF = c[2*i]*mach + c[2*i+1]


    # calculate the MEF at the given mach number from the coefficients more efficiently
    # find the first index where mach_exp is greater than mach




    

    if plot:
        # plot the piecewise linear function
        mach_test = np.linspace(0.25, 2.6, 10000)
        mach_exp = np.array([0.25,0.5,0.65,0.7,0.8,0.9,0.95,1,1.1,1.2,1.3,1.4,1.5,1.6,1.8,1.9,2,2.2,2.5,2.6])
        MEF_exp = np.array([1,1,1,0.98,0.9,0.72,0.66,0.76,0.9,0.96,0.99,0.999,0.992,0.98,0.9,0.85,0.82,0.75,0.65,0.62])

        MEF_test = np.zeros((len(mach_test), 1))
        for i in range(nsegs):
            for j in range(len(mach_test)):
                if mach_test[j] >= mach_exp[i] and mach_test[j] <= mach_exp[i+1]:
                    MEF_test[j] = c[2*i]*mach_test[j] + c[2*i+1]
        plt.figure()
        plt.plot(mach_test, MEF_test, 'mediumpurple', linewidth=2, label='Linear Splines')
        plt.plot(mach_exp, MEF_exp, 'ko', label='Experimental Data')
        plt.plot(mach, MEF, 'blue', label='Interpolated Data', marker='o', markersize=10)
        plt.title('Drag Coefficient Modifier Function')
        plt.xlabel('Mach Number')
        plt.ylabel('Mach Efficiency Factor')
        plt.legend()
        plt.grid()
        plt.savefig('DCM.png')
        plt.tight_layout()
        plt.show()

    # # Testing different interpolation methods, not used in final code

    # mach_fit_quad = sp.interp1d(mach_exp, MEF_exp, kind='quadratic')
    # mach_fit_cubic = sp.interp1d(mach_exp, MEF_exp, kind='cubic')
    # mach_fit_slinear = sp.interp1d(mach_exp, MEF_exp, kind='slinear')
    # mach_fit_linear = sp.interp1d(mach_exp, MEF_exp, kind='linear')

    # plt.figure()
    # plt.plot(mach_test, mach_fit_quad(mach_test), 'grey', linewidth=2, label='Quadratic')
    # plt.plot(mach_test, mach_fit_cubic(mach_test), 'r', linewidth=2, label='Cubic')
    # # plt.plot(mach_test, mach_fit_slinear(mach_test), 'b', linewidth=2, label='Slinear')
    # # plt.plot(mach_test, mach_fit_linear(mach_test), 'g', linewidth=2, label='Linear')
    # plt.plot(mach_exp, MEF_exp, 'ko', label='Experimental Data')
    # plt.title('Drag Coefficient Modifier Function')
    # plt.xlabel('Mach Number')
    # plt.ylabel('Mach Efficiency Factor')
    # plt.legend()
    # plt.grid()
    # plt.savefig('CubicQuadTest.png')
    # plt.show()

    # # Error Analysis

    # quad = mach_fit_quad(mach_test)
    # cubic = mach_fit_cubic(mach_test)

    # # relative error between cubic and linear
    # rel_err_cubic = np.zeros((len(mach_test), 1))
    # rel_err_quad = np.zeros((len(mach_test), 1))
    # for i in range(len(mach_test)):
    #     rel_err_cubic[i] = (abs(cubic[i] - MEF_test[i]) / MEF_test[i]) * 100
    #     rel_err_quad[i] = (abs(quad[i] - MEF_test[i]) / MEF_test[i]) * 100

    # # plot relative error between cubic/linear and quad/linear in 2x1 subplot
    # plt.figure()
    # plt.subplot(2, 1, 1)
    # plt.plot(mach_test, rel_err_cubic, 'r', linewidth=2, label='Cubic')
    # plt.title('Relative Error Between Cubic/Linear and Quad/Linear')
    # plt.xlabel('Mach Number')
    # plt.ylabel('Relative Error (%)')
    # plt.legend()
    # plt.grid()
    # plt.subplot(2, 1, 2)
    # plt.plot(mach_test, rel_err_quad, 'g', linewidth=2, label='Quadratic')
    # plt.xlabel('Mach Number')
    # plt.ylabel('Relative Error (%)')
    # plt.legend()
    # plt.grid()
    # plt.savefig('RelError.png')
    # plt.show()

    return float(MEF)

c = __get_coefficients_for_drag_model()
__get_MEF = ft.partial(__full_get_MEF, c=c)

def get_MEF(v, alt, plot=False):
    return __get_MEF(v, alt, plot=plot)

if __name__ == '__main__':

    print('starting')

    MEF = get_MEF(400, 0)
    print(MEF)

