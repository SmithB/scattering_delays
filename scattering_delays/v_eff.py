#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Wed July 6 11:30:02 2022

@author: ben
"""

'''
    Calculate the speed of light in snow based on a linear mixing law

    Parameters
    ----------
    λ: numeric
        wavelength, in m
    ρ: numeric
        density, in kg m^-3

    Returns:
    --------
    numeric: speed of light, in m/s

'''

from .ice_n_of_lambda import ice_n_of_lambda
import numpy as np

def v_eff(λ, ρ):

    ρ_ice=917;
    f_ice=ρ/ρ_ice;
    n_ice=ice_n_of_lambda(λ);
    n_air=1;

    n_snow=f_ice*n_ice + (1-f_ice)*n_air;

    return 3.e8/np.real(n_snow)
