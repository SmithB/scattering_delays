#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed July 6 11:15:22 2022

@author: ben
"""


import numpy as np
from pathlib import Path
import os
import h5py
import scipy.interpolate as si

def optical_properties_pure_snow(r, λ, ρ_s, LUT=None, skip_calc=False):
    '''
    Calculate scattering properties for pure snow


    Parameters
    ----------
    r : numeric
        grain radius, in m.
    λ : numeric
        wavelength, in m.
    ρ_s : numeric
        snow density, in kg m^-3.
    LUT : dict, optional
        Dictionary in which the look-up tables are stored.  If provided as an
        empty dict, it can be reused in subsequent calculations. The default is None.
    skip_calc : boolean, optional
        If true, the look-up table is generated, but no properties are calculated. The default is False.

    Returns
    -------
    Kext : numeric
        Extinction coefficient, m^-1.
    Ksca : numeric
        Scattering coefficient, m^-1.
    Kabs : numeric
        Absorption coefficient, m^-1.
    g : numeric
        Asymmetry parameter.

    Note: this script was adapted from optical_properties_snow_bc by Alex
    Gardner (JPL)

    '''


    ρ_i=910 # kg m^-3, this is what Alex used
    if LUT is None or 'WL' not in LUT:
        path = str(Path(__file__).parent.parent.absolute())
        thefile=os.path.join(path, 'gardner_scattering', 'MieIce_LookUp.mat')
        if not isinstance(LUT, dict):
            LUT={}
        with h5py.File(thefile,'r') as h5f:
            r0=np.array(h5f['r'])
            WL=np.array(h5f['WL'])
            #print({'r':[np.min(r0), np.max(r0)],
            #       'WL':[np.min(WL), np.max(WL)]})
            for field in ['qext','qsca','g']:
                temp=np.array(h5f[field])
                #print([r0.shape, WL.shape, temp.shape])
                LUT[field]=si.RectBivariateSpline(r0, WL, temp, kx=1, ky=1)
    if skip_calc:
        return
    qext = LUT['qext'].ev(r*1000, λ*1.e6)
    qsca = LUT['qsca'].ev(r*1000, λ*1.e6)
    g = LUT['g'].ev(r*1.e3, λ*1.e6)

    # number of snow grains [#/m3]
    NS = (3 * ρ_s) / (4 * ρ_i * np.pi * (r**3))
    # snow cross sect area [m2/m3]
    caS =  np.pi * NS * (r**2);
    Kext =  qext * caS
    Ksca =  qsca * caS
    Kabs  = Kext-Ksca
    return Kext, Ksca, Kabs, g
