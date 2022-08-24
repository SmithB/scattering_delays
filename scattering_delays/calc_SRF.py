#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  3 17:45:25 2022

@author: ben
"""

import numpy as np
from . import read_MC_results
from . import  v_eff
from . import optical_properties_pure_snow
from . import R_from_MC
from . import read_MC_results


def gaussian(t, μ, σ):
    return 1/(σ*np.sqrt(2*np.pi))*np.exp(-1/2*((t-μ)/σ)**2)

def calc_srf(t_wf, mu_s, mu_a, g, v, D_MC, sigma_wf=0.1*1.e-9, τ_exp=np.Inf,
             skip_tail=False):
    '''
    Calculate a single set of surface-response function values based on a set
    of bulk scattering properties.

    Parameters
    ----------
    t_wf : np.array
        time values for which to genereate the waveform The default is None.
    μ_s : float
        Scattering cross-section.
    μ_a : float
        Absorption cross-section.
    g : float
        Asymmetry factor.
    v : float
        velocity in medium
    D_MC : np.array
        Monte-carlo output data.
    sigma_wf : float, optional
        standard deviation to use when averaging monte-carlo data. The default is 0.1*1.e-9.
    τ_range_tail : iterable, optional
        range of optical thicknesses over which a t^-3/2 curve should be fit to the 
        waveform (see R_from_MC.py). Defautls to [200, np.Inf],

    Returns
    -------
    numpy array
        Power values for time values.
    float
        Power returned beyond the end of t_WF

    '''

    # choose a time vector based on the scattering coefficients
    L_s=1/mu_s;
    t_s=L_s/v;
    t_win = 2*t_wf[-1]
    
    t1=np.logspace(np.log10(0.025*t_s), np.log10(t_win), 3000)
    t1=np.concatenate((t1, t1[-1]+np.arange(0.1, 40.1, 0.1)*1e-9))

    # calculate the energy in the return,
    R, R_tail = R_from_MC( t1, mu_s, mu_a, g, v, D_MC, τ_exp=τ_exp )
    Rbar=(R[1:]+R[0:-1])/2
    dt=np.diff(t1)
    E_tot=np.sum(Rbar*dt)+R_tail

        # Algebraic convolution between the IRF and the SRF
    delta_t=np.concatenate((t1[1:]-t1[0:-1], [t1[-1]-t1[-2]]))
    G=np.zeros((len(t_wf), len(t1)));
    for k, ti in enumerate(t_wf):
        these=np.abs(t1-ti) < 12*sigma_wf;
        if np.any(these):
            # the kth row of G multiplies the attenuated normalized photon counts
            # by a Gaussian centered on t_wf(k) evaluated at t1, times the delta_t values associated with t1 (dt)
            G[k,these]=gaussian(t1[these], ti, sigma_wf)*delta_t[these];

    return G.dot(R), E_tot


def make_srf_catalog( MC_dir=None, N_MC=None, D_MC=None, t_wf=None, r_vals=None, μ_s=None, μ_a=0, g=0, v=None, ρ=None, λ=532e-9, τ=None ):
    '''
    Calculate the surface response function for a set of optical properties

    Parameters
    ----------
    MC_dir : str, optional
        Directory containing MC results.
    N_MC : int, optional
        Number of photons launched in each MC file.
    D_MC: np.array, optional
        Data read from an MC file
    t_wf : np.array, optional
        time values for which to genereate the waveform The default is None.
    r_vals : TYPE, optional
        effective grain-size radii.  Used if scattering properties are not specified. The default is None.
    μ_s : iterable, optional
        Scattering cross-section values. The default is None.
    μ_a : iterable, optional
        Absorption cross-section values. The default is None.
    g : iterable, optional
        Asymmetry factors. The default is None.
    v : float, optional
        velocity in medium. The default is None.
    ρ : float, optional
        density, The default is None.
    λ : float, optional
        wavelength. The default is 532e-9.
    z : float, optional
        Optical thickness of the layer. The default is None.

    Returns
    -------
    SRF: dictionary whose keys specify either scaattering coefficients or gain sizes,
    and whose entries specify times, normalized power values, and post-tail power

    '''
    if D_MC is None:
        D_MC=read_MC_results(MC_dir, N_MC, τ=τ)
    if t_wf is None:
        dt=0.25e-9
        t_wf=np.arange(-2.e-8, 2.e-8+dt, dt)


    if v is None:
        if ρ is None:
            v = 3.e8
        else:
            v = v_eff(λ, ρ)

    # If we have not specified a scattering length, calculate it from r
    if μ_s is None:
        LUT={}
        if r_vals is None:
            r_vals=np.logspace(np.log10(2.5e-5), np.log10(4.e-2), 71)
        K0=r_vals

        μ_s, μ_a, g = [[], [], []]
        for r in r_vals:
            mu_e, mu_s, mu_a, gi = optical_properties_pure_snow(r, λ, ρ, LUT=LUT)
            μ_s+=[mu_s]
            μ_a+=[mu_s]
            g += [gi]
    else:
        K0=μ_s

    SRF={}
    for k_i, mu_s, mu_a, gi in zip(K0, μ_s, μ_a, g):
        SRF[k_i]={}
        # calculate the optical properties
        P, E_tail=calc_srf(t_wf, mu_s, mu_a, gi, v, D_MC)
        # save the output
        SRF[k_i]['t']=t_wf
        SRF[k_i]['P']=P
        SRF[k_i]['E_tot']=E_tail

    return SRF

