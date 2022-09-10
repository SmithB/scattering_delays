import numpy as np
import scipy.special as spspec


def fit_exponential_tail(tt, R, N, L_scale, τ_exp):

    τ_out=tt/L_scale

    # power where the tail diverges from the parent curve
    P0_exp = R[np.argmin(τ_out < τ_exp )]

    exp_fit_ind = np.flatnonzero( (τ_out > τ_exp ) & (R<0.5*P0_exp) & (N>50) & np.isfinite(R))
    if len(exp_fit_ind) < 5:
        return R.copy(), None
    G_exp = np.c_[np.ones(len(exp_fit_ind)), tt[exp_fit_ind]/1.e-9]
    m_exp = np.linalg.inv(G_exp.T.dot(G_exp)).dot(G_exp.T.dot(np.log(R[exp_fit_ind])))
    # undo the xcaling on m[1]
    m_exp[1] /= 1.e-9

    # Fill in the low-precision elements of the tail (including the zeros)
    exp_fill_ind=np.flatnonzero(((N<200) | np.isnan(R)) & (τ_out > τ_exp ))
    R1=R.copy()
    R1[exp_fill_ind] = np.exp(m_exp[0]+m_exp[1]*tt[exp_fill_ind])
    R_tail = -np.exp(m_exp[0])/(m_exp[1])*np.exp(m_exp[1]*tt[-1])

    return R1, R_tail

def R_from_MC(t_out, mu_s, mu_a, g, c, D_MC, τ_tail=200, τ_exp=None, skip_tail=False, L_max=1900, DEBUG=False):

    """
    Created on Wed July 6 13:10:51 2022

    @author: ben
    """

    '''
        Calculate the time distribution returning from a scattering layer

        Parameters
        ----------
        tt: numeric
            time vector
        mu_s, mu_a: numeric
            scattering and absorption coefficients in m^-1
        g: numeric
            scattering asymmetry parameter
        c: numeric
            speed of light, m/s
        D_MC: dict
            monte-carlo results, from read_MC_results
        τ_tail:
            Fit the t^-3/2 power model to the tail beyond this time; defaults to 200
        τ_exp:
            Beyond this time, expect the tail fit to follow an exponential decay, defaults to np.Inf
        skip_tail:
            If true, skip fitting the tail
        L_max:
            maximum expected time in the monte carlo

        Returns:
        --------
        R: vector
            rate of photon returns for the times in tt
        R_tail:
            Integrated energy in the return after the end of tt

    '''
    #print("-----")
    if τ_exp is None:
        τ_exp = np.Inf

    L_MC=D_MC['lc']/((1-g)*mu_s);
    t_MC=L_MC/c;
    # scales l_mc to time: a very small number.
    L_scale=1/c/((1-g)*mu_s)
    #print(f'mu_s={mu_s}\nL_scale={L_scale}\nlast t ={t_out[-1]}')

    # waveform time expressed as optical thickness
    τ_out=t_out/L_scale

    # count corrected for absorption;
    E_abs = D_MC['E']*np.exp(-mu_a*L_MC);

    P_abs = E_abs/(D_MC['lb']-D_MC['la'])
    # interpolate to get the power as a function of time
    R = np.interp(t_out, t_MC, P_abs, right=np.NaN) / L_scale
    # get the data count
    #N = np.interp(t_out, t_MC, D_MC['E']/D_MC['norm'])

    R[t_out<t_MC[0]]=0

    if skip_tail:
        R[t_out>t_MC[-1]]=0
        R_tail=0
        return R, R_tail

    R1=R.copy()

    # Powerlaw tail fit
    tail = (D_MC['lb'] > τ_tail) & (D_MC['la'] < L_max) & (D_MC['lb'] < τ_exp)
    if np.sum(tail)<2:
        A_tail=None
    else:
        G_tail = np.exp(-mu_a * L_MC[tail])*L_MC[tail]**(-3/2);
        d_tail = P_abs[tail];  # power per unit (nondimensional) length
        A_tail = 1/(G_tail.dot(G_tail))*(G_tail.T.dot(d_tail))
        power_fill_els = ((τ_out > D_MC['lc'][-1]) | np.isnan(R)) & (τ_out <= τ_exp)
        L_fill=c*t_out[power_fill_els]
        R1[power_fill_els] = A_tail * np.exp(-mu_a * L_fill) * L_fill**(-3/2)/L_scale

    if τ_exp is not None and np.isfinite(τ_exp) and τ_exp < τ_out[-1]:
        N=np.interp(t_out, t_MC, D_MC['E']/D_MC['norm'], right=np.NaN)
        R1, R_tail = fit_exponential_tail(t_out, R1, N, L_scale, τ_exp)
        if R_tail is not None:
            return R1, R_tail


    # Calculate the total energy in the tail beyond the end of the requested
    # waveform based on the attenuated power law fit
    b0 = t_out[-1]*c;
    R_tail = A_tail*(1-g)*mu_s*(2*np.exp(-b0*mu_a)/np.sqrt(b0) -
                                2*np.sqrt(np.pi*mu_a)* spspec.erfc(np.sqrt(b0*mu_a)))

    if np.isfinite(τ_exp):
        # if there is a finite exponential decay, extrapolate its slope from the
        # the tangent to the T^-3/2 curve
        # the tail continues on a tangent on semilog axes
        # log(R1) = log(A_tail) + (-mu_a * c * t)  + -3/2 * log(ct) - log(L_scale)
        # d/dt (log(R1)) = -mu_a c - 3/2*(1/t)
        t_exp = τ_exp*L_scale
        L_exp = τ_exp/((1-g)*mu_s)
        K_exp = -c*mu_a - 1.5/t_exp
        A_exp = A_tail*np.exp(-mu_a * L_exp)*L_exp**(-3/2)/L_scale
        #print(f'last tau_out = {τ_out[-1]}\n, tau_exp={τ_exp}\n t_exp={t_exp}\n L_exp={L_exp}\nK_exp={K_exp}\nA_exp={A_exp}')
        if τ_exp < np.max(τ_out):
            late_els = (t_out > t_exp) & ( t_out > t_MC[-1])
            R1[late_els] = A_exp * np.exp(K_exp*(t_out[late_els]-t_exp))
            R_tail = -A_exp/K_exp*np.exp(K_exp*(t_out[-1]-t_exp))
            return R1, R_tail
        else:
            R_tail_exp = -A_exp/K_exp
        # if the tail goes exponential at some late time, subtract the power in the
        # t^-3/2 curve beyond this time, and add the power in the exponential curve
        b_exp = t_exp*c
        R_tail_corr = A_tail*(1-g)*mu_s*(2*np.exp(-b_exp*mu_a)/np.sqrt(b_exp) -
                                    2*np.sqrt(np.pi*mu_a)* spspec.erfc(np.sqrt(b_exp*mu_a)))
        #print(f"R_tail = {R_tail} for b0={b0}; b_exp={b_exp}, subtracting {R_tail_corr}")
        R_tail -= R_tail_corr
        R_tail += R_tail_exp
    return R1, R_tail

