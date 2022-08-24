function [w, g, t, Kext, Kabs, Ksca, wvl] = opticalProperties_snowBC(r, d, c)
%% THIS FUNCTION IS BASTERDIZED FROM WORK DONE IN:
%
% Gardner, A. S., and M. J. Sharp (2010), A review of snow and ice albedo
% and the development of a new physically based broadband albedo
% parameterization, J. Geophys. Res., 115(F1), F01009,
% doi:10.1029/2009jf001444 

% last updated: 18-08-06
%   - function modified to be more generic

% last updated: 27-02-09
%   - function returns same results as before but has been modified to
%     match the methods used in inDISORT_iceBC

%% Determines optical depths, single scattering albedo and asymmetry
% parameter for a specified planer layer of snow for wavelengths
% wvl = 0.20:0.01:4 [um]

% Inputs:
%   d = layer depth [m];
%   r = effective radi of snow grains[mm]
%   c = black carbon concentration [ppmw] [OPTIONAL]

% Outputs:
%   t = optical depth
%   w = extinction/scattering efficiency (single scattering albedo)
%   g = asymmetry parameter (function of wavelength)
%   Kext = extinction coefficients
%   Kabs = absorption coefficients
%   Ksca = scattering coefficients
%   wvl = wavelengths [um]
%% Checks
if length(d) ~= length(r)
    error('# of depths (z) must equal number of effective radii (r)')
end

if nargin > 2
    if length(d)~= length(c)
        error('# of depths (z) must equal number of BC concentrations (c)')
    end
    
    if sum(c ~= 0) == 0;
        c = NaN;
    end
    
else
    % set black carbon concentration to 'NaN'
    c = NaN;
end

%% INITALIZE
% wavelengths (DO NOT CHANGE FROM wvl = 0.20:0.01:4)
wvl = 0.20:0.01:4;
ds = 1E6;                               % snow density [g m-3]
di = 0.91E6;                            % density of pure ice [g m-2]

m = length(d);
n = length(wvl);
qext = zeros(m,n);
qsca = qext;
g = qext;
t = qext;


%% Account for any black carbon (BC) in snow
% Carbon partical size, depends on particle generation, ageing, and
% contact with moisture
%   -> typically grain size ~0.1um (Henson, 2004)
rBC = 0.1 * 1E-3;                        % particle radius [mm]

[qextBC, qscaBC, gBC] = MieBC_LookUp(rBC);

%% account for internal mixing
% Empirical measurements suggest that absorption by light-absorbing
% carbon are generally larger than predicted by models (Warren and
% Wiscombe 1980) . We take a middle of the road approach to the two
% extremes suggested by Hansen and Nazarenko (2004) and increase the
% effective absorbance of light-absorbing carbon by multiplying the
% absorption efficiency by 1.5.
% qabsBC = (qextBC - qscaBC) * 1.5;

% do calculations without adjusted absorption efficiency               changed 16-03-09
qabsBC = (qextBC - qscaBC);
qextBC = qabsBC + qscaBC;

%% determine BC weighting coefficinets [cross-sectional weighting]
% set ds = 1000 [kg m-3] for snow water eqivalent
vBC = 4/3 * pi * (rBC*1E-3)^3;          % BC sphere volume [m3]
dBC = 1.8 * 1E6;                        % BC density [g/m3](Bond, 2008)
mBC = 1E-6 * c * ds;                    % mass of BC [g/(m3 of snow)]
NS = (3 * ds) ./ ...                    % number of snow grains [#/m3]
    (4 * di * pi * ((r*1E-3).^3));
NBC = mBC/(vBC * dBC);                  % number of BC particles [#/m3]
caBC = NBC * pi * (rBC*1E-3).^2;        % BC cross section area [m2/m3]
caS =  pi * NS .* (r*1E-3).^2;          % snow cross sect area [m2/m3]
fS = caS ./ (caS + caBC);               % snow weighting factor
fBC = 1 - fS;                           % BC weighting factor

% initialize optical coefficients
Kext =  zeros(size(qext));
Ksca =  Kext;

for i = 1:m
    % Kext = qext * cross sectional area per bubble * # of particles
    [qextS, qscaS, gS] = MieIce_LookUp(r(i));
    Kext(i,:) =  qextS .* caS(i) + qextBC .* caBC(i);
    Ksca(i,:) =  qscaS .* caS(i) + qscaBC .* caBC(i);
    g(i,:) = fS(i) * gS + fBC(i) * gBC;
    t(i,:) = Kext(i,:) .* d(i);      % optical depth
    
    % Above is the same as ...
    % qext(i,:) = fS(i) * qextS + fBC(i) * qextBC;
    % qsca(i,:) = fS(i) * qscaS + fBC(i) * qscaBC;
    % g(i,:) = fS(i) * gS + fBC(i) * gBC;
    % t(i,:) = (3*d(i)*qext(i,:)*ds)./(4*(r(i)/1E3)*di);
end

w = Ksca./Kext;     % single scattering albedo
Kabs = Kext - Ksca;