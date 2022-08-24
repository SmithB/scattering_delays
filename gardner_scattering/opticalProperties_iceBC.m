function [w, g, t, Kext, Kabs, Ksca, wvl] = opticalProperties_iceBC(rb, N, d, c)
%% THIS FUNCTION IS BASTERDIZED FROM WORK DONE IN:
%
% Gardner, A. S., and M. J. Sharp (2010), A review of snow and ice albedo
% and the development of a new physically based broadband albedo
% parameterization, J. Geophys. Res., 115(F1), F01009,
% doi:10.1029/2009jf001444 

% last updated: 18-08-06
%   - function modified to be more generic

% last updated: 26-02-09
%   -  modified to include black carbon in ice

%% Determines optical depths, single scattering albedo and asymmetry
% parameter for a specified planer layer of bubbly ice for wavelengths
% WL = 0.20:0.01:4

% Inputs:
%   d = layer depth [m];
%   rb = effective radi of air bubbles in ice [mm]
%   N = bubble concentraton [# mm-3];
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
if length(d) ~= length(N)
    error('# of depths (z) must equal numer of bubble concentrations (N)')
end

if length(d) ~= length(rb)
    error('# of depths (z) must equal numer of effective bubble radi (rb)')
end

if nargin > 3
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
% wavelengths (DO NOT CHANGE FROM WL = 0.20:0.01:4)
wvl = 0.20:0.01:4; % [um]

% set ds = 1000 [kg m-3] for snow water eqivalent
ds = 1E6;                               % snow density [g m-3]
di = 0.91E6;                            % density of pure ice [g m-2]

m = length(d);
n = length(wvl);
qsca = zeros(m,n);
Ksca = qsca;
g = qsca;
t = qsca;

%% absorption coefficient
% look up the wavelength dependent refractive index of ice
[Re, Im] = refracICE(wvl);

%% determine the effective bubble concentrations
% volume fraction of air in ice sample
VFa = (4/3 * pi * rb .^ 3 .* N );

%% check to see if volume of air exceeds unit volume
if sum(VFa > 1) > 0
    error ('effective bubble concentration (VFa) > unity')
end

%% absorption coefficient for w.e. depth
% Number of air bubbles per unit volume [m w.e]
NB = ds ./(di .* (1-VFa)) .* (N * 1E9);
caB = pi * NB .*  (rb*1E-3).^2;         % bubble cross section area [m2/m3]

% The reciprocal of the absorption coefficient is the average distance a
% photon travels in pure bubble-free ice before being absorbed.
Kabs = (4 * pi * Im') ./ (wvl * 1E-6);
Kabs = ((ds ./ di) .* ones(m,n)) .* (ones(m,1) * Kabs);

%% Look up Mie scattering and determine scattering coefficinets
for i = 1:m
    % look-up bubble scattering efficiencies and asymetery factors
    % for all wavelength (0.2-4 um)
    [X, qsca(i,:), g(i,:)] = MieIceBubble_LookUp(rb(i));
    % Qsca = 2.0;  % set constant (this approx. introduces minimal error)
    
    % The scattering coefficient ksca is the reciprocal of the average
    % distance a photon travels through the bubbly ice before it is
    % scattered [mm-4]
    % Ksca = qsca * cross sectional area per bubble * # of bubble
    Ksca(i,:) = qsca(i,:) * caB(i);
end


%% Account for any black carbon (BC) in ice
if ~isnan(c)
    % Carbon partical size, depends on particle generation, ageing, and
    % contact with moisture
    %   -> typically grain size ~0.1um (Henson, 2004)
    rBC = 0.1 * 1E-3;                        % particle radius [mm]
    [qextBC qscaBC gBC] = MieBC_LookUp(rBC);
    
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
   
    %% determine BC weighting coefficinets [cross-sectional weighting]
    % set ds = 1000 [kg m-3] for snow water eqivalent
    vBC = 4/3 * pi * (rBC*1E-3)^3;          % BC sphere volume [m3]
    dBC = 1.8 * 1E6;                        % BC density [g/m3](Bond, 2008)
    mBC = 1E-6 * c * ds;                    % mass of BC [g/(m3 of snow)]
    NBC = mBC/(vBC * dBC);                  % number of BC particles [#/m3]
    caBC = NBC * pi * (rBC*1E-3)^2;         % BC cross section area [m2/m3]
    caI = ds./di;                           % ice cross section area [m2/m3]
    fB = caB ./ (caB + caBC);               % bubble weighting factor
    fBC = 1 - fB;                           % BC weighting factor

    % convert BC single scattering properties to coefficients
    % Ksca = qsca * cross sectional area per bubble * # of particles
    KscaBC = (caBC * ones(1,n)) .* (ones(m,1) * qscaBC');
    KabsBC = (caBC * ones(1,n)) .* (ones(m,1) * qabsBC');
    gBC = gBC';
    
    % modify extinction and scattering efficiencies to account for BC
    for i = 1:m
        Kabs(i,:) = Kabs(i,:) + KabsBC(i,:);
        Ksca(i,:) = Ksca(i,:) + KscaBC(i,:);
        g(i,:) = fB(i) * g(i,:) + fBC(i) * gBC;

    end
end

%% Extinction coef., single scattering and optical thickness

Kext = Ksca + Kabs;                 % extinction coefficient
w = Ksca ./ Kext;                   % single scattering albedo
for i = 1:m; 
    t(i,:) = d(i) .* Kext(i,:);     % optical depth
end

end
