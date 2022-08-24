function [qext0, qsca0, g0] = MieIceBubble_LookUp(r0, WL0)
%% Mie qext qsca g look-up function for bubbles in ice

% Look up parameters
%   r = effective grain radius [mm]
%   WL = wavelength [um] -> OPTIONAL

% Return parameters
%   qext = Mie extinction efficiencies of size [m n]
%   qsca = Mie scatter efficiencies of size [m n]
%   g = asymmetery factor of size [m n]

% NOTE:
%   If WL is not specified then function returns values for wavelengths
%   WL = 0.2:0.01:4
%
%   function linearly interpolates between effective grain sizes 
%   and wavelengths (if specified)

%% CHECK NUMBER OF FUNCTION INPUTS
if nargin == 2, allWL = 0; else allWL = 1; end

%% LOAD IN MIE VALUES
load MieIceBubble_LookUp

%% CHECK IF LOO-UP VALUES ARE OUT OF TABLE RANGE
if r0 > max(r) || r0 < min(r)
    error('effective grain radius r is out of MieIceBubble_LookUp range')
elseif nargin == 2 
    if WL0 > max(WL) || WL0 < min(WL)
     error('wavelength WL is out of MieIceBubble_LookUp range')
    end
end

%% LOOK-UP MIE VALUES
n = find(r0==r); % check if r0 is an exact look-up grain size
if isempty(n), D = 0; else D = 1; end
switch allWL
    case 0 % extract values for a single wavelength
        m = find(roundn(WL0, -3)==roundn(WL, -3)); % check if WL0 is an exact look-up wavelength
        if isempty(m), W = 0; else W = 1; end
        
        % r0 & WL0 are exact look-up values
        if D == 1 && W == 1 
            qext0 = qext(m,n);
            qsca0 = qsca(m,n);
            g0 = g(m,n);
            
        % WL0 is an exact look-up value
        elseif D == 0 && W == 1 
                n = find(r0<r, 1, 'first');

                % linear weightings
                j = (r0 - r(n-1))/(r(n) - r(n-1));
                k = (r(n)-r0)/(r(n) - r(n-1));

                % linear weighted values
                qext0 = qext(m,n-1)*k + qext(m,n)*j;
                qsca0 = qsca(m,n-1)*k + qsca(m,n)*j;
                g0 = g(m,n-1)*k + g(m,n)*j;
               
        % r0 is an exact look-up value   
        elseif D == 1 && W == 0 
                m = find(WL0<WL, 1, 'first');

                % linear weightings
                h = (WL0 - WL(m-1))/(WL(m) - WL(m-1));
                i = (WL(m)-WL0)/(WL(m) - WL(m-1));

                % linear weighted values
                qext0 = qext(m-1,n)*i + qext(m,n)*h;
                qsca0 = qsca(m-1,n)*i + qsca(m,n)*h;
                g0 = g(m-1,n)*i + g(m,n)*h;
                
        % r0 & WL0 are both not exact look-up values    
        else
                m = find(WL0<WL, 1, 'first');
                n = find(r0<r, 1, 'first');

                % linear weightings
                h = (WL0 - WL(m-1))/(WL(m) - WL(m-1));
                i = (WL(m)-WL0)/(WL(m) - WL(m-1));
                j = (r0 - r(n-1))/(r(n) - r(n-1));
                k = (r(n)-r0)/(r(n) - r(n-1));
                
                % linear weighted values
                % weighted for wavelength
                qext1 = qext(m-1,n-1)*i + qext(m,n-1)*h;
                qext2 = qext(m-1,n)*i + qext(m,n)*h;
                qsca1 = qsca(m-1,n-1)*i + qsca(m,n-1)*h;
                qsca2 = qsca(m-1,n)*i + qsca(m,n)*h;
                g1 = g(m-1,n-1)*i + g(m,n-1)*h;
                g2 = g(m-1,n)*i + g(m,n)*h;

                % weighted for effective grain size
                qext0 = qext1*k + qext2*j;
                qsca0 = qsca1*k + qsca2*j;
                g0 = g1*k + g2*j;

                % extract values for all wavelengths
        end

        case 1 % extract values for all wavelengths
            % r0 is an exact look-up values
            if D == 1
                qext0 = qext(:,n);
                qsca0 = qsca(:,n);
                g0 = g(:,n);
                
            % r0 is not an exact look-up value
            else
                n = find(r0<r, 1, 'first');

                % linear weightings
                j = (r0 - r(n-1))/(r(n) - r(n-1));
                k = (r(n)-r0)/(r(n) - r(n-1));

                % linear weighted values
                qext0 = qext(:,n-1)*k + qext(:,n)*j;
                qsca0 = qsca(:,n-1)*k + qsca(:,n)*j;
                g0 = g(:,n-1)*k + g(:,n)*j;
            end
end
end


