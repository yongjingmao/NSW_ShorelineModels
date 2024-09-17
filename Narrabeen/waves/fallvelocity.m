function w = fallvelocity(D, Tw)
    % Calculates dimensionless fall velocity assuming T = 15 deg C
    % Args:
    %   D: d50 in meters
    %   Tw: Water temperature in degrees Celsius (default: 15)
    % Returns:
    %   dimensionless fall velocity in m/s
    
    if nargin < 2
        Tw = 15; % Default temperature
    end
    
    D = D * 100; % Convert diameter to cm
    ROWs = 2.75; % Density of sand
    g = 981; % Gravity in cm/s^2
    
    % Define temperature and corresponding velocity and density arrays
    T = [5, 10, 15, 20, 25];
    v = [0.0157, 0.0135, 0.0119, 0.0105, 0.0095];
    ROW = [1.028, 1.027, 1.026, 1.025, 1.024];
    
    % Interpolate to get vw and ROWw for the given temperature Tw
    vw = interp1(T, v, Tw, 'linear');
    ROWw = interp1(T, ROW, Tw, 'linear');
    
    % Calculate dimensionless fall velocity
    A = ((ROWs - ROWw) * g * (D^3)) / (ROWw * (vw^2));
    
    if A < 39
        w = ((ROWs - ROWw) * g * (D^2)) / (18 * ROWw * vw);
    elseif A < 10^4
        w = ((((ROWs - ROWw) * g / ROWw)^0.7) * (D^1.1)) / (6 * (vw^0.4));
    else
        w = sqrt(((ROWs - ROWw) * g * D) / (0.91 * ROWw));
    end
    
    % Convert to SI units (m/s)
    w = w / 100;
end
