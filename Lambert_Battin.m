% Battin's Lambert's Problem Solver

%...................................................................................................

%                                        DESCRIPTION:

% This function uses Battin's for solving Lambert's Problem. Given two orbital positions and a time
% of flight, the terminal velocities are calculated. The method is more robust than other older
% methhods and does not suffer from above 180 degree transfers. The use of continued fractions is
% present and the accuracy chosen should be sufficient. 

% The codes is an adaptation of the methodology featured in:
% Vallado (2013), Lambert Solution -- Battin Method, "Fundamentals of Astrodynamics and 
% Applications, 4 ed.", 7.6:493-497

%...................................................................................................

%                                        INFORMATION:

% Title: Lambert_Battin
% Author: Oliver K. Morrison (adaptation)
% Date: Fri, December 27, 2021
% Email: oliverkwii@gmail.com

%...................................................................................................

%                                        PARAMETERS:

% Variable  | Description                                      | Units      | Dimensions

% Inputs:
% r_vec_i   | Orbital Position Vector at i (i = 1,2)           | {km}       | [1x3]
% TOF       | Total Time of Flight of Trajectory               | {s}        | [-]
% t_m       | Short-way (=1) or Long-way (=-1)                 | {-}        | [-]
% mu        | Standard Gravitational Parameter of Orbited Body | {km^3/s^2} | [-]

% Outputs:
% v_vec_i   | Orbital Velocity Vectors at i (i = 1,2)          | {km/s}     | [1x3]
% iter      | Number of Iterations Required                    | {-}        | [-]

% Relevant Terms:
% Dnu       | Change in True Anomaly                           | {rad}      | [-]
% r_op      | Parabolic Mean Point Radius                      | {km}       | [-]
% a         | Semi-Major Axis                                  | {km}       | [-]
% f,g,g_dot | Kepler Functions                                 | {-}        | [-]

%...................................................................................................

%                                        FUNCTIONS:

% Function       | Description                        | Inputs                  | Outputs

% Main Function:
% Lambert_Battin | Finds Terminal Velocity Vectors    | [r_vec_i, TOF, t_m, mu] | [v_vec_i, iter]

% Helper Functions:
% Findksi        | Computes Continued Fraction ksi(x) | [eta, x]                | [ksi]
% FindK          | Computes Continued Fraction K(U)   | [U]                     | [K]

%...................................................................................................
                                         
%                                          NOTES:

% It is possible to obtain better results through two methods: 
% (1) Calculate the continued fraction out to further places.
% (2) Solve for y through the cubic: y^3 - y^2 - h_1*y - h_2 = 0

% For this method, elliptical orbits are assumed, however, you can change the initial value of x to
% 0 to converge to a hyperbolic solution if that is desired. 

%...................................................................................................


%===================================================================================================

%                                        MAIN FUNCTION

function [v_vec_1, v_vec_2, iter] = Lambert_Battin(r_vec_1, r_vec_2, TOF, t_m, mu)
    
%---------------------------------------------------------------------------------------------------
%                                1 - Compute Useful Constants
    
    r_1 = norm(r_vec_1);
    r_2 = norm(r_vec_2);
    
    % Change in True Anomaly
    cos_Dnu = dot(r_vec_1, r_vec_2)/(r_1*r_2);
    sin_Dnu = t_m*sqrt(1 - cos_Dnu^2);
    Dnu = atan2(sin_Dnu, cos_Dnu);
    
    % Chords
    c = sqrt(r_1^2 + r_2^2 - 2*r_1*r_2*cos_Dnu);
    s = (r_1 + r_2 + c)/2;
    epsilon = (r_2 - r_1)/r_1;
    
    % Parabolic Mean Point and Angle
    tan_sqr_2w = ((epsilon^2)/4)/(sqrt(r_2/r_1) + (r_2/r_1)*(2 + sqrt(r_2/r_1)));
    r_op = sqrt(r_1*r_2)*(cos(Dnu/4)^2 + tan_sqr_2w);
    
%---------------------------------------------------------------------------------------------------
%                        2 - Iterate on x and y to find Semi-Major Axis
    
    % Loop Setup Variables
    if Dnu > 0 && Dnu < pi
        
        l = (sin(Dnu/4)^2 + tan_sqr_2w)/(sin(Dnu/4)^2 + tan_sqr_2w + cos(Dnu/2));
        
    else
        
        l = (cos(Dnu/4)^2 + tan_sqr_2w - cos(Dnu/2))/(cos(Dnu/4)^2 + tan_sqr_2w);
        
    end
    
    m = (mu*(TOF^2))/(8*(r_op^3));
    
    x = l; % Let x = l if the orbit is elliptical, 0.0 otherwise
    err = 1;
    iter = 0;
    
    % Iterate x and y
    while err > 1e-8
        
        eta = x/((sqrt(1 + x) + 1)^2);
        ksi = Findksi(eta, x);
        
        h_1 = (((l + x)^2)*(1 + 3*x + ksi))/((1 + 2*x + l)*(4*x + ksi*(3 + x)));
        h_2 = (m*(x - l + ksi))/((1 + 2*x + l)*(4*x + ksi*(3 + x)));
        
        B = (27*h_2)/(4*(1 + h_1)^3);
        U = B/(2*(sqrt(1 + B) + 1));
        K = FindK(U);
        
        y = ((1 + h_1)/3)*(2 + sqrt(1 + B)/(1 + 2*U*K^2));
        x_new = sqrt(((1 - l)/2)^2 + m/(y^2)) - (1 + l)/2;
        
        err = abs(x_new - x);
        x = x_new;
        iter = iter + 1;
        
    end
    
    % Semi-Major Axis
    a = (mu*(TOF^2))/(16*(r_op^2)*x*(y^2));
    
%---------------------------------------------------------------------------------------------------
%               3 - Determine Elliptic or Hyperbolic and Find Terminal Velocities
    
    % If Elliptical
    if a > 0
        
        % Use Minimum Energy Method
        alpha_min = pi;
        beta_min = 2*asin(sqrt((s - c)/s));
        
        if Dnu > pi
            
            beta_min = -beta_min;
            
        end
        
        a_min = s/2;
        t_min = sqrt((a_min^3)/mu)*(alpha_min - beta_min + sin(beta_min));
        
        % Intermediate Angles
        alpha_e = 2*asin(sqrt(s/(2*a)));
        beta_e = 2*asin(sqrt((s - c)/(2*a)));
        
        if TOF > t_min
            
            alpha_e = 2*pi - alpha_e;
            
        end
        
        DE = alpha_e - beta_e;
        
        % Kepler's Functions
        f = 1 - (a/r_1)*(1 - cos(DE));
        g = TOF - sqrt((a^3)/mu)*(DE - sin(DE));
        g_dot = 1 - (a/r_2)*(1 - cos(DE));
    
    % If Hyperbolic
    else
        
        % Intermediate Angles
        alpha_h = 2*asinh(sqrt(s/(-2*a)));
        beta_h = 2*asinh(sqrt((s - c)/(-2*a)));
        DH = alpha_h - beta_h;
        
        % Kepler's Functions
        f = 1 - (a/r_1)*(1 - cosh(DH));
        g = TOF - sqrt((-a^3)/mu)*(sinh(DH) - DH);
        g_dot = 1 - (a/r_2)*(1 - cosh(DH));
        
    end
    
    % Terminal Velocities
    v_vec_1 = (r_vec_2 - f*r_vec_1)/g;
    v_vec_2 = (g_dot*r_vec_2 - r_vec_1)/g;
    
%---------------------------------------------------------------------------------------------------
    
end


%===================================================================================================

%                                         HELPER FUNCTIONS

%---------------------------------------------------------------------------------------------------
%                            1 - Calculate ksi(x) | Continued Fraction

function ksi = Findksi(eta, x)
    
    
    ksi_old = 1;
    
    % Iterate on Fraction
    for n = 100:4
        
        c_eta = (n^2)/((2*n)^2 - 1);
        ksi_new = 1 + c_eta*eta/ksi_old;
        ksi_old = ksi_new;
        
    end
    
    ksi = 8*(sqrt(1 + x) + 1)/(3 + 1/(5 + eta + ((9/7)*eta)/(1 + ((16/63)*eta)/ksi_old)));
    
end

%---------------------------------------------------------------------------------------------------
%                            2 - Calculate K(U) | Continued Fraction

function K = FindK(U)
    
    K_old = 1;
    
    % Iterate on Fraction
    for n = 100:0
        
        if mod(n, 2) == 0
            
            c_U = (2*(3*n + 2)*(6*n + 1))/(9*(4*n + 1)*(4*n + 3));
            
        else
            
            c_U = (2*(3*n + 1)*(6*n - 1))/(9*(4*n - 1)*(4*n + 1));
            
        end
        
        K_new = 1 + c_U*U/K_old;
        K_old = K_new;
        
    end
    
    K = (1/3)/(1 + ((4/27)*U)/(1 + ((8/27)*U)/(1 + K_old)));
    
end

%---------------------------------------------------------------------------------------------------

%===================================================================================================


        
        
    
    
        
        
    
    
    
    
        
        
            
            
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    