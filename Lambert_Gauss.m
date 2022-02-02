% Gauss' Lambert's Problem Solver

%...................................................................................................

%                                        DESCRIPTION:

% This function uses Gauss' method for solving Lambert's Problem. Given two orbital
% positions and a time between them, the terminal velocities are calculated. This method is
% rudamentary and is only capable of simple solutions. However, it is capable of both elliptic and
% hyperbolic differentiation.

% The codes is an adaptation of the methodology featured in:
% Vallado (2013), Lambert -- Gauss's Solution, "Fundamentals of Astrodynamics and Applications, 
% 4 ed.", 7.6:475-479

%...................................................................................................

%                                        INFORMATION:

% Title: Lambert_Gauss
% Author: Oliver K. Morrison (adaptation)
% Date: Fri, December 24, 2021
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
% p         | Semi-Parameter                                   | {km}       | [-]
% f,g,g_dot | Kepler Functions                                 | {-}        | [-]

%...................................................................................................

%                                        FUNCTIONS:

% Function      | Description                     | Inputs                  | Outputs

% Main Function:
% Lambert_Gauss | Finds Terminal Velocity Vectors | [r_vec_i, TOF, t_m, mu] | [v_vec_i, iter]

% Helper Function: 
% Findx2        | Find Kepler's Common Terms      | [x_1]                   | [x_2]


%...................................................................................................
                                         
%                                          NOTES:

% Since Gauss's method is geometric, solutions over 180 degrees are difficult to calculate and the
% provided method only converges "well" for vectors that are less than 90 degrees apart. This limits
% the solution space considerably. 

%...................................................................................................


%===================================================================================================

%                                        MAIN FUNCTION

function [v_vec_1, v_vec_2, iter] = Lambert_Gauss(r_vec_1, r_vec_2, TOF, t_m, mu)
    
%---------------------------------------------------------------------------------------------------
%                        1 - Define delta-nu and Temporary Variables
    
    r_1 = norm(r_vec_1);
    r_2 = norm(r_vec_2);
    
    % Change in True Anomaly
    cos_Dnu = dot(r_vec_1, r_vec_2)/(r_1*r_2);
    sin_Dnu = t_m*sqrt(1 - cos_Dnu^2);
    Dnu = atan2(sin_Dnu, cos_Dnu);
    
    % Temporary Variables
    l = (r_1 + r_2)/(4*sqrt(r_1*r_2)*cos(Dnu/2)) - 1/2;
    m = (mu*(TOF^2))/((2*sqrt(r_1*r_2)*cos(Dnu/2))^3);

%---------------------------------------------------------------------------------------------------
%                                      2 - Iterate on y
   
    % Loop Setup
    y = 1;
    err = 1;
    iter = 0;
    
    % Find x_1, x_2 until y Stops Changing
    while err > 1e-8
        
        x_1 = m/(y^2) - l;
        x_2 = Findx2(x_1);
        y_new = 1 + x_2*(l + x_1);
        
        err = abs(y_new - y);
        y = y_new;
        iter = iter + 1;
        
    end
    
%---------------------------------------------------------------------------------------------------
%                              3 - Calculate Terminal Velocities
    
    % Semi-Parameter
    cos_DE2 = 1 - 2*x_1;
    p = (r_1*r_2*(1 - cos_Dnu))/(r_1 + r_2 - 2*sqrt(r_1*r_2)*cos(Dnu/2)*cos_DE2);
    
    % Use Kepler's Functions to Solve for Velocities
    f = 1 - (r_2/p)*(1 - cos_Dnu);
    g = (r_1*r_2*sin_Dnu)/sqrt(mu*p);
    g_dot = 1 - (r_1/p)*(1 - cos_Dnu);
    
    v_vec_1 = (r_vec_2 - f*r_vec_1)/g;
    v_vec_2 = (g_dot*r_vec_2 - r_vec_1)/g;
    
%---------------------------------------------------------------------------------------------------
    
end


%===================================================================================================

%                                         HELPER FUNCTION

%---------------------------------------------------------------------------------------------------
%                                 Finding Temporary Variable

function x_2 = Findx2(x_1)
    
    x_2_array = zeros(100, 1);
    x_2_array(1) = 4/3;
    
    % Define Series for Summation
    for i = 2:100
        
        x_2_array(i) = x_2_array(i - 1)*x_1*((2*i + 2)/(2*i + 1));
        
    end
    
    x_2 = sum(x_2_array);
    
end

%---------------------------------------------------------------------------------------------------

%===================================================================================================


        
        
    