%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NACA2412.m
%
% Written by: Daniel Fust
%
% This script demonstrates example usage of the OGridGenerator on a 
% NACA2412 airfoil. The OGridGenerator creates a structured mesh based on
% elliptic grid generation.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Reset Workspace
close all
clear
clc

import sggpack.OGridGenerator
%% Generate Mesh

%Define Airfoil
naca2412 = @(z) naca_airfoil(z,0.02,0.4,0.12,1);


%Create an instance of the OGridGenerator
mesher = OGridGenerator(naca2412,...          %Parametric curve for body
                        'DomainRadius',3,...  %Radius of outer boundary
                        'BodyNodes', 150,...  %Number of nodes along body
                        'RadialNodes', 50,... %Number of nodes between body and outer boundary
                        'Verbose',true,...    %Report on convergence
                        'Relaxation',0.8,...  %Solver relaxation parameter (between 0 and 2). Decrease for greater stability, increase for covergence speed. 
                        'Pins',[0.0,0.5]);    %Pin nodes to points (as specified by parametric coordinate) on the body to capture important features (e.g. here t=0.0 is trailing edge, t=0.5 is leading edge)
%Generate mesh
mesher.mesh();

%Display mesh
mesher.show_mesh();




%% Local Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function r = naca_airfoil(z,y_mc,x_mc,t_m,c)
%
% Parametric equation in 'z' for NACA series airfoil.
%
% Parameters:
%     z         Parametric coordinate 0<z<1
%     y_mc      Maximum camber (absolute not percent)
%     x_mc      Location of maximum camber (absolute not tenths)
%     t_m       Maximum thickness (absolute not percent)
%     c         Cord length
% Returns:
%     r   (2,1) [x,y] coordinates of airfoil surface parameterized by z
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function r = naca_airfoil(z,y_mc,x_mc,t_m,c)
    arguments
        z    (1,1) double
        y_mc (1,1) double
        x_mc (1,1) double
        t_m  (1,1) double
        c    (1,1) double
    end
    
    %Lower surface
    if (z<=0.5)
        x = -2*c*z + c;
    else %Upper surface
        x = 2*c*(z-0.5);
    end
    
    %aArgument given relative to chord length
    x_mc = x_mc*c;
    
    %Camber line
    if (x < x_mc)
        y_c = y_mc * ((2*(x/x_mc) - (x/x_mc)^2));
    else
        y_c = y_mc * ((2*((c-x)/(c-x_mc)) - ((c-x)/(c-x_mc))^2));
    end
    
    %Thickness
    t = t_m*(2.969*sqrt(x/c) - 1.260*x/c - 3.516*(x/c).^2 ...
           + 2.843*(x/c).^3 - 1.015*(x/c).^4);
    if (z < 0.5)
        r = [x-0.5*c, -0.5*t + y_c];
    else
        
        r = [x-0.5*c, 0.5*t + y_c];
    end
    
end

