%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Square.m
%
% Written by: Daniel Fust
%
% This script demonstrates example usage of the OGridGenerator on a 
% square cross section.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Reset Workspace
close all
clear
clc

import sggpack.OGridGenerator

%% Generate Mesh

%Create an instance of the OGridGenerator
mesher = OGridGenerator(@square,...                  %Parametric curve for body (pass function handle)
                        'DomainRadius',5,...         %Radius of outer boundary
                        'BodyNodes',   80,...        %Number of nodes along body
                        'RadialNodes', 30,...        %Number of nodes between body and outer boundary
                        'Verbose',     true,...      %Report on convergence
                        'Pins',[0.0,0.25,0.5,0.75]); %Pin nodes to the corners (t=0.0,0.25,0.5,0.75) to capture sharp features
%Generate mesh
mesher.mesh();

%Display mesh
mesher.show_mesh();




%% Local Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function r = square(t)
%
% Parametric equation in 't' for a square cross section.
%
% Parameters:
%     t         Parametric coordinate 0<t<1
% Returns:
%     r   (2,1) [x,y] coordinates of airfoil surface parameterized by t
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function r = square(t)
    a = 1; r = zeros(1,2);
    if (t < 0.25)
        r(1) = -4*a*t + a/2;
        r(2) = -a/2;
    elseif (t < 0.5)
        r(1) = -a/2;
        r(2) = 4*a*(t-0.25) - a/2;
    elseif(t < 0.75)
        r(1) = 4*a*(t-0.5) - a/2;
        r(2) = a/2;
    else
        r(1) = a/2;
        r(2) = -4*a*(t-0.75) + a/2;
    end
end
