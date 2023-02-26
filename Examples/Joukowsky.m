%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Joukowski.m
%
% Written by: Daniel Fust
%
% This script demonstrates example usage of the OGridGenerator on a 
% Joukosky airfoil.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Reset Workspace
close all
clear
clc

import sggpack.OGridGenerator
%% Generate Mesh
N=1000;
t = linspace(0,1,N);
joukowsky_airfoil = @(t) joukowsky_transform(t, 1.25, [-0.2, 0.2]);



%%
%Create an instance of the OGridGenerator
mesher = OGridGenerator(joukowsky_airfoil,... %Parametric curve for body
                        'DomainRadius',5,...  %Radius of outer boundary
                        'BodyNodes', 80,...   %Number of nodes along body
                        'RadialNodes', 30,... %Number of nodes between body and outer boundary
                        'Relaxation',0.8,...  %Solver relaxation parameter (between 0 and 2). Decrease for greater stability, increase for covergence speed. 
                        'Pins',0.0);          %Pin nodes to points (as specified by parametric coordinate) on the body to capture important features (e.g. here t=0.0 is trailing edge)
%Generate mesh
mesher.mesh();

%Display mesh
mesher.show_mesh();




%% Local Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function r = joukowsky_transform(t,R,zeta0)
%
% Parametric equation in 'z' for NACA series airfoil.
%
% Parameters:
%     t           Parametric coordinate 0<t<1
%     R           Radius of polar coordinate in untransformed coordinates
%     zeta0 (1,2) Offset of origin in complex plane
% Returns:
%     r     (1,2) [x,y] coordinates of transformed point
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function r = joukowsky_transform(t,R,zeta0)
    arguments
        t     (1,1) double
        R     (1,1) double
        zeta0 (1,2) double
    end
    
    theta = -2*pi*t;
    i = sqrt(-1);
    zeta = R*exp(i*theta) + zeta0(1) + i*zeta0(2);
    
    z = zeta + 1/zeta;
    r = [real(z), imag(z)];
end

