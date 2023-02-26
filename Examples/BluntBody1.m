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
import sggpack.Spline
import sggpack.ClosedSpline

%% Generate Mesh
 
%Create a point set
amplitude = @(t) 1 + 0.2*sin(12*pi*t);
noisy_circle = @(t) [amplitude(t).*cos(2*pi*t),...
                     amplitude(t).*sin(-2*pi*t)];
n = 10;
t = linspace(0,0.8,n)';
r = noisy_circle(t);

%Create a closed curve from the point set
curve = ClosedSpline(r);

%Show Curve
n = 1000;
t = linspace(0,1,n)';
f = curve(t);
figure(1)
set(gcf,'Color','w')
hold on
plot(f(:,1),f(:,2),'-k','LineWidth',2);
plot(r(:,1),r(:,2),'or');
hold off
grid on
xlabel('$x$','Interpreter','latex','Fontsize',14)
ylabel('$y$','Interpreter','latex','Fontsize',14)
legend({'curve', 'point set'},'Interpreter','latex','Fontsize',14)
title('Blunt Body','Interpreter','latex','Fontsize',20)


%Create an instance of the OGridGenerator
mesher = OGridGenerator(@curve.interpolate,...   %Parametric curve for body (pass function handle)
                        'DomainRadius',5,...     %Radius of outer boundary
                        'BodyNodes',   80,...    %Number of nodes along body
                        'RadialNodes', 30);      %Number of nodes between body and outer boundary
                       
%Generate mesh
mesher.mesh();

%Display mesh in figure 2
fig = mesher.show_mesh(2); 



