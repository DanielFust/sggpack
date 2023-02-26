%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file is part of sggpack
%
% sggpack is free software: you can redistribute it and/or modify it 
% under the terms of the GNU General Public License as published by the 
% Free Software Foundation, either version 3 of the License, or (at your 
% option) any later version.
% This program is distributed in the hope that it will be useful, but 
% WITHOUT ANY WARRANTY; without even the implied warranty of 
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General 
% Public License for more details.
% You should have received a copy of the GNU General Public License along 
% with this program. If not, see <https://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% classdef Spline < handle
%
% Written by: Daniel Fust
%
% Class representing a cubic spline curve in R2. Provides parametric
% mapping [x(t),y(t)] on 0<t<1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef Spline < handle
    
    properties (Access = private)
        N
        a
        b
        c
        d
    end
    
    methods (Access = public)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % function obj = Spline(r)
        %
        % Creates an open cubic spline in 2D space
        %
        % Parameters:
        %    r (:,1)    [x,y] coordinates of interpolated points
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = Spline(r)
            arguments
                r (:,2) 
            end
            
            %Initialize member variables/arrays
            [obj.N,~] = size(r);
            obj.a = zeros(obj.N,2);
            obj.b = zeros(obj.N,2);
            obj.c = zeros(obj.N,2);
            obj.d = zeros(obj.N,2);
            
            %Compute coefficient values
            [obj.a(:,1), obj.b(:,1), obj.c(:,1), obj.d(:,1)] ...
                = spline_coefficients(r(:,1));
            [obj.a(:,2), obj.b(:,2), obj.c(:,2), obj.d(:,2)] ...
                = spline_coefficients(r(:,2));
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % function obj = interpolate(obj,t)
        %
        % Interpolates along the curve
        %
        % Parameters:
        %     t         Parametric value of curve
        % Returns:
        %     r (:,2)   [x,y] coordinates at parametric coordinate t
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function r  = interpolate(obj,t)
            arguments
                obj
                t   (:,1) double
            end
            
            n = length(t);
            r = zeros(n,2);
            for k=1:n
                %Index and Internal parametric value
                dist = t(k)*(obj.N-1);
                i = floor(dist)+1;
                z = dist-(i-1);

                x = obj.a(i,1) + obj.b(i,1)*z + obj.c(i,1)*z^2 + obj.d(i,1)*z^3;
                y = obj.a(i,2) + obj.b(i,2)*z + obj.c(i,2)*z^2 + obj.d(i,2)*z^3;

                r(k,:) = [x,y];
            end 
        end
        
        
        
    
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % function obj = subsref(obj,t)
        %
        % operator() overload
        %
        % Parameters:
        %     t         Parametric value of curve
        % Returns:
        %     r (:,2)   [x,y] coordinates at parametric coordinate t
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function varargout = subsref(obj,s)
            arguments
                obj
                s 
            end
            
            switch s(1).type
                case '()'
                    t = s(1).subs{:};
                    [varargout{1:nargout}] = obj.interpolate(t);
                otherwise
                    [varargout{1:nargout}] = builtin('subsref',obj,s);
            end
        end
        
    end %Methods
        
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    Helper Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [a,b,c,d] = spline_coefficients(f)
%
% Computes spline coefficients for a 1D cubic spline
%
% Parameters:
%    r (:,1)    [x,y] coordinates of interpolated points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [a,b,c,d] = spline_coefficients(f)
    N = length(f);
    a = ones(N,1);
    b = 4*ones(N,1); b(1)=2; b(end)=2;
    c = ones(N,1);
    d = zeros(N,1);
    
    rhs = zeros(N,1);
    rhs(2:N-1) = 3*(f(3:N)-f(1:N-2));
    rhs(1) = 3*(f(2)-f(1)); 
    rhs(N) = 3*(f(N)-f(N-1));
    
    b = thomas3(a,b,c,rhs);
    a = f;
    
    c(1:N-1) = 3*(f(2:N)-f(1:N-1)) - 2*b(1:N-1) - b(2:N);
    d(1:N-1) = 2*(f(1:N-1)-f(2:N)) + b(1:N-1) + b(2:N);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function x = thomas3( a, b, c, d )
%
% Direct tridiagonal matrix solver
%
% Parameters:
%     a (N,1)    Lower diagonal
%     b (N,1)    Diagonal
%     c (N,1)    Upper diagonal
%     d (N,1)    Right hand side
% Returns:
%     x         Solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x = thomas3( a, b, c, d )
    n = length(b);
    x = zeros(n,1);
    cp = zeros(n,1);
    dp = zeros(n,1);

    cp(1) = c(1)/b(1);
    for i=2:(n-1)
      cp(i) = c(i)/(b(i) - a(i)*cp(i-1));
    end

    dp(1) = d(1)/b(1);
    for i=2:n
     dp(i) = (d(i) - a(i)*dp(i-1))/(b(i) - a(i)*cp(i-1));
    end

    x(n) = dp(n);
    for i=(n-1):-1:1
     x(i) = dp(i)-cp(i)*x(i+1);
    end 
end



