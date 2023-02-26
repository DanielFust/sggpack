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
% classdef OGridGenerator < handle
%
% Written by: Daniel Fust
%
% OGridGenerator is a body-fitted 2D structured O-grid generator based on
% elliptic mesh generation techniques. An instance of OGridGenerator is
% created by passing a function handle to a parametric curve r(t), where
% r(t) = [x(t), y(t)], for 0<t<1 and by specifying any desired solver and 
% grid parameters:
%       BodyNodes...................Nodes on body
%       RadialNodes.................Nodes between body and boundary
%       DomainRadius................Radius of outer domain 
%       Tolerance...................Solution tolerance
%       MaxIterations...............Maximum allowed iterations
%       Relaxation..................Over or Under-relaxation parameter
%       Verbose.....................Toggle for some warnings and messages
%       ReportFrequency.............Iterations between covergence reporting
%       Pins........................Parametric locations to strictly
%                                   enforce body conformity without regard
%                                   for orthogonality (recommended for 
%                                   sharp features)
%
% Example usages may be found in the 'Examples' directory
%    
%
% References:
%     Thompson, J. F.; Warsi, Z. U. A.; Mastin, C. W. (1985), Numerical 
%     Grid Generation: Foundations and Applications, North-Holland, 
%     Elsevier.
%
%    Joe F. Thompson, Elliptic grid generation, Applied Mathematics and 
%    Computation, Volumes 10–11, 1982, Pages 79-105, ISSN 0096-3003,
%    https://doi.org/10.1016/0096-3003(82)90188-6.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef OGridGenerator < handle
   
    
    properties (Access = private)
        imax
        jmax
        body_curve = @(t) 0;
        tolerance
        max_iterations
        relaxation
        x
        y
        radius
        verbose
        report_freq
        theta0
        pins
    end
    
    properties (Constant, Access = private)
        % Basically an enumeration for utility, but Matlab syntax makes you
        % either define a different class which I don't want to do, frankly
        % or not restrict access. Just don't touch these
        x_equation = 1;
        y_equation = 2;
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                        Public Methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Access = public)
        
        
        
        function obj = OGridGenerator(body_contour, options)
            arguments
                body_contour      (1,1) function_handle
                options.BodyNodes (1,1) double = 40
                options.RadialNodes (1,1) double = 20
                options.DomainRadius      (1,1) double = 1.0 
                options.Tolerance   (1,1) double = 1.0E-05 
                options.MaxIterations (1,1) double = 100
                options.Relaxation (1,1) double = 0.9
                options.Verbose    (1,1) logical = false
                options.ReportFrequency (1,1) double = 1;
                options.Pins             (1,:) double = [];
            end
            
            obj.imax = floor(options.BodyNodes);
            obj.jmax = floor(options.RadialNodes);
            obj.radius = options.DomainRadius;
            obj.tolerance = options.Tolerance;
            obj.max_iterations = options.MaxIterations;
            obj.relaxation = options.Relaxation;
            obj.verbose = options.Verbose;
            obj.report_freq = options.ReportFrequency;
            
            obj.set_body_curve(body_contour);
            obj.set_pins(options.Pins);
        end
        
        
        function obj = set(obj,options)
            arguments
                obj
                options.BodyNodes       (1,1) double  = obj.imax
                options.RadialNodes     (1,1) double  = obj.jmax
                options.DomainRadius    (1,1) double  = obj.radius 
                options.Tolerance       (1,1) double  = obj.tolerance 
                options.MaxIterations   (1,1) double  = obj.max_iterations
                options.Relaxation      (1,1) double  = obj.relaxation
                options.Verbose         (1,1) logical = obj.verbose
                options.ReportFrequency (1,1) double  = obj.report_freq;
                options.Pins            (1,:) double  = obj.pins(:,1)
                options.Body            (1,1) function_handle = obj.body_curve;
            end
            
            obj.imax           = options.BodyNodes;
            obj.jmax           = options.RadialNodes;
            obj.radius         = options.DomainRadius;
            obj.tolerance      = options.Tolerance;
            obj.max_iterations = options.MaxIterations;
            obj.verbose        = options.Verbose;
            obj.report_freq    = options.ReportFrequency;
            
            
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                    Getter/Setter Methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = set_body_nodes(obj,body_nodes)
            obj.imax = body_nodes;
        end
        
        function obj = set_radial_nodes(obj,radial_nodes)
            obj.jmax = radial_nodes;
        end
        
        function obj = set_domain_radius(obj,domain_radius)
            obj.radius = domain_radius;
        end
        
        function obj = set_tolerance(obj,tol)
            obj.tolerance = tol;
        end
        
        function obj = set_body_curve(obj, curve)
            arguments
                obj
                curve (1,1) function_handle
            end
            %Set body curve
            obj.body_curve = curve;
            r = obj.body(0);
            obj.theta0 = atan2(r(2),r(1));
        end
        
        function obj = set_pins(obj,pin_locations)
            arguments
                obj
                pin_locations (:,1) double
            end
            %Initialize Pins
            [n_pins, ~] = size(pin_locations);
            if (n_pins > 0)
                obj.pins = zeros(n_pins,2);
                obj.pins(:,1) = pin_locations(:);
            end
        end
        
        function obj = set_verbose(obj,is_verbose)
            arguments
                obj
                is_verbose (1,1) logical
            end
            obj.verbose = is_verbose;
        end
        
        function obj = set_report_frequency(obj,freq)
            obj.report_freq(freq);
        end
        
        function [X,Y] = get_mesh(obj)
            X = obj.x;
            Y = obj.y;
        end
        
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % function [x,y] = body(obj,t)
        % 
        % Calculates [x,y] coordinates for a parameter value 0<t<1,
        % accounting for wrapping. 
        % 
        % Parameters:
        %     t      Parametric coordinate for curve
        % Returns:
        %     r     (2,1)    [x,y] Coordinates of body curve at 't'
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function r = body(obj,t)
            if(t>1)
                t=t-floor(t);
            elseif(t<0)
                t=1-ceil(t)+t;
            end
            r = obj.body_curve(t);
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % function [x,y] = boundary(obj,t)
        % 
        % Calculates [x,y] coordinates for a parameter value 0<t<1,
        % accounting for wrapping. 
        % 
        % Parameters:
        %     t      Parametric coordinate for curve
        % Returns:
        %     [x,y]  Coordinates of boundary curve at 't'
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function r = boundary(obj,t)
            r = zeros(1,2);
            r(1) =  obj.radius * cos(2*pi*t - obj.theta0);
            r(2) = -obj.radius * sin(2*pi*t - obj.theta0);
        end
        
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % function [X,Y] = mesh(obj)
        %
        % This method is called to compute the mesh based on the set
        % parameters of the solver
        % 
        % Parameters:
        % Returns:
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [X,Y] = mesh(obj)
            import sggpack.SIS
            
            %Set aliases because Matlab OOP syntax sucks
            m = obj.imax;
            n = obj.jmax;
            N = m*n;
            
            %Solver parameters
            tol = obj.tolerance;
            iter_max = obj.max_iterations;
            beta = obj.relaxation;
            
            %Initialize x,y and contour curves with reasonable guess
            [t_body,t_boundary] = obj.initial_guess();
            s_body(:,:)     = [obj.x(:,1), obj.y(:,1)];
            s_boundary(:,:) = [obj.x(:,n), obj.y(:,n)];
            
            %obj.show_mesh(1);
            
            %Solver loop
            res = inf;
            iter = 0;
            while (res > tol)
                
                %Generate x-system
                [a_sw, a_w, a_nw, a_s, a_p, a_n, a_se, a_e, a_ne, b_p] ...
                = obj.build_linear_system(s_body,s_boundary,obj.x_equation);
            
                %Compute new x values with relaxation and measure convergence 
                f_new = SIS(a_sw,a_s,a_se,a_w,a_p,a_e,a_nw,a_n,a_ne,obj.x(:),b_p,n,m);
                res = norm(f_new-obj.x(:))/N;
                obj.x(:) = beta*f_new + (1-beta)*obj.x(:);
                
                %obj.show_mesh(1);
                
                %Estimate new boundary points based on updated x values
                if (res < 10*tol)  
                    for i=1:m %i=1 is always pinned
%                         r = [obj.x(i,1), obj.y(i,1)];
%                         t_body(i) = obj.closest_point(r, @obj.body, t_body(i));
                     
                        r = [obj.x(i,n), obj.y(i,n)];
                        t_boundary(i) = obj.closest_point(r, @obj.boundary, t_boundary(i));
                     
%                         s_body(i,:) = obj.body(t_body(i));
                        s_boundary(i,:) = obj.boundary(t_boundary(i));
                    end  
                end
                       
            
                %Generate y-system
                [a_sw, a_w, a_nw, a_s, a_p, a_n, a_se, a_e, a_ne, b_p] ...
                = obj.build_linear_system(s_body,s_boundary,obj.y_equation);
            
                %Compute new x values with relaxation and measure convergence 
                f_new = SIS(a_sw,a_s,a_se,a_w,a_p,a_e,a_nw,a_n,a_ne,obj.y(:),b_p,n,m);
                res = max(norm(f_new-obj.y(:))/N, res);
                obj.y(:) = beta*f_new + (1-beta)*obj.y(:);
               
                %Estimate new boundary points based on updated y values
                if (res < 10*tol)   
                    for i=1:m %i=1 is always pinned
%                         r = [obj.x(i,1), obj.y(i,1)];
%                         t_body(i) = obj.closest_point(r, @obj.body, t_body(i));

                        r = [obj.x(i,n), obj.y(i,n)];
                        t_boundary(i) = obj.closest_point(r, @obj.boundary, t_boundary(i));

%                         s_body(i,:) = obj.body(t_body(i));
                        s_boundary(i,:) = obj.boundary(t_boundary(i));
                    end   
                 end
                
%                  obj.show_mesh(1); 
                  
                %Check for divergence
                if (isnan(res)||isinf(res))
                    warning("Solution diverged. Failed to generate mesh.")
                    return
                end
                
                %Iteration Cap
                iter = iter+1;
                if (iter >= iter_max)
                    fprintf("Maximum iterations reached. Aborting mesher with residual = %d",res);
                    break;
                end
                
                
                %Report convergence progress to user
                if (mod(iter,obj.report_freq)==0)
                    fprintf("Iteration: %i, Residual: %d\n",iter,res);
                end
                
            end
            
            %Ensure that boundary points fall exactly on body/perimeter
            
            for i=1:m 
                r = [obj.x(i,1), obj.y(i,1)];
                t_body(i) = obj.closest_point(r, @obj.body, t_body(i));

                r = [obj.x(i,n), obj.y(i,n)];
                t_boundary(i) = obj.closest_point(r, @obj.boundary, t_boundary(i));

                s_body(i,:) = obj.body(t_body(i));
                s_boundary(i,:) = obj.boundary(t_boundary(i));
            end   
                 
            for i=1:m
                obj.x(i,1) = s_body(i,1);
                obj.y(i,1) = s_body(i,2);
                
                obj.x(i,n) = s_boundary(i,1);
                obj.y(i,n) = s_boundary(i,2);
            end
            X=obj.x; Y=obj.y;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % function show_mesh(obj)
        %
        % Plots the mesh. Seam is shown as connected (i.e. as if indices 
        % wrap such that for index i, m+1 = 1)
        %
        % Parameters:
        %     fig_num   Figure number
        % Returns:
        %     fig       figure handle
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function fig = show_mesh(obj,fig_num)
            arguments
                obj
                fig_num (1,1) double = -1
            end
            %If Mesh hasn't been generated, return early
            if (isempty(obj.x) == 1)
                warning("OGridGenerator.show_mesh(). Mesh has not been generated. Call OGridGenerator.mesh() method to generate mesh.");
                return;
            end
            
            %Plotting parameters
            lw = 1; %LineWidth
            text_size = 16;
            title_size = 20;
            
            m = obj.imax;
            n = obj.jmax;
            
            if (fig_num<0)
                fig = figure();
            else
                fig = figure(fig_num);
                clf(fig);
            end
            
            set(gcf,'Color','w')
            hold on
            %Radial Grid Lines
            for i=1:m
                plot(obj.x(i,:),obj.y(i,:),'-k','LineWidth',lw)
            end
            %Circumferential Grid Lines
            for j=1:n
                plot([obj.x(:,j); obj.x(1,j)], ...
                     [obj.y(:,j); obj.y(1,j)],...
                     '-k', 'LineWidth',lw)
            end
            hold off
            xlabel('$x$','Interpreter','latex','FontSize',text_size)
            ylabel('$y$','Interpreter','latex','FontSize',text_size)
            title('Mesh','Interpreter','latex','FontSize',title_size)
            
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                        Private Methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Access = private)
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % function [s_body, s_perimeter] = initial_guess(obj)
        %
        % Initializes and establishes x,y grids. Provides guess for
        % (x,y) locations on body curve and perimeter curve at i indices 
        %
        % Parameters:
        % Returns:
        %     s_body      (m,2) [x,y] coordinates at i indices for body
        %                 curve
        %     s_boundary  (m,2) [x,y] coordinates at i indices for
        %                  perimeter curve
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [t_body, t_boundary] = initial_guess(obj)
            % Initialize grids
            m = obj.imax;
            n = obj.jmax;
            obj.x = zeros(m,n);
            obj.y = zeros(m,n);
            
            % Parametric samples
            t_body     = uniform_sample(@obj.body,m);
            t_boundary = uniform_sample(@obj.boundary,m);
            s_body     = zeros(m,2);
            s_boundary = zeros(m,2);
            
       
            % Correct parametric locations for pins
            [n_pins, ~] = size(obj.pins);
            if (n_pins > 0)
                pin_index = 1;
                prev_pt_pinned = false;
                for i=2:m
                    t_pin = obj.pins(pin_index,1);
                    % If there is a pin between indices
                    if (t_body(i)>t_pin)
                        if (((t_pin - t_body(i-1)) < (t_body(i) - t_pin)) && (~prev_pt_pinned))
                            obj.pins(pin_index,2) = i-1;
                            prev_pt_pinned = false;
                        else
                            obj.pins(pin_index,2) = i;
                            prev_pt_pinned = true;
                        end
                        pin_index = pin_index+1;
                        if (pin_index > n_pins)
                            break;
                        end
                    else
                        prev_pt_pinned = false;
                    end
                end
            end
            
            % Get body and perimeter grid locations
            for i=1:m
                s_body(i,:)     = obj.body(t_body(i));
                s_boundary(i,:) = obj.boundary(t_boundary(i));
            end
            
            % Linearly Interpolate initial grid
            for i=1:m
                r0 = s_body(i,:);
                rf = s_boundary(i,:);
                obj.x(i,:) = linspace(r0(1),rf(1),n);
                obj.y(i,:) = linspace(r0(2),rf(2),n);
            end
        end
        
        
        
          
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % function [a_sw, a_w, a_nw, a_s, a_p, a_n, a_se, a_e, a_ne, b_p]
        %     = build_linear_system(x,y,s_body,s_perimeter)
        %
        % Generate linear system for x for implict 9-point 2D stencil
        %
        % Parameters:
        %     s_body      (m,2) array of (x,y) coordinates of body contour 
        %                 at xi (i-indexed) locations
        %     s_perimeter (m,2) array of (x,y) coordinates of perimeter 
        %                 contour at xi (i-indexed) locations
        %     eqn_number  (int) 1 for x-equation 2 for y-equation
        % Returns:
        %     (m*n,1) coefficient arrays of the form: sum(a*u) = b
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [a_sw, a_w, a_nw, a_s, a_p, a_n, a_se, a_e, a_ne, b_p] ...
                = build_linear_system(obj,s_body,s_perimeter,eqn_number)
            %Set aliases because Matlab OOP syntax sucks
            m = obj.imax;
            n = obj.jmax;
            N = m*n;
            
            %Allocate space for coefficients
            a_sw = zeros(N,1);
            a_w  = zeros(N,1);
            a_nw = zeros(N,1);
            a_s  = zeros(N,1);
            a_p  = zeros(N,1);
            a_n  = zeros(N,1);
            a_se = zeros(N,1);
            a_e  = zeros(N,1);
            a_ne = zeros(N,1);
            b_p  = zeros(N,1);
            
            %Lower O-grid seam
            i=1;
            for j=2:(n-1)
                k=i+m*(j-1);
                %k=j+n*(i-1);
                
                %a_p(k) = 1;
                %b_p(k) = (s_perimeter(i,eqn_number)-s_body(i,eqn_number)) * (j-1) / (n-1) + s_body(i,eqn_number);
                
                x_sw = obj.x(m,j-1); y_sw = obj.y(m,j-1);
                x_w  = obj.x(m,j);   y_w  = obj.y(m,j);
                x_nw = obj.x(m,j+1); y_nw = obj.y(m,j+1);
                
                x_s = obj.x(i,j-1); y_s = obj.y(i,j-1);
                x_p = obj.x(i,j);   y_p = obj.y(i,j);
                x_n = obj.x(i,j+1); y_n = obj.y(i,j+1);
                
                x_se = obj.x(i+1,j-1); y_se = obj.y(i+1,j-1);
                x_e  = obj.x(i+1,j);   y_e  = obj.y(i+1,j);
                x_ne = obj.x(i+1,j+1); y_ne = obj.y(i+1,j+1);
                
                x_loc = [x_sw, x_w, x_nw, x_s, x_p, x_n, x_se, x_e, x_ne]';
                y_loc = [y_sw, y_w, y_nw, y_s, y_p, y_n, y_se, y_e, y_ne]';
                
                % Compute Coefficients
                coefs = get_coefs(x_loc,y_loc);
                a_sw(k) = coefs(1); a_w(k) = coefs(2); a_nw(k) = coefs(3);
                a_s(k)  = coefs(4); a_p(k) = coefs(5); a_n(k)  = coefs(6);
                a_se(k) = coefs(7); a_e(k) = coefs(8); a_ne(k) = coefs(9);
                switch eqn_number
                    case 1
                        b_p(k) = coefs(10) - (a_sw(k)*x_sw + a_w(k)*x_w + a_nw(k)*x_nw); %Lagging the seam BC as it doesn't conform with SIS
                        a_sw(k) = 0; a_w(k) = 0; a_nw(k) = 0;
                    case 2
                        b_p(k) = coefs(10) - (a_sw(k)*y_sw + a_w(k)*y_w + a_nw(k)*y_nw); %Lagging the seam BC as it doesn't conform with SIS
                        a_sw(k) = 0; a_w(k) = 0; a_nw(k) = 0;
                end
            end
            
            
            
            
            %Upper O-grid seam
            i=m;
            for j=2:(n-1)
                k=i+m*(j-1);
                %k=j+n*(i-1);
                
                %a_p(k) = 1;
                %b_p(k) = (s_perimeter(i,eqn_number)-s_body(i,eqn_number)) * (j-1) / (n-1) + s_body(i,eqn_number);
                x_sw = obj.x(i-1,j-1); y_sw = obj.y(i-1,j-1);
                x_w  = obj.x(i-1,j);   y_w  = obj.y(i-1,j);
                x_nw = obj.x(i-1,j+1); y_nw = obj.y(i-1,j+1);
                
                x_s = obj.x(i,j-1); y_s = obj.y(i,j-1);
                x_p = obj.x(i,j);   y_p = obj.y(i,j);
                x_n = obj.x(i,j+1); y_n = obj.y(i,j+1);
                
                x_se = obj.x(1,j-1); y_se = obj.y(1,j-1);
                x_e  = obj.x(1,j);   y_e  = obj.y(1,j);
                x_ne = obj.x(1,j+1); y_ne = obj.y(1,j+1);
                
                x_loc = [x_sw, x_w, x_nw, x_s, x_p, x_n, x_se, x_e, x_ne]';
                y_loc = [y_sw, y_w, y_nw, y_s, y_p, y_n, y_se, y_e, y_ne]';
                
                % Compute Coefficients
                coefs = get_coefs(x_loc,y_loc);
                a_sw(k) = coefs(1); a_w(k) = coefs(2); a_nw(k) = coefs(3);
                a_s(k)  = coefs(4); a_p(k) = coefs(5); a_n(k)  = coefs(6);
                a_se(k) = coefs(7); a_e(k) = coefs(8); a_ne(k) = coefs(9);
                switch eqn_number
                    case 1
                        b_p(k) = coefs(10) - (a_se(k)*x_se + a_e(k)*x_e + a_ne(k)*x_ne); %Lagging the seam BC as it doesn't conform with SIS
                    case 2
                        b_p(k) = coefs(10) - (a_se(k)*y_se + a_e(k)*y_e + a_ne(k)*y_ne); %Lagging the seam BC as it doesn't conform with SIS
                end
                a_se(k) = 0; a_e(k) = 0; a_ne(k) = 0;
            end
            
%             
%             %Body Contour
%             j=1;
%             for i=1:m
%                 k=i+m*(j-1);
%                 
%                 x_p = obj.x(i,j);   y_p = obj.y(i,j);
%                 x_n = obj.x(i,j+1); y_n = obj.y(i,j+1);
%                 x_s  = 2*s_body(i,1) - x_n;    y_s  = 2*s_body(i,2) - y_n;
%                 
%                 %Lower O-grid seam
%                 if (i==1)
%                     x_nw = obj.x(m,j+1); y_nw = obj.y(m,j+1);
%                     x_w  = obj.x(m,j);   y_w  = obj.y(m,j);
%                     x_sw = 2*s_body(m,1) - x_nw; y_sw = 2*s_body(m,2) - y_nw;
%                 else
%                     x_nw = obj.x(i-1,j+1); y_nw = obj.y(i-1,j+1);
%                     x_w  = obj.x(i-1,j);   y_w  = obj.y(i-1,j);
%                     x_sw = 2*s_body(i-1,1) - x_nw; y_sw = 2*s_body(i-1,2) - y_nw;
%                 end
%                 
%                 %Upper O-grid seam
%                 if (i==m)
%                     x_e  = obj.x(1,j);   y_e  = obj.y(1,j);
%                     x_ne = obj.x(1,j+1); y_ne = obj.y(1,j+1);
%                     x_se = 2*s_body(1,1) - x_ne; y_se = 2*s_body(1,2) - y_ne;
%                 else
%                     x_e  = obj.x(i+1,j);   y_e  = obj.y(i+1,j);
%                     x_ne = obj.x(i+1,j+1); y_ne = obj.y(i+1,j+1);
%                     x_se = 2*s_body(i+1,1) - x_ne; y_se = 2*s_body(i+1,2) - y_ne;
%                 end
%                 
%                 
%                 x_loc = [x_sw, x_w, x_nw, x_s, x_p, x_n, x_se, x_e, x_ne]';
%                 y_loc = [y_sw, y_w, y_nw, y_s, y_p, y_n, y_se, y_e, y_ne]';
%                 
%                 % Compute Coefficients
%                 coefs = get_coefs(x_loc,y_loc);
%                 a_sw(k) = coefs(1); a_w(k) = coefs(2); a_nw(k) = coefs(3);
%                 a_s(k)  = coefs(4); a_p(k) = coefs(5); a_n(k)  = coefs(6);
%                 a_se(k) = coefs(7); a_e(k) = coefs(8); a_ne(k) = coefs(9);
%                 b_p(k) = coefs(10);
%                 
%                 
%                 
%                 %Upper O-grid seam
%                 if (i==m)
%                     %b_p(k) = b_p(k) - 2*a_se(k)*s_body(1,eqn_number);
%                     switch eqn_number
%                         case 1
%                             b_p(k) = b_p(k) - (a_se(k)*x_se + a_e(k)*x_e + a_ne(k)*x_ne);
%                         case 2
%                             b_p(k) = b_p(k) - (a_se(k)*y_se + a_e(k)*y_e + a_ne(k)*y_ne);
%                     end
%                     a_se(k)=0; a_e(k)=0; a_ne(k)=0;
%                 else
%                     a_ne(k) = a_ne(k) - a_se(k);
%                     b_p(k) = b_p(k) - 2*a_se(k)*s_body(i+1,eqn_number);
%                     a_se(k) = 0;
%                 end
%                     
%                 %Lower O-grid seam    
%                 if (i==1)
%                     %b_p(k) = b_p(k) - 2*a_sw(k)*s_body(m,eqn_number);
%                     switch eqn_number
%                         case 1
%                             b_p(k) = b_p(k) - (a_sw(k)*x_sw + a_w(k)*x_w + a_nw(k)*x_nw);
%                         case 2
%                             b_p(k) = b_p(k) - (a_sw(k)*y_sw + a_w(k)*y_w + a_nw(k)*y_nw);
%                     end
%                     a_sw(k)=0; a_w(k)=0; a_nw(k)=0;
%                 else
%                     a_nw(k) = a_nw(k) - a_sw(k);
%                     b_p(k) = b_p(k) - 2*a_sw(k)*s_body(i-1,eqn_number);
%                     a_sw(k) = 0;
%                 end
%                 b_p(k) = b_p(k) - 2*a_s(k)*s_body(i,eqn_number);
%                 a_n(k) = a_n(k) - a_s(k); a_s(k) = 0;
%             end



            %Body Contour
            j=1;
            for i=1:m
                k=i+m*(j-1);
                
                x_p = obj.x(i,j);   y_p = obj.y(i,j);
                x_n = obj.x(i,j+1); y_n = obj.y(i,j+1);
                x_s  = 2*s_body(i,1) - x_n;    y_s  = 2*s_body(i,2) - y_n;
                
                %Lower O-grid seam
                if (i==1)
                    x_nw = obj.x(m,j+1); y_nw = obj.y(m,j+1);
                    x_w  = obj.x(m,j);   y_w  = obj.y(m,j);
                else
                    x_nw = obj.x(i-1,j+1); y_nw = obj.y(i-1,j+1);
                    x_w  = obj.x(i-1,j);   y_w  = obj.y(i-1,j);
                end
                
                %Upper O-grid seam
                if (i==m)
                    x_e  = obj.x(1,j);   y_e  = obj.y(1,j);
                    x_ne = obj.x(1,j+1); y_ne = obj.y(1,j+1);
                else
                    x_e  = obj.x(i+1,j);   y_e  = obj.y(i+1,j);
                    x_ne = obj.x(i+1,j+1); y_ne = obj.y(i+1,j+1);
                end
                
                delta = [(y_n-y_s), -(x_n-x_s)];
                x_se = 2*s_body(i,1) + delta(1) - x_ne;
                y_se = 2*s_body(i,2) + delta(2) - y_ne;
                x_sw = 2*s_body(i,1) - delta(1) - x_nw;
                y_sw = 2*s_body(i,2) - delta(2) - y_nw;
                
                x_loc = [x_sw, x_w, x_nw, x_s, x_p, x_n, x_se, x_e, x_ne]';
                y_loc = [y_sw, y_w, y_nw, y_s, y_p, y_n, y_se, y_e, y_ne]';
                
                % Compute Coefficients
                coefs = get_coefs(x_loc,y_loc);
                a_sw(k) = coefs(1); a_w(k) = coefs(2); a_nw(k) = coefs(3);
                a_s(k)  = coefs(4); a_p(k) = coefs(5); a_n(k)  = coefs(6);
                a_se(k) = coefs(7); a_e(k) = coefs(8); a_ne(k) = coefs(9);
                b_p(k) = coefs(10);
                
                
                %Upper O-grid seam
                if (i==m)
                    %b_p(k) = b_p(k) - 2*a_se(k)*s_body(1,eqn_number);
                    switch eqn_number
                        case 1
                            b_p(k) = b_p(k) - (a_se(k)*x_se + a_e(k)*x_e + a_ne(k)*x_ne);
                        case 2
                            b_p(k) = b_p(k) - (a_se(k)*y_se + a_e(k)*y_e + a_ne(k)*y_ne);
                    end
                    a_se(k)=0; a_e(k)=0; a_ne(k)=0;
                else
                    a_ne(k) = a_ne(k) - a_se(k);
                    b_p(k) = b_p(k) - 2*a_se(k)*(s_body(i,eqn_number) + 0.5*delta(eqn_number));
                    a_se(k) = 0;
                end
                    
                %Lower O-grid seam    
                if (i==1)
                    %b_p(k) = b_p(k) - 2*a_sw(k)*s_body(m,eqn_number);
                    switch eqn_number
                        case 1
                            b_p(k) = b_p(k) - (a_sw(k)*x_sw + a_w(k)*x_w + a_nw(k)*x_nw);
                        case 2
                            b_p(k) = b_p(k) - (a_sw(k)*y_sw + a_w(k)*y_w + a_nw(k)*y_nw);
                    end
                    a_sw(k)=0; a_w(k)=0; a_nw(k)=0;
                else
                    a_nw(k) = a_nw(k) - a_sw(k);
                    b_p(k) = b_p(k) - 2*a_sw(k)*(s_body(i,eqn_number) - 0.5*delta(eqn_number));
                    a_sw(k) = 0;
                end
                b_p(k) = b_p(k) - 2*a_s(k)*s_body(i,eqn_number);
                a_n(k) = a_n(k) - a_s(k); a_s(k) = 0;
            end

            
            
            
            
            
            
            %Perimeter Contour
            j=n;
            for i=1:m
                k=i+m*(j-1);
                
                x_p = obj.x(i,j);                y_p = obj.y(i,j);
                x_s = obj.x(i,j-1);              y_s = obj.y(i,j-1);
                x_n  = 2*s_perimeter(i,1) - x_s; y_n  = 2*s_perimeter(i,2) - y_s;
                
                %Lower O-grid seam
                if (i==1)
                    x_sw = obj.x(m,j-1); y_sw = obj.y(m,j-1);
                    x_w  = obj.x(m,j);   y_w  = obj.y(m,j);
                    x_nw = 2*s_perimeter(m,1) - x_sw;   y_nw = 2*s_perimeter(m,2) - y_sw;
                else
                    x_sw = obj.x(i-1,j-1); y_sw = obj.y(i-1,j-1);
                    x_w  = obj.x(i-1,j);   y_w  = obj.y(i-1,j);
                    x_nw = 2*s_perimeter(i-1,1) - x_sw; y_nw = 2*s_perimeter(i-1,2) - y_sw;
                end
                
                %Upper O-grid seam
                if (i==m)
                    x_e  = obj.x(1,j);   y_e  = obj.y(1,j);
                    x_se = obj.x(1,j-1); y_se = obj.y(1,j-1);
                    x_ne = 2*s_perimeter(1,1) - x_se;   y_ne = 2*s_perimeter(1,2) - y_se;
                else
                    x_e  = obj.x(i+1,j);   y_e  = obj.y(i+1,j);
                    x_se = obj.x(i+1,j-1); y_se = obj.y(i+1,j-1);
                    x_ne = 2*s_perimeter(i+1,1) - x_se; y_ne = 2*s_perimeter(i+1,2) - y_se;
                end
                
                
                x_loc = [x_sw, x_w, x_nw, x_s, x_p, x_n, x_se, x_e, x_ne]';
                y_loc = [y_sw, y_w, y_nw, y_s, y_p, y_n, y_se, y_e, y_ne]';
                
                % Compute Coefficients
                coefs = get_coefs(x_loc,y_loc);
                a_sw(k) = coefs(1); a_w(k) = coefs(2); a_nw(k) = coefs(3);
                a_s(k)  = coefs(4); a_p(k) = coefs(5); a_n(k)  = coefs(6);
                a_se(k) = coefs(7); a_e(k) = coefs(8); a_ne(k) = coefs(9);
                b_p(k) = coefs(10);
                
                % Implicit BC Mapping
                a_s(k) = a_s(k) - a_n(k);
                b_p(k) = b_p(k) - 2*a_n(k)*s_perimeter(i,eqn_number);
                a_n(k) = 0;
                %Upper O-grid seam
                if (i==m)
                    %b_p(k) = b_p(k) - 2*a_ne(k)*s_perimeter(1,eqn_number);
                    switch eqn_number
                        case 1
                            b_p(k) = b_p(k) - (a_se(k)*x_se + a_e(k)*x_e  + a_ne(k)*x_ne);
                        case 2
                            b_p(k) = b_p(k) - (a_se(k)*y_se + a_e(k)*y_e  + a_ne(k)*y_ne);
                    end
                    a_se(k)=0; a_e(k)=0; a_ne(k)=0;         
                else
                    a_se(k) = a_se(k) - a_ne(k);
                    b_p(k) = b_p(k) - 2*a_ne(k)*s_perimeter(i+1,eqn_number);
                    a_ne(k) = 0;
                end
                %Lower Seam
                if (i==1)
                    %b_p(k) = b_p(k) - 2*a_nw(k)*s_perimeter(m,eqn_number);
                    switch eqn_number
                        case 1
                            b_p(k) = b_p(k) - (a_sw(k)*x_sw + a_w(k)*x_w  + a_nw(k)*x_nw);
                        case 2
                            b_p(k) = b_p(k) - (a_sw(k)*y_sw + a_w(k)*y_w  + a_nw(k)*y_nw);
                    end
                       a_sw(k)=0; a_w(k)=0; a_nw(k)=0;      
                else
                    a_sw(k) = a_sw(k) - a_nw(k);
                    b_p(k) = b_p(k) - 2*a_nw(k)*s_perimeter(i-1,eqn_number);
                    a_nw(k) = 0;
                end
                
                
            end
            
            
            %Pins
            [n_pins,~] = size(obj.pins);
            j=1;
            for p=1:n_pins
                t = obj.pins(p,1);
                i = obj.pins(p,2);
                
                k=i+m*(j-1);
                r = obj.body(t);
                a_p(k)=1; b_p(k) = r(eqn_number);
                a_sw(k)=0; a_w(k)=0; a_nw(k)=0; a_s(k)=0; a_n(k)=0; 
                a_se(k)=0; a_e(k)=0; a_ne(k)=0;
            end
            
            


            %Interior Points
            for i=2:(m-1)
                for j=2:(n-1)
                    k=i+m*(j-1);
                    %k=j+n*(i-1);

                    x_sw = obj.x(i-1,j-1); y_sw = obj.y(i-1,j-1);
                    x_w  = obj.x(i-1,j);   y_w  = obj.y(i-1,j);
                    x_nw = obj.x(i-1,j+1); y_nw = obj.y(i-1,j+1);

                    x_s = obj.x(i,j-1); y_s = obj.y(i,j-1);
                    x_p = obj.x(i,j);   y_p = obj.y(i,j);
                    x_n = obj.x(i,j+1); y_n = obj.y(i,j+1);

                    x_se = obj.x(i+1,j-1); y_se = obj.y(i+1,j-1);
                    x_e  = obj.x(i+1,j);   y_e  = obj.y(i+1,j);
                    x_ne = obj.x(i+1,j+1); y_ne = obj.y(i+1,j+1);

                    x_loc = [x_sw, x_w, x_nw, x_s, x_p, x_n, x_se, x_e, x_ne]';
                    y_loc = [y_sw, y_w, y_nw, y_s, y_p, y_n, y_se, y_e, y_ne]';

                    % Compute Coefficients
                    coefs = get_coefs(x_loc,y_loc);
                    a_sw(k) = coefs(1); a_w(k) = coefs(2); a_nw(k) = coefs(3);
                    a_s(k)  = coefs(4); a_p(k) = coefs(5); a_n(k)  = coefs(6);
                    a_se(k) = coefs(7); a_e(k) = coefs(8); a_ne(k) = coefs(9);
                    b_p(k)  = coefs(10);
                end
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % function t = closest_point(r, curve, t0)
        %
        % Finds the parametric value for the curve that minimizes the distance to 
        % point 'r'
        %
        % Parameters:
        %     r         Point [x,y]
        %     curve     s(t)->[x,y] parametric curve
        %     t0        Initial parameter value guess
        % Returns:
        %     t         Parameter value that minimizes distance from the curve to
        %               the point 'r'
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function t = closest_point(obj, r, curve, t0)
            arguments
                obj
                r     (1,2) double
                curve (1,1) function_handle 
                t0    (1,1) double          = 0
            end

            %Distance function
            dist = @(z) norm(r - curve(z));

            %Loop initialization
            tol = 1.0E-05;
            h = 1.0E-04;
            iter = 0;

            %Move limits
            move_limit = 1/(obj.imax-1);

            %Initial Evaluation Points
            z_e = t0+h;
            z_p = t0;
            z_w = t0-h;
            
            d_p = dist(z_p);
            %Check if early return is possible
            if (d_p<1.0e-06) %Point is on the curve
                t=t0;
                return;
            end
            d_w = dist(z_w);
            d_e = dist(z_e);
            res = abs(d_e - d_w)/(2*h);

            while (res > tol)
                iter = iter+1;
                if (iter > 100)
                    if (obj.verbose)
                        warning('Failed to find new point on curve during iteration of OGridGenerator.mesh(). Commencing with previous point');
                    end
                    t=t0;
                    return
                end
                %Newton's method
                slope = (d_e - d_w)/(2*h);
                curvature = (d_e - 2*d_p + d_w)/(h^2);
                %Locally convex -> Newton's method
                if (curvature > 0)
                    dz = 0.5*(-slope/curvature);
        %                 %Choose new centroid for stencil
                    if (dz>0)
                        dz = min(dz,move_limit); 
                    else
                        dz = max(dz,-move_limit);
                    end
                    d_new = dist(z_p+dz);
                    
                    %Gradient Descent Fallback
                    if (d_new > d_p)
                        dz = -sign(slope)*move_limit;
                        d_new = dist(z_p+dz);
                        sub_iter = 0;
                        while (d_new > d_p)
                            dz = 0.5*dz;
                            %move_limit = 0.9*move_limit;
                            d_new = dist(z_p+dz);
                            sub_iter = sub_iter+1;
                            if (sub_iter>10)
                                break;
                            end
                        end
                    end
                    z_p = z_p + dz;
                    d_p = d_new;
                %Gradient descent subiteration    
                else
                    dz = -sign(slope)*move_limit;
                    d_new = dist(z_p+dz);
                    sub_iter = 0;
                    while (d_new > d_p)
                        dz = 0.5*dz;
                        %move_limit = 0.9*move_limit;
                        d_new = dist(z_p+dz);
                        sub_iter = sub_iter+1;
                        if (sub_iter>10)
                            break;
                        end
                    end
                    z_p = z_p+dz;
                    d_p = d_new;
                end

                %New evaluation points
                z_e = z_p + h;
                z_w = z_p - h;

                %New Sample points
                d_w = dist(z_w);
                d_e = dist(z_e);

                %Calculate residual at new point
                res = abs(d_e -  d_w)/(2*h);
            end
            % Calculate return value
            t = z_p;
        end
    end
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        Helper Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function coefs = get_coefs(x,y)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get 9-point stencil coefficients for an interior discretizatin of the
% Poisson equation. Inputs and outputs are stored as 1-D vectors of 9 
% (10 for output) values in the Cardinal locality labling:
%     (1)  SW
%     (2)  W
%     (3)  NW
%     (4)  S
%     (5)  P
%     (6)  N
%     (7)  SE
%     (8)  E
%     (9)  NE
%     (10) b 
%
% Parameters:
%     x        9x1 vector of x-values from previos iteration
%     y        9x1 vector of y-values from previos iteration
% Returns:
%     coefs    10x1 vector of coefficient values (assuming the equation is
%              solved implicitly, i.e. sign of a_p is not flipped) and b_p
%              value.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    arguments
        x (9,1) double
        y (9,1) double
    end
    
    %Local Variables (for convenience since matlab doesn't do aliasing)
    a_sw = 0.0;
    a_w  = 0.0;
    a_nw = 0.0;
    a_s  = 0.0;
    a_p  = 0.0;
    a_n  = 0.0;
    a_se = 0.0;
    a_e  = 0.0;
    a_ne = 0.0;
    b_p  = 0.0;
    
    %Unpack parameters
    x_sw = x(1);
    x_w  = x(2);
    x_nw = x(3);
    x_s  = x(4);
    x_p  = x(5);
    x_n  = x(6);
    x_se = x(7);
    x_e  = x(8);
    x_ne = x(9);
    
    y_sw = y(1);
    y_w  = y(2);
    y_nw = y(3);
    y_s  = y(4);
    y_p  = y(5);
    y_n  = y(6);
    y_se = y(7);
    y_e  = y(8);
    y_ne = y(9);
    
    % Tablulated values
    x_xi  = 0.5*(x_e - x_w);
    x_eta = 0.5*(x_n - x_s);
    y_xi  = 0.5*(y_e - y_w);
    y_eta = 0.5*(y_n - y_s);
    alpha = x_eta^2 + y_eta^2;
    beta  = x_xi*x_eta + y_xi*y_eta;
    gamma = x_xi^2 + y_xi^2;
    %J     = x_xi*x_eta - y_xi*y_eta; % Needed for wall control (not currently implemented)
   
    % alpha * x_xi_xi
    a_e = a_e + alpha;
    a_p = a_p - 2*alpha;
    a_w = a_w + alpha;
    
    % -2*beta * x_xi_eta
    a_ne = a_ne - 0.5*beta;
    a_nw = a_nw + 0.5*beta;
    a_se = a_se + 0.5*beta;
    a_sw = a_sw - 0.5*beta;
    
    % gammma * x_eta_eta
    a_n = a_n + gamma;
    a_p = a_p - 2*gamma;
    a_s = a_s + gamma;
    
    
    % Aggregate coefficients
    coefs = [a_sw, a_w, a_nw, a_s, a_p, a_n, a_se, a_e, a_ne, b_p]';
end








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function t = uniform_sample(curve,N,t0,b)
%
% Generates a (approximately) uniformally spaced sampling of a curve.
% Parametric curve is not assumed to have a 'natural' parameterization.
% Assumes a closed curve [x,y] = s(t), t0 < t < tf, s(t0)=s(tf) 
% 
% Dev Note: This is pretty crudely implemented. Should be upgraded later 
% after proof-of-concept.
%
% Parameters:
%     curve     s(t)->(x,y) parametric curve
%     N         Number of divisions
%     t0        Initial parameter value
%     tf        Final parameter value
% Returns:
%     ts        (N,1) Parametric point locations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ts = uniform_sample(curve,N,t0,tf)
    arguments
        curve function_handle
        N  (1,1) double  %Actually integer, but Matlab weak typing has weird rules
        t0 (1,1) double = 0
        tf (1,1) double = 1
    end
    ts = zeros(N,1);
    
    %Step size
    m = 50*N; %20 steps between grid points
    dt = (tf-t0)/(m-1);
    
    %Initial position
    t = t0;
    r = curve(t0);
    
    %Get total length of curve
    L = 0;
    for i=2:m
        t=t+dt;
        r_new = curve(t);
        L = L + norm(r_new-r);
        r = r_new;
    end
    
    %Length between gridpoints
    dl = L/N;
    
    %Reset Initial position
    t = t0;
    r = curve(t0);
    
    segment_length = 0;
    index = 2;
    for i=2:m
        t=t+dt;
        r_new = curve(t);
        segment_length = segment_length + norm(r_new-r);
        r = r_new;
        if (segment_length > dl)
            ts(index) = t;
            index = index+1;
            segment_length = 0;
            if (index > N)
                break;
            end
        end
    end 
end

















