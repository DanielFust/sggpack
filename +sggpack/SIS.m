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
% function theta1 = SIS(A,B,C,D,E,F,G,L,theta1,theta2,N,m,n)
%
% Written by: Daniel Fust 12/19/2018
%
% The purpose of this function is to carry out the Strongly Implicit Solver
% (SIS) algorithm of Shong-Leih Lee (1990) [1] for the iterative solution 
% of a banded diagonal matrix created by the solution of elliptic PDEs in 
% 2D with a 9-point stencil. Source algorithm assumes row-major indexing.
%
% A(m*n,1)...........SW coefficient
% B(m*n,1)...........W coefficient
% C(m*n,1)...........NW coefficient
% D(m*n,1)...........S coefficient
% E(m*n,1)...........P coefficient
% F(m*n,1)...........N coefficient
% G(m*n,1)...........SE coefficient
% L(m*n,1)...........E coefficient
% theta1(m*n,1)......NE coefficient
% theta2(m*n,1)......Initial Guess
% N(m*n,1)...........RHS vector
% m..................Number of rows in matrix
% n..................Number of columns in matrix
%
% Possible calling statement for row major ordering:
%     SIS(a_sw,a_w,a_nw,a_s,a_p,a_n,a_ne,a_e,a_ne,x0,b_p,m,n)
%
% Possible calling statement for column major ordering:
%     SIS(a_sw,a_s,a_se,a_w,a_p,a_e,a_nw,a_n,a_ne,x0,b_p,n,m)
%
% [1] Lee, Shong-Leih, "A Strongly Implicit Solver for Two-Dimensional
% Elliptic Differential Equations," Numerical Heat Transfer, Part B
% Fundamentatls, 16:2, 161-178, DOI:10.1080/10407798908944933
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function theta1 = SIS(A,B,C,D,E,F,G,L,theta1,theta2,N,m,n)
    alpha = zeros(m*n,1);
    beta = zeros(m*n,1);
    H = zeros(m*n,1);
    phi = zeros(m*n,1);


    %% LU decomposition
    for k = 1:(m*n)
        %save unmodified coefficients
        aSW = A(k);
        aW = B(k);
        aNW = C(k);
        aS = D(k);
        aP = E(k);
        aN = F(k);
        aSE = G(k);
        aE = L(k);
        aNE = theta1(k);

        if (k==1)
            A(k) = aSW;
            B(k) = aW;
            C(k) = aS;   
            alpha(k) = -aNW;
            beta(k) = -aSE;
            D(k) = aP;
            E(k) = aN/D(k);
            F(k) = aE/D(k);
            G(k) = aNE/D(k);
        elseif ((k-n)<1)
            A(k) = aSW;
            B(k) = aW;
            C(k) = aS;  
            alpha(k) = -aNW;
            beta(k) = C(k)*F(k-1) - aSE;
            D(k) = aP - C(k)*E(k-1);
            E(k) = aN/D(k);
            F(k) = (aE - C(k)*G(k-1))/D(k);
            G(k) = aNE/D(k);
        elseif ((k-n-1)<1)
            A(k) = aSW;
            B(k) = aW;
            C(k) = aS;  
            alpha(k) = B(k)*E(k-n) - aNW;
            beta(k) = C(k)*F(k-1) - aSE;
            D(k) = aP - B(k)*F(k-n) - C(k)*E(k-1);
            E(k) = (aN - B(k)*G(k-n))/D(k);
            F(k) = (aE - C(k)*G(k-1))/D(k);
            G(k) = aNE/D(k);
        else
            A(k) = aSW;
            B(k) = aW - A(k)*E(k-n-1);
            C(k) = aS - A(k)*F(k-n-1);  
            alpha(k) = B(k)*E(k-n) - aNW;
            beta(k) = C(k)*F(k-1) - aSE;
            D(k) = aP - A(k)*G(k-n-1) - B(k)*F(k-n) - C(k)*E(k-1);
            E(k) = (aN - B(k)*G(k-n))/D(k);
            F(k) = (aE - C(k)*G(k-1))/D(k);
            G(k) = aNE/D(k);
        end
    end



    % Initial Guess
    theta1 = theta2;

    %% Iteration
    iter = 0;
    res = 1;
    SOR = 0.65;
    TOL = 1e-6; %tolerance 
    while (res > TOL)
        iter = iter+1;

        %Get new H vector
        for k = 1:(m*n)
            if (k-n+1)<1
                H(k) = N(k) + beta(k)*theta1(k+n-1);
            elseif (k+n-1)>(m*n)
                H(k) = N(k) + alpha(k)*theta1(k-n+1);
            else
                H(k) = N(k) + alpha(k)*theta1(k-n+1) + beta(k)*theta1(k+n-1);
            end
        end

        %save last guess
        theta2 = theta1;

        % Lower diagonal system [L]{phi} = {H} 
        phi(1) = H(1)/D(1);
        for k = 2:(m*n)
            if k<(n+1)
                phi(k) = (H(k) - C(k)*phi(k-1))/D(k);
            elseif k==(n+1)
                phi(k) = (H(k) - C(k)*phi(k-1) - B(k)*phi(k-n))/D(k);
            else
                phi(k) = (H(k) - C(k)*phi(k-1) - B(k)*phi(k-n) - A(k)*phi(k-n-1))...
                          /D(k);
            end
        end

        % Upper diagonal system [U]{theta1} = {phi} 
        theta1(m*n) = phi(m*n);
        for k = (m*n - 1):-1:1
            if k>(m*n-n) 
                theta1(k) = phi(k) - E(k)*theta1(k+1);
            elseif k==(m*n-n)
                theta1(k) = phi(k) - E(k)*theta1(k+1) - F(k)*theta1(k+n);
            else
                theta1(k) = phi(k) - E(k)*theta1(k+1) - F(k)*theta1(k+n) ...
                            - G(k)*theta1(k+n+1);
            end
        end

        %Successive Overrelaxation
        theta1 = theta2 + SOR*(theta1 - theta2);

        %calculate normalized residual
        res = max(abs((theta1-theta2)/(max(theta1)-min(theta1))));
    %     fprintf('iter = %i, res = %f \n',iter,res)
    end

end %function SIS()





