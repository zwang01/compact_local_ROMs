function [psi, xn, yn, err] = mac_streamfunction(par, u, v)
% MAC_STREAMFUNCTION  Discrete stream function of a 2D MAC velocity field on
% the cell corners, used for the streamline panels of the figures.
%
%   [psi, xn, yn, err] = mac_streamfunction(par, u, v)
%
% u is (Nx+1) x Ny on the x-faces and v is Nx x (Ny+1) on the y-faces. The
% stream function (dpsi/dy = u, dpsi/dx = -v) integrates exactly from corner to
% corner,
%   psi(i,j+1) - psi(i,j) = +h*u(i,j),   psi(i+1,j) - psi(i,j) = -h*v(i,j),
% because the discrete divergence vanishes in every cell (the continuity
% equation of the Brinkman model is the unweighted div u = 0), so the result
% is path independent. psi is (Nx+1) x (Ny+1) on the corners xn, yn with
% psi = 0 at the first corner; err is the largest cell-loop mismatch, i.e. h
% times the largest discrete divergence (a check that the field is solenoidal).
Nx = par.Nx; Ny = par.Ny; h = par.h;
xn = par.x0 + (0:Nx)*h;
yn = par.y0 + (0:Ny)*h;
err = max(abs( u(2:Nx+1,:) - u(1:Nx,:) + v(:,2:Ny+1) - v(:,1:Ny) ), [], 'all');
psi = zeros(Nx+1, Ny+1);
psi(1, 2:Ny+1) = cumsum(h*u(1,:));              % up the left edge
psi(2:Nx+1, :) = -h*cumsum(v, 1);               % then rightwards along each row
psi(2:Nx+1, :) = psi(2:Nx+1, :) + psi(1, :);
end
