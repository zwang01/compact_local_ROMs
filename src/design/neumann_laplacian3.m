function L = neumann_laplacian3(Nx, Ny, Nz, h)
% Cell-centered 7-point Laplacian on an Nx x Ny x Nz grid with homogeneous
% Neumann boundary conditions (3D Allen-Cahn design update).
T = @(n) spdiags([ones(n,1) -2*ones(n,1) ones(n,1)], [-1 0 1], n, n) ...
    + sparse([1 n],[1 n],[1 1],n,n);
Ix = speye(Nx); Iy = speye(Ny); Iz = speye(Nz);
L = (kron(Iz, kron(Iy, T(Nx))) + kron(Iz, kron(T(Ny), Ix)) ...
   + kron(T(Nz), kron(Iy, Ix)))/h^2;
end
