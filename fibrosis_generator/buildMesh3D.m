function mesh = buildMesh3D(Nx, Ny, Nz, dx)
% Builds a 3D mesh with voxel centers and spacing dx
%
% Inputs:
%   Nx, Ny, Nz - number of volumes in x, y, z directions
%   dx         - voxel spacing
%
% Output:
%   mesh struct with fields:
%       .points   - point coordinates of the mesh
%       .N_points - [Nx, Ny, Nz]
%       .dx       - voxel side length

    % Coordinates of voxel centers
    xv = linspace(dx/2, dx*(Nx - 0.5), Nx);
    yv = linspace(dx/2, dx*(Ny - 0.5), Ny);
    zv = linspace(dx/2, dx*(Nz - 0.5), Nz);

    % 3D grid of centers
    [Y,X,Z] = meshgrid(yv, xv, zv);
    mesh.points = [X(:) Y(:) Z(:)];

    mesh.Nx = Nx;
    mesh.Ny = Ny;
    mesh.Nz = Nz;

    mesh.dx = dx;
    mesh.dy = dx;
    mesh.dz = dx;
end
