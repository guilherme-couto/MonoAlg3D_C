function mesh = buildMesh(domain_config)
% This function simply builds a basic struct that contains the mesh
% information for use inside other functions
%
% Usage:
%     mesh = buildMesh(domain_config)
%
% INPUTS:
%
%     domain_config: A struct containing the information:
%       - Lx : Width of the domain (cm)
%       - Ly : Height of the domain (cm)
%       - Lz : Depth of the domain (cm) [Optional, default = 0]
%       - dx : Spatial resolution (cm)
%
% OUPUTS:
%
%   mesh: A struct with fields "points" (point locations), "Nx", "Ny", "Nz" (number of
%         points in each direction), "dx" (spatial resolution), and "is3D" (boolean).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extract parameters
Lx = domain_config.Lx;
Ly = domain_config.Ly;
dx = domain_config.dx;

% Check for optional Lz
if isfield(domain_config, 'Lz')
    Lz = domain_config.Lz;
else
    Lz = 0; % Default to 0 (2D)
end

% Compute number of points in each direction
Nx = ceil(Lx/dx);
Ny = ceil(Ly/dx);

% Create a set of points that fall in centres of pixels
if Lz == 0
    % 2D case
    xv = linspace(dx/2, dx*(Nx - 0.5), Nx);
    yv = linspace(dx/2, dx*(Ny - 0.5), Ny);
    [X,Y] = meshgrid(xv,yv);
    points = [X(:),Y(:)];
    mesh.Nz = 1;
    mesh.is3D = false;
else
    % 3D case
    Nz = ceil(Lz/dx);
    xv = linspace(dx/2, dx*(Nx - 0.5), Nx);
    yv = linspace(dx/2, dx*(Ny - 0.5), Ny);
    zv = linspace(dx/2, dx*(Nz - 0.5), Nz);
    [X,Y,Z] = meshgrid(xv,yv,zv);
    points = [X(:),Y(:),Z(:)];
    mesh.Nz = Nz;
    mesh.is3D = true;
end

% Store in mesh struct
mesh.points = points;
mesh.Nx = Nx;
mesh.Ny = Ny;
mesh.dx = dx;
end
