#!/usr/bin/octave -qf

% ------------------------ START OF MESH PARAMETER REGION -------------------- %

% Foil geometry
c = 1;                % Geometric chord length
alpha = 0.1745;       % Angle of attack (in radians)
NACA = [4 4 1 2];     % NACA 4-digit designation as a row vector;

% Mesh dimensions
scale = 1;            % Scaling factor
H = 8;                % *Half* height of channel
W = 0.5;              % *Half* depth of foil (y-direction)
D = 16;               % Length of downstream section

% Mesh resolution parameters
Ni = 400;             % Number of interpolation points along the foil
Nx = 250;             % Number of mesh cells along the foil
ND = 150;             % Number of cells in the downstream direction
NT = 100;             % Number of cells the transverse direction
NW = 1;               % Number of cells in the y-direction (along the foil axis)

% Expansion rates
ExpT = 500;           % Expansion rate in transverse direction
ExpD = 100;           % Expansion rate in the downstream direction
ExpArc = 50;          % Expansion rate along the inlet arc

% ------------------------- END OF MESH PARAMETER REGION --------------------- %


% ---------------------------------- LICENCE  -------------------------------- %
%                                                                              %
%     Copyrighted 2011, 2012 by Håkon Strandenes, hakostra@stud.ntnu.no        %
%                                                                              % 
%     This program is free software: you can redistribute it and/or modify     %
%     it under the terms of the GNU General Public License as published by     %
%     the Free Software Foundation, either version 3 of the License, or        %
%     (at your option) any later version.                                      %
%                                                                              %
%     This program is distributed in the hope that it will be useful,          %
%     but WITHOUT ANY WARRANTY; without even the implied warranty of           %
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            %
%     GNU General Public License for more details.                             %
%                                                                              %
%     You should have received a copy of the GNU General Public License        %
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.    %
% ---------------------------------------------------------------------------- %


% Create a vector with x-coordinates, camber and thickness
beta=linspace(0,pi,Ni);
x = c*(0.5*(1-cos(beta)));
z_c = zeros(size(x));
z_t = zeros(size(x));
theta = zeros(size(x));


% Values of m, p and t
m = NACA(1)/100;
p = NACA(2)/10;
t = (NACA(3)*10 + NACA(4))/100;


% Calculate thickness
% The upper expression will give the airfoil a finite thickness at the trailing
% edge, witch might cause trouble. The lower expression is correxted to give 
% zero thickness at the trailing edge, but the foil is strictly speaking no
% longer a proper NACA airfoil.
%
% See http://turbmodels.larc.nasa.gov/naca4412sep_val.html
%     http://en.wikipedia.org/wiki/NACA_airfoil

%y_t = (t*c/0.2) * (0.2969.*(x/c).^0.5 - 0.1260.*(x/c) - 0.3516.*(x/c).^2 + 0.2843.*(x/c).^3 - 0.1015.*(x/c).^4);
z_t = (t*c/0.2) * (0.2969.*(x/c).^0.5 - 0.1260.*(x/c) - 0.3516.*(x/c).^2 + 0.2843.*(x/c).^3 - 0.1036.*(x/c).^4);


% Calculate camber
if (p > 0)
  % Calculate camber
  z_c = z_c + (m.*x/p^2) .* (2*p - x/c) .* (x < p*c);
  z_c = z_c + (m.*(c-x)/(1-p)^2) .* (1 + x/c - 2*p) .* (x >= p*c);


  % Calculate theta-value
  theta = theta + atan( (m/p^2) * (2*p - 2*x/c) ) .* (x < p*c);
  theta = theta + atan( (m/(1-p)^2) * (-2*x/c + 2*p) ) .* (x >= p*c);
end


% Calculate coordinates of upper surface
Xu = x - z_t.*sin(theta);
Zu = z_c + z_t.*cos(theta);


% Calculate coordinates of lower surface
Xl = x + z_t.*sin(theta);
Zl = z_c - z_t.*cos(theta);


% Rotate foil to reach specified angle of attack
upper = [cos(alpha), sin(alpha); -sin(alpha), cos(alpha)] * [Xu ; Zu];
lower = [cos(alpha), sin(alpha); -sin(alpha), cos(alpha)] * [Xl ; Zl];

Xu = upper(1,:)';
Zu = upper(2,:)';
Xl = lower(1,:)';
Zl = lower(2,:)';


if (p > 0)
  % Find index i of max. camber
  C_max_idx = find(z_c == max(z_c));
else
  % Otherwise use location of max. thickness
  C_max_idx = find(z_t == max(z_t));
end


% Move point of mesh "nose"
NoseX = (-H+Xu(C_max_idx))*cos(alpha);
NoseZ = -(-H+Xu(C_max_idx))*sin(alpha);


% Calculate the location of the vertices on the positive y-axis and put them in a matrix
vertices = zeros(12, 3);

vertices(1,:)  = [NoseX,            W,            NoseZ];
vertices(2,:)  = [Xu(C_max_idx),    W,                H];
vertices(3,:)  = [Xu(Ni),           W,                H];
vertices(4,:)  = [D,                W,                H];
vertices(5,:)  = [0,                W,                0];
vertices(6,:)  = [Xu(C_max_idx),    W,    Zu(C_max_idx)];
vertices(7,:)  = [Xl(C_max_idx),    W,    Zl(C_max_idx)];
vertices(8,:)  = [Xu(Ni),           W,           Zu(Ni)];
vertices(9,:)  = [D,                W,           Zu(Ni)];
vertices(10,:) = [Xl(C_max_idx),    W,               -H];
vertices(11,:) = [Xu(Ni),           W,               -H];
vertices(12,:) = [D,                W,               -H];

% Create vertices for other side (negative y-axis)
vertices = [vertices; vertices(:,1), -vertices(:,2), vertices(:,3)]; 


% Edge 4-5 and 16-17
pts1 = [Xu(2:C_max_idx-1), W*ones(size(Xu(2:C_max_idx-1))), Zu(2:C_max_idx-1)]; 
pts5 = [pts1(:,1), -pts1(:,2), pts1(:,3)];

% Edge 5-7 and 17-19
pts2 = [Xu(C_max_idx+1:Ni-1), W*ones(size(Xu(C_max_idx+1:Ni-1))), Zu(C_max_idx+1:Ni-1)]; 
pts6 = [pts2(:,1), -pts2(:,2), pts2(:,3)];

% Edge 4-6 and 16-18
pts3 = [Xl(2:C_max_idx-1), W*ones(size(Xl(2:C_max_idx-1))), Zl(2:C_max_idx-1)]; 
pts7 = [pts3(:,1), -pts3(:,2), pts3(:,3)];

% Edge 6-7 and 18-19
pts4 = [Xl(C_max_idx+1:Ni-1), W*ones(size(Xl(C_max_idx+1:Ni-1))), Zl(C_max_idx+1:Ni-1)]; 
pts8 = [pts4(:,1), -pts4(:,2), pts4(:,3)];

% Edge 0-1 and 12-13
pts9 = [-H*cos(pi/4)+Xu(C_max_idx), W, H*sin(pi/4)];
pts11 = [pts9(:,1), -pts9(:,2), pts9(:,3)];

% Edge 0-9 and 12-21
pts10 = [-H*cos(pi/4)+Xu(C_max_idx), W, -H*sin(pi/4)];
pts12 = [pts10(:,1), -pts10(:,2), pts10(:,3)];


% Calculate number of mesh points along 4-5 and 4-6
%Nleading = (C_max_idx/Ni)*Nx;
Nleading = (x(C_max_idx)/c)*Nx;

% Calculate number of mesh points along 5-7 and 6-7
Ntrailing = Nx-Nleading;


% Open file
fo = fopen('constant/polyMesh/blockMeshDict', 'w');


% Write file
fprintf(fo, '/*--------------------------------*- C++ -*----------------------------------*\\ \n');
fprintf(fo, '| =========                 |                                                 | \n');
fprintf(fo, '| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           | \n');
fprintf(fo, '|  \\\\    /   O peration     | Version:  2.1.0                                 | \n');
fprintf(fo, '|   \\\\  /    A nd           | Web:      www.OpenFOAM.com                      | \n');
fprintf(fo, '|    \\\\/     M anipulation  |                                                 | \n');
fprintf(fo, '\\*---------------------------------------------------------------------------*/ \n');
fprintf(fo, 'FoamFile                                                                        \n');
fprintf(fo, '{                                                                               \n');
fprintf(fo, '    version     2.0;                                                            \n');
fprintf(fo, '    format      ascii;                                                          \n');
fprintf(fo, '    class       dictionary;                                                     \n');
fprintf(fo, '    object      blockMeshDict;                                                  \n');
fprintf(fo, '}                                                                               \n');
fprintf(fo, '// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * // \n');
fprintf(fo, '\n');
fprintf(fo, 'convertToMeters %f; \n', scale);
fprintf(fo, '\n');
fprintf(fo, 'vertices \n');       
fprintf(fo, '( \n');
fprintf(fo, '    (%f %f %f)\n', vertices');
fprintf(fo, '); \n');
fprintf(fo, '\n');
fprintf(fo, 'blocks \n');
fprintf(fo, '( \n');
fprintf(fo, '    hex (4 5 1 0 16 17 13 12)     (%i %i %i) edgeGrading (1 %f %f 1 %f %f %f %f 1 1 1 1) \n', Nleading, NT, NW, 1/ExpArc, 1/ExpArc, ExpT, ExpT, ExpT, ExpT);
fprintf(fo, '    hex (5 7 2 1 17 19 14 13)     (%i %i %i) simpleGrading (1 %f 1) \n', Ntrailing, NT, NW, ExpT);
fprintf(fo, '    hex (7 8 3 2 19 20 15 14)     (%i %i %i) simpleGrading (%f %f 1) \n', ND, NT, NW, ExpD, ExpT);
fprintf(fo, '    hex (16 18 21 12 4 6 9 0)     (%i %i %i) edgeGrading (1 %f %f 1 %f %f %f %f 1 1 1 1) \n', Nleading, NT, NW, 1/ExpArc, 1/ExpArc, ExpT, ExpT, ExpT, ExpT);
fprintf(fo, '    hex (18 19 22 21 6 7 10 9)    (%i %i %i) simpleGrading (1 %f 1) \n', Ntrailing, NT, NW, ExpT);
fprintf(fo, '    hex (19 20 23 22 7 8 11 10)   (%i %i %i) simpleGrading (%f %f 1) \n', ND, NT, NW, ExpD, ExpT);

fprintf(fo, '); \n');
fprintf(fo, '\n');
fprintf(fo, 'edges \n');
fprintf(fo, '( \n');

fprintf(fo, '    spline 4 5 \n');
fprintf(fo, '        ( \n');
fprintf(fo, '            (%f %f %f) \n', pts1');
fprintf(fo, '        ) \n');

fprintf(fo, '    spline 5 7 \n');
fprintf(fo, '        ( \n');
fprintf(fo, '            (%f %f %f)\n', pts2');
fprintf(fo, '        ) \n');

fprintf(fo, '    spline 4 6 \n');
fprintf(fo, '        ( \n');
fprintf(fo, '            (%f %f %f)\n', pts3');
fprintf(fo, '        ) \n');

fprintf(fo, '    spline 6 7 \n');
fprintf(fo, '        ( \n');
fprintf(fo, '            (%f %f %f)\n', pts4');
fprintf(fo, '        ) \n');

fprintf(fo, '    spline 16 17 \n');
fprintf(fo, '        ( \n');
fprintf(fo, '            (%f %f %f)\n', pts5');
fprintf(fo, '        ) \n');

fprintf(fo, '    spline 17 19 \n');
fprintf(fo, '        ( \n');
fprintf(fo, '            (%f %f %f)\n', pts6');
fprintf(fo, '        ) \n');

fprintf(fo, '    spline 16 18 \n');
fprintf(fo, '        ( \n');
fprintf(fo, '            (%f %f %f)\n', pts7');
fprintf(fo, '        ) \n');

fprintf(fo, '    spline 18 19 \n');
fprintf(fo, '        ( \n');
fprintf(fo, '            (%f %f %f)\n', pts8');
fprintf(fo, '        ) \n');

fprintf(fo, '    arc 0 1 (%f %f %f) \n', pts9');
fprintf(fo, '    arc 0 9 (%f %f %f) \n', pts10');
fprintf(fo, '    arc 12 13 (%f %f %f) \n', pts11');
fprintf(fo, '    arc 12 21 (%f %f %f) \n', pts12');

fprintf(fo, '); \n');
fprintf(fo, '\n');
fprintf(fo, 'boundary \n');
fprintf(fo, '( \n');

fprintf(fo, '    inlet \n');
fprintf(fo, '    { \n');
fprintf(fo, '        type patch; \n');
fprintf(fo, '        faces \n');
fprintf(fo, '        ( \n');
fprintf(fo, '            (1 0 12 13) \n');
fprintf(fo, '            (0 9 21 12) \n');
fprintf(fo, '        ); \n');
fprintf(fo, '    } \n');
fprintf(fo, '\n');

fprintf(fo, '    outlet \n');
fprintf(fo, '    { \n');
fprintf(fo, '        type patch; \n');
fprintf(fo, '        faces \n');
fprintf(fo, '        ( \n');
fprintf(fo, '            (11 8 20 23) \n');
fprintf(fo, '            (8 3 15 20) \n');
fprintf(fo, '        ); \n');
fprintf(fo, '    } \n');
fprintf(fo, '\n');

fprintf(fo, '    topAndBottom \n');
fprintf(fo, '    { \n');
fprintf(fo, '        type patch; \n');
fprintf(fo, '        faces \n');
fprintf(fo, '        ( \n');
fprintf(fo, '            (3 2 14 15) \n');
fprintf(fo, '            (2 1 13 14) \n');
fprintf(fo, '            (9 10 22 21) \n');
fprintf(fo, '            (10 11 23 22) \n');
fprintf(fo, '        ); \n');
fprintf(fo, '    } \n');
fprintf(fo, '\n');

fprintf(fo, '    airfoil \n');
fprintf(fo, '    { \n');
fprintf(fo, '        type wall; \n');
fprintf(fo, '        faces \n');
fprintf(fo, '        ( \n');
fprintf(fo, '            (5 4 16 17) \n');
fprintf(fo, '            (7 5 17 19) \n');
fprintf(fo, '            (4 6 18 16) \n');
fprintf(fo, '            (6 7 19 18) \n');
fprintf(fo, '        ); \n');
fprintf(fo, '    } \n');
fprintf(fo, '); \n');
fprintf(fo, ' \n');
fprintf(fo, 'mergePatchPairs \n');
fprintf(fo, '( \n');
fprintf(fo, '); \n');
fprintf(fo, ' \n');
fprintf(fo, '// ************************************************************************* // \n');

% Close file
fclose(fo);
