%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computed Imaging Systems - ELEC- 6810 %
% Programming Project                   %
% Author: Janzaib Masood                %
% Auburn University MRI Research Center %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%---- Sinogram ----%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load the input kspace projections and angles
% Radon Transform is not required because we already have the projections.
load project_1_new_data.mat;
projections = kspace;
numProjections = size(projections,2);


%%%%%%%%%%%%%%%%%%%%---- Filtered Backprojection ----%%%%%%%%%%%%%%%%%%%%%
N = size(projections, 1); % Number of elements in a projection
N_ = 2*(2^nextpow2(N)); % Number of elements in a zero padded projection

disp(N_);

T = 180/(numProjections); 

newProjections = zeros(size(projections));

% Loop over the angles from 1 to thetaMax
for index=1:1:numProjections
    angle = theta(index);
    
    % Converting the Kspace Projections to Spatial Domain with an inverse
    % fourier transform.
    G = projections(:, index);
    g = abs(fftshift(ifft(G)));
    newProjections(:, index) = g;

end

P = 362;
img = zeros(P);

[gdat] = gridkb(kspacelocations, spiraldata, dcf, 256, 1.5, 2);
