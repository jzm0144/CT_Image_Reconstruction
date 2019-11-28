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
    % Padding the Projections
    g(N_) = 0;

    % Calculate DFT of the Padded Projection
    fg = fft(g);

    % Calculate the frequencies for projection and shift them to align with
    % the DFT of projection
    p = fftshift(-1/(2*T):1/(N_*T):(N_/2 -1)/(N_*T));

    % Using magnitude of p to calculate our  rho filter.
    rhoFilter = abs(p)';

    % Applying the filter in frequency domain and calculate the inverse DFT
    % of the result.
    gFilter = ifft(fg .* rhoFilter);
    
    
    % Truncate the filtered projection to get back to the original size of
    % the g projection.
    gFilter((N+1):end) = [];
    gTruncated = gFilter';
    
    % Save the filtered projection for the corresponding angle.
    newProjections(:, index) = gTruncated;

end


% If you want to plot the Rho Filter
% figure(1);
% plot(-1/(2*T):1/(N_*T):(N_/2 -1)/(N_*T), abs(-1/(2*T):1/(N_*T):(N_/2 -1)/(N_*T)));
% title('rhoFilter');


P = 2*floor(size(projections, 1) / (2*(sqrt(2))));%362;
disp(P);

img = zeros(P);


% Creating x and y grid to help with the backward-transformation or
% backprojection
[x, y] = meshgrid((1:P) - P/2);
y = flipud(y);


% Looping over the number of projections or number of angles.
for i=1:length(theta)
    gf1 = newProjections(:,i);
    % Applying the Transformation using x, y and angle theta. This step
    % calculates the actual location of each point of the projection. And
    % the loop running over it sums for all the projections at the
    % corresponding angles to build the image.
    t = round((x*cos(theta(i)*(pi/180))) + y*sin(theta(i)*(pi/180)));
    
    % Fill the output image the backward projection of each angle
    Q = gf1(t+ ceil(N/2));
    img = img+ Q;
    disp(theta(i));
end

% Normalize the output image after all the iterations. This step is
% required because we summed for all the projections to build the image. So
% we divide by the number of thetas.
img = pi * img/length(theta);
figure(2);
imagesc(abs(fliplr(img)));
title('Filtered Backprojection');
colormap(gray);
