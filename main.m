%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computed Imaging Systems - ELEC- 6810 %
% Programming Project                   %
% Author: Janzaib Masood                %
% Auburn University MRI Research Center %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%---- Sinogram ----%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make the list of angles
thetaMax = 180;
theta = 1:1:thetaMax;

% Load the input image
image = phantom(256);
%image = imread('cameraman.tif');

% Calcuate all the projections using the Radon Transform
projections = myRadon(image, theta);
xp = size(projections,1);
xp = ceil(-xp/2):1:ceil(xp/2);
disp(xp);

% Display the Sinogram
subplot(411);
imagesc(theta, xp, projections);
title('Sinogram')
colormap(gray);



%%%%%%%%%%%%%%%%%%%%---- Filtered Backprojection ----%%%%%%%%%%%%%%%%%%%%%
N = size(projections, 1); % Number of elements in a projection

N_ = 2*(2^nextpow2(N)); % Number of elements in a zero padded projection

disp(N_);
T = 1; 

newProjections = zeros(size(projections));

% Loop over the angles from 1 to thetaMax
for angle=1:1:thetaMax
    
    g = projections(:,angle);
    % Padd the Projection
    g(N_) = 0;

    % Calculate DFT of the Padded Projection
    fg = fft(g);

    % Calculate the frequencies for projection and shift them to align with
    % the DFT of projection
    p = fftshift(-1/(2*T):1/(N_*T):(N_/2 -1)/(N_*T));

    % Using magnitude of p as our filter.
    rhoFilter = abs(p)';

    % Applying the filter in frequency domain and calculate the inverse DFT
    % of the result.
    gFilter = ifft(fg .* rhoFilter);
    
    
    % Truncate the filtered projection
    gFilter((N+1):end) = [];
    gTruncated = gFilter';
    
    % Save the filtered projection for the angle.
    newProjections(:, angle) = gTruncated;

end

subplot(412);
imagesc(theta, xp, newProjections);
title('filtered Sinogram');
colormap(gray);


subplot(413);
plot(-1/(2*T):1/(N_*T):(N_/2 -1)/(N_*T), abs(-1/(2*T):1/(N_*T):(N_/2 -1)/(N_*T)));
title('rhoFilter');


P = size(image,1);
img = zeros(P);


% Creating x and y grid to help with the backward-transformation or
% backprojection
[x, y] = meshgrid((1:P) - P/2);
y = flipud(y);


% Looping over the number of projections or number of angles.
for i=1:length(theta)
    gf1 = newProjections(:,i);
    % Applying the Transformation using x, y and angle theta. This step
    % calculates the 
    t = round((x*cos(theta(i)*(pi/180))) + y*sin(theta(i)*(pi/180)));
    
    % Fill the output image the backward projection of each angle
    img = img + gf1(t + ceil(N/2));
end

% Normalize the output image after all the iterations.
img = pi * img/length(theta);



subplot(414);
imagesc(img);
title('Reconstructed Image');
colormap(gray);