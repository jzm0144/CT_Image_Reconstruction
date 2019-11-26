%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computed Imaging Systems - ELEC- 6810 %
% Programming Project                   %
% Author: Janzaib Masood                %
% Auburn University MRI Research Center %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%---- Sinogram ----%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Load the input image
load project_1_new_data.mat;


projections = kspace;
numProjections = size(projections,2);
%%%%%%%%%%%%%%%%%%%%---- Filtered Backprojection ----%%%%%%%%%%%%%%%%%%%%%
N = size(projections, 1); % Number of elements in a projection

N_ = 2*(2^nextpow2(N)); % Number of elements in a zero padded projection

disp(N_);
T = 180/numProjections; 


newProjections = zeros(size(projections));

% Loop over the angles from 1 to thetaMax
for index=1:1:numProjections
    angle = theta(index);
    
    G = projections(:, index);
    g = abs(fftshift(ifft(G)));
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
    newProjections(:, index) = gTruncated;

end


subplot(223);
plot(-1/(2*T):1/(N_*T):(N_/2 -1)/(N_*T), abs(-1/(2*T):1/(N_*T):(N_/2 -1)/(N_*T)));
title('rhoFilter');


P = 362;
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
    
    Q = gf1(t+ ceil(N/2));
    img = img+ Q;
    disp(theta(i));
end

% Normalize the output image after all the iterations.
img = pi * img/length(theta);

subplot(224);

imagesc(abs(fliplr(img)));
title('Reconstructed Image');
colormap(gray);
