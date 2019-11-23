function projection = myRadon(image, theta)
    %Inputs
        % theta: a vector of angles or a single angle.
        % image: a 2d image
    %Output
        % projection: matrix where every coloum contain the projection of image
        %             at a given angle

    [rowLength, colLength]= size(image);
    
    disp('Input Image Shape');size(image)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%---- Padding ----%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Padding the Image to avoid losing information on the edges when we rotate.
    diagLength = sqrt(rowLength^2 + colLength^2);
    % Let us pad the image so that its width and height is a little wider than
    % the original width of the image.
    rowPad = ceil(diagLength - rowLength) + 3;
    colPad = ceil(diagLength - colLength) + 3;
    % The new frame with the padded dimensions
    newImage = zeros(rowLength + rowPad, colLength + colPad);
    % Put the image in the center of the frame
    newImage(ceil(rowPad/2):(ceil(rowPad/2)+rowLength-1),...
             ceil(colPad/2):(ceil(rowPad/2)+rowLength-1)) = image;
    image = newImage;
    %%%%%%%%%%%%%%%---- Getting the Projections at all angles ----%%%%%%%%%%%%%

    x = linspace(-1, 1, size(image, 1));
    [X, Y] = meshgrid(x, x);
    %disp('size of meshgrid');
    %disp(size(X));
    %disp(size(Y));


    numAngles = length(theta);
    projection = zeros(size(image,2), numAngles);

    for index=1:numAngles
        angle = (90-theta(index))*pi/180; cos(angle)*X -sin(angle)*Y;

        xVals = cos(angle)*X - sin(angle)*Y;
        yVals = sin(angle)*X + cos(angle)*Y;

        % Interpolate using interp2 method from matlab
        % We do interpolation here to fill the indices
        % where no data is sitting.
        interpImg = interp2(X,Y,image, xVals, yVals);
        % Replace NaN values with zero.
        interpImg(isnan(interpImg)) = 0;
        % Take the sum of the rotated images(with interpolations)
        % in the the axis of the of the detector. At every run of
        % the loop, this step fills the projection matrix for an
        % angle.
        projection(:, index) = (sum(interpImg))'; %'
    end

    disp('Sinogram Shape');size(projection)

end
