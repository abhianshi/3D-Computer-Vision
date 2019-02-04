function [result] = Assignment3(directoryName)
    % Load images.
    Scene = imageDatastore(directoryName);

    % Display images to be stitched
    %$montage(Scene.Files)
    
    % Read the first image from the image set.
    I = readimage(Scene, 1);
    
    % Change to grayScale if image is RGB
    if size(I,3)==3
        grayImg = rgb2gray(I);
    end
    
    % Detect and Extract Features for start Image
    points = detectSURFFeatures(grayImg);
    [features,points] = extractFeatures(grayImg,points);
    
    % Total Number of images in this directory
    numImages = numel(Scene.Files);
   
    % Initialize variable to hold image sizes.
    imageSize = zeros(numImages,2);
    tforms(numImages) = projective2d(eye(3));
    
    % Point Correspondence between the adjacent image pairs
    for n = 2:numImages

        % Store points and features for the previous image.
        pointsPrevious = points;
        featuresPrevious = features;

        % Read I(n).
        I = readimage(Scene, n);

        % Convert image to grayscale.
        if size(I,3)==3
            grayImg = rgb2gray(I);
        end
        
        % Save image size.
        imageSize(n,:) = size(grayImg);
    
        % Detect and extract SURF features for the image.
        points = detectSURFFeatures(grayImg);
        [features, points] = extractFeatures(grayImg, points);

        % Find correspondences between the adjacent images.
        indexPairs = matchFeatures(features, featuresPrevious, 'Unique', true);

        matchedPoints = points(indexPairs(:,1), :) ;
        matchedPointsPrev = pointsPrevious(indexPairs(:,2), :);
        
        % Use Ransac to get the Homography
        [H, inliers] = ransacfithomography((matchedPoints.Location)', (matchedPointsPrev.Location)', 0.005);
        H = double(H);
        tforms(n) = projective2d(H');

        % Compute T(n) * T(n-1) * ... * T(1)
        tforms(n).T = tforms(n).T * tforms(n-1).T;
    end
    
    % Compute the output limits  for each transform
    for i = 1:numel(tforms)
        [xlim(i,:), ylim(i,:)] = outputLimits(tforms(i), [1 imageSize(i,2)], [1 imageSize(i,1)]);
    end

    avgXLim = mean(xlim, 2);

    [~, idx] = sort(avgXLim);

    centerIdx = floor((numel(tforms)+1)/2);

    centerImageIdx = idx(centerIdx);

    Tinv = invert(tforms(centerImageIdx));

    for i = 1:numel(tforms)
        tforms(i).T = tforms(i).T * Tinv.T;
    end

    for i = 1:numel(tforms)
        [xlim(i,:), ylim(i,:)] = outputLimits(tforms(i), [1 imageSize(i,2)], [1 imageSize(i,1)]);
    end
    
    maxImageSize = max(imageSize);

    % Find the minimum and maximum output limits
    xMin = min([1; xlim(:)]);
    xMax = max([maxImageSize(2); xlim(:)]);

    yMin = min([1; ylim(:)]);
    yMax = max([maxImageSize(1); ylim(:)]);

    % Width and height of panorama.
    width  = round(xMax - xMin);
    height = round(yMax - yMin);

    % Initialize the "empty" panorama.
    panorama = zeros([height width 3], 'like', I);

    % Combines two images, overlays one image over another
    blender = vision.AlphaBlender('Operation', 'Binary mask', ...
        'MaskSource', 'Input port');
       
    % Create a 2-D spatial reference object defining the size of the panorama.
    xLimits = [xMin xMax];
    yLimits = [yMin yMax];
    panoramaView = imref2d([height width], xLimits, yLimits);

    % Create the panorama.
    for i = 1:numImages

        I = readimage(Scene, i);

        % Transform I into the panorama.
        warpedImage = imwarp(I, tforms(i), 'OutputView', panoramaView, 'interp', 'bilinear');
        
        % Generate a binary mask.
        mask = imwarp(true(size(I,1),size(I,2)), tforms(i), 'OutputView', panoramaView, 'interp', 'bilinear');

        % Overlay the warpedImage onto the panorama.
        panorama = step(blender, panorama, warpedImage, mask);
    end

    figure
    imshow(panorama)
    
end