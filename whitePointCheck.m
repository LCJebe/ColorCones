clear;
close all;

%% create solid color image
img_size = 512;

color = [1, 1, 1];
img = ones(img_size, img_size, 3);
for i = 1:3
    img(:, :, i) = img(:, :, i)*color(i);
end

%% create scene and stuff
fov = 4;
[s, cMosaic] = createSceneFromImage(img, fov);

%% reconstruct image and get white balance
rec1 = reconstructImage(cMosaic);

% get the white point (asusming whole image is white)
whiteRGB = reshape(rec1, [], 3);
avgWhiteRGB = mean(whiteRGB, 1);

% white balance
rec1_bal = rec1;
for i = 1:3
    rec1_bal(:, :, i) = rec1(:, :, i) / avgWhiteRGB(i);
end
rec1_bal(rec1_bal>1) = 1;
rec1_bal(rec1_bal<0) = 0;

figure();
imshow(rec1_bal);

%% reconstruction function
function rgb = reconstructImage(cMosaic)
    disp('Reconstructing Image from Cone Absorptions...')
    
    search_radii = [2, 3, 5]; % in pixels
    [im] = imageFromConeAbs(cMosaic, search_radii);
    imd = im2double(im);
    %imd = (imd - min(imd(:))) / (max(imd(:)) - min(imd(:)));
    imd = imd / max(imd(:));

    rgb = lms2srgb(imd);
    
    figure();
    imshow(rgb);
    title('Reconstructed Image');
end


%% create scenes and absorptions
function [s, cMosaic] = createSceneFromImage(img, fov)
    % create scene from file
    [s, ~] = sceneFromFile(img, 'rgb');
    s = sceneSet(s, 'fov', fov);

    ieAddObject(s);
    sceneWindow;

    %% Build oi (retinal image) 
    oi = oiCreate;
    oi = oiCompute(oi, s);

    ieAddObject(oi);

    %% Build a default cone mosaic and compute isomerizatoins

    % Create the coneMosaic object
    cMosaic = coneMosaic;

    % Set size to show a fraction of the scene. Speeds things up.
    val = 1; % 1 means full scene
    cMosaic.setSizeToFOV(val * sceneGet(s, 'fov'));

    %disp(sceneGet(s,'fov'));

    %% Set integration time 
    cMosaic.integrationTime = 0.05;

    %% Generate a sequence of 100 eye posistions.
    % cMosaic.emGenSequence(100);

    %% Compute isomerizations for each eye position.
    cMosaic.compute(oi);

    %% Bring up a window so that we can look at things.
    %
    % Using the pull down in the window, you can look at
    % the mosaic, the isomerizations for one fixation, or
    % the movie of fixations.
    cMosaic.window;
end
