%% Initialize and clear
addpath(genpath('/Users/Joe1/isetbio/'));
addpath(genpath('/Users/Joe1/Documents/MATLAB/'));
clc; clear; 
close all;
ieInit;

%% create our own scene form a colorful spoke target

disp('Setting up Scene...');
% create spoke
col1 = [1, 0, 0];
col2 = [0, 0, 1];
spoke_size = 512;
spoke = createSpoke(col1, col2, 25, spoke_size);


%% Modify original scene

% Option 1: Blur spoke 
sigma = 3;
% spoke_blur = imgaussfilt(spoke, sigma);

% Option 2: Don't change the spoke (used for checking noise) 
spoke_blur = spoke;

% Option 3: Switch colors 
%
spoke_blur = createSpoke(col2, col1, 10, spoke_size);


spoke = im2double(spoke);
spoke_blur = im2double(spoke_blur);


%% First do S-CIELAB calculations between scene and modified scene

%% define spacial calibration (samples per degree of visual angle)
dpi = 72; % dots per inch
dist = 18; % viewing distances in inches

sampPerDeg = round(dpi / ((180/pi)*atan(1/dist)));


%% define display / color calibration
% load given display spectral power density and the cone absorption curves
load displaySPD
load SmithPokornyCones 

% get transformation from rgb to lms
rgb2lms = cones'* displaySPD; 

% define or load display gamma table
load displayGamma

% define white point (neccessary for S-CIELAB)
rgbWhite = [1 1 1]';
whitepoint = rgb2lms* rgbWhite;


%% set up image data for S-CIELAB
spokeXYZ = rgb2xyz(spoke);
spoke_blurXYZ = rgb2xyz(spoke_blur);
%imgRGB = dac2rgb(img01, gammaTable);
%imgLMS = changeColorSpace(imgRGB, rgb2lms);

% define image input format for S-CIELAB
%imageformat = 'lms';


%% Perform S-CIELAB calculation
%filtOpp = scielab(sampPerDeg, imgLMS, imageformat); 
%filtXYZ = changeColorSpace(filtOpp, cmatrix('opp2xyz')); 

% get S-CIELAB delta E as difference to the white point of the image
whitepointRGB = [1, 1, 1] * 0.5;
whitepointXYZ = rgb2xyz(whitepointRGB);

errorImage = scielab(sampPerDeg, spokeXYZ, spoke_blurXYZ, whitepointXYZ, 'xyz');


%% analyze S-CIELAB image (it's not really an error image, just a delta E
    % image with respect to a reference white image

%% show images
fig_title = sprintf('Original & blurred scenes, sigma = %.1f', sigma);
figure('name', fig_title)
subplot(1, 3, 1);
imshow(spoke);
subplot(1, 3, 2);
imshow(spoke_blur);
subplot(1, 3, 3);
imshow(errorImage, [0, max(errorImage(:))]);
colorbar;


%% Now calculate cone responses from each image 

%% Create scene from original image first
% create scene from file
[s1, ~] = sceneFromFile(spoke, 'rgb');
s1 = sceneSet(s1, 'fov', 2);

ieAddObject(s1);
sceneWindow;
% 
% [s2, ~] = sceneFromFile(spoke_blur, 'rgb');
% s2 = sceneSet(s2, 'fov', 2);
% 
% ieAddObject(s2);
% sceneWindow;

%% Build oi (retinal image) 
oi = oiCreate;
oi = oiCompute(oi, s1);

ieAddObject(oi);

%% Build a default cone mosaic and compute isomerizatoins

% Create the coneMosaic object
cMosaic = coneMosaic;

% Set size to show a fraction of the scene. Speeds things up.
val = 1; % 1 means full scene
cMosaic.setSizeToFOV(val * sceneGet(s1, 'fov'));

%% Set integration time 
cMosaic.integrationTime = 0.1;

%% Compute isomerizations for each eye position.
cMosaic.compute(oi);

%% Bring up a window so that we can look at things.
%
% Using the pull down in the window, you can look at
% the mosaic, the isomerizations for one fixation, or
% the movie of fixations.
cMosaic.window;

%% Get image from cone absorptions
disp('Reconstructing Image from Cone Absorptions...')
search_radii = [2, 3, 5]; % in pixels
[im1] = imageFromConeAbs(cMosaic, search_radii);

figure()
imd1 = im2double(im1);
imd1 = (imd1 - min(imd1(:))) / (max(imd1(:)) - min(imd1(:)));

% 
% imd1 = imgaussfilt(imd1, sigma);


imshow(imd1);
title('Reconstructed Original Image');


%% Repeat process for blurred scene image

%% Create scene from blurred image
[s2, ~] = sceneFromFile(spoke_blur, 'rgb');
s2 = sceneSet(s2, 'fov', 2);

ieAddObject(s2);
sceneWindow;

%% Build oi (retinal image) 
oi = oiCreate;
oi = oiCompute(oi, s2);

ieAddObject(oi);

%% Build a default cone mosaic and compute isomerizatoins

% Create the coneMosaic object
cMosaic = coneMosaic;

% Set size to show a fraction of the scene. Speeds things up.
val = 1; % 1 means full scene
cMosaic.setSizeToFOV(val * sceneGet(s2, 'fov'));

%% Set integration time 
cMosaic.integrationTime = 0.1;

%% Compute isomerizations for each eye position.
cMosaic.compute(oi);

%% Bring up a window so that we can look at things.
%
% Using the pull down in the window, you can look at
% the mosaic, the isomerizations for one fixation, or
% the movie of fixations.
cMosaic.window;

%% Get image from cone absorptions
disp('Reconstructing Image from Cone Absorptions...')
search_radii = [2, 3, 5]; % in pixels
[im2] = imageFromConeAbs(cMosaic, search_radii);

figure()
imd2 = im2double(im2);
imd2 = (imd2 - min(imd2(:))) / (max(imd2(:)) - min(imd2(:)));

%
% imd2 = imgaussfilt(imd2, sigma);

imshow(imd2);
title('Reconstructed Blurred Image');

%% Find delta E between reconstructed images 
disp('Getting Delta E between the images');

% get delta E in CIELAB
whiteRGB = [1, 1, 1] * 0.5;

img1_XYZ = rgb2xyz(imd1);
img2_XYZ = rgb2xyz(imd2);
whiteXYZ = rgb2xyz(whiteRGB);
[dE, lab1, lab2] = deltaLab(img1_XYZ,img2_XYZ, whiteXYZ);
dE_LAB = im2double(dE);

% get Delta E in S-CIELAB
fov = sceneGet(s1, 'horizontal fov'); % in degree
sampPerDeg = sceneGet(s1, 'cols') / fov;
img1_XYZ = rgb2xyz(imd1);
img2_XYZ = rgb2xyz(imd2);
dE = scielab(sampPerDeg, img1_XYZ, img2_XYZ, whiteXYZ, 'xyz');
dE_SLAB = im2double(dE);

fig_title = sprintf('Delta E Image, sigma = %.1f', sigma);
figure('name', fig_title)
subplot(2, 2, 1);
imshow(imd1);
xlabel('Original');
subplot(2, 2, 2);
imshow(imd2);
xlabel('Blurred');
subplot(2, 2, 3);
imshow(dE_LAB, [min(dE_LAB(:)), max(dE_LAB(:))]);
xlabel(['Delta E. LAB.', 'Max: ', num2str(max(dE_LAB(:)))]);
colorbar;
subplot(2, 2, 4);
imshow(dE_SLAB, [min(dE_SLAB(:)), max(dE_SLAB(:))]);
xlabel(['Delta E. S-LAB.', 'Max: ', num2str(max(dE_SLAB(:)))]);
colorbar;

%% 