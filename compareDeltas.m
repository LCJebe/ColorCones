%% Initialize and clear
%addpath(genpath('/Users/Joe1/isetbio/'));
%addpath(genpath('/Users/Joe1/Documents/MATLAB/'));
clc; clear; 
close all;
ieInit;

%% create our own scene image

disp('Setting up Scene...');

scene = 'square'; % or 'square'
scene_size = 512;

col1 = [1,1,0];
col2 = [1,1,1];

wbal_on = 1; % set to 1 to do white balancing 

if strcmp(scene, 'spoke')
    % create spoke
    scene = createSpoke(col1, col2, 10, scene_size);
elseif strcmp(scene, 'square')
    % create square
    scene = createSquare(col1, scene_size);
end

    

%% Modify original scene

mod_scene = 4; 

if mod_scene == 1
%   Option 1: Blur scene
    sigma = 3;
    scene_mod = imgaussfilt(spoke, sigma);
elseif mod_scene == 2
    % Option 2: Don't change the scene (used for checking noise) 
    scene_mod = scene;
elseif mod_scene == 3
    % Option 3: Switch colors of spoke
    scene_mod = createSpoke(col2, col1, 10, spoke_size);
elseif mod_scene == 4
    % Option 4: Create new uniform color square
    scene_mod = createSquare(col2, scene_size);
end


scene = im2double(scene);
scene_mod = im2double(scene_mod);


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
sceneXYZ = rgb2xyz(scene);
scene_modXYZ = rgb2xyz(scene_mod);
%imgRGB = dac2rgb(img01, gammaTable);
%imgLMS = changeColorSpace(imgRGB, rgb2lms);

% define image input format for S-CIELAB
%imageformat = 'lms';


%% Perform S-CIELAB calculation
%filtOpp = scielab(sampPerDeg, imgLMS, imageformat); 
%filtXYZ = changeColorSpace(filtOpp, cmatrix('opp2xyz')); 

% get S-CIELAB delta E as difference to the white point of the image
whitepointRGB = [1, 1, 1];
whitepointXYZ = rgb2xyz(whitepointRGB);

errorImage = scielab(sampPerDeg, sceneXYZ, scene_modXYZ, whitepointXYZ, 'xyz');


%% analyze S-CIELAB image (it's not really an error image, just a delta E
    % image with respect to a reference white image

%% show images
fig_title = sprintf('Original & modified scenes ([%.1f %.1f %.1f],[%.1f %.1f %.1f])',col1(1), col1(2), col1(3), col2(1), col2(2), col2(3));
figure('name', fig_title)
subplot(2, 2, 1);
imshow(scene);
subplot(2, 2, 2);
imshow(scene_mod);
subplot(2, 2, 3);
imshow(errorImage, [0, max(errorImage(:))]);
mean_dE_S = mean(mean(errorImage(:)));
xlabel(['Delta E. S-LAB,', ' Avg: ', num2str(mean_dE_S)]);

colorbar;

subplot(2, 2, 4);
% delta E (cielab) here
whiteRGB = [1, 1, 1];
img1_XYZ = rgb2xyz(scene);
img2_XYZ = rgb2xyz(scene_mod);
whiteXYZ = rgb2xyz(whiteRGB);
[dE, lab1, lab2] = deltaLab(img1_XYZ,img2_XYZ, whiteXYZ);
dE_LAB = im2double(dE);

imshow(dE_LAB, [0, max(dE_LAB(:))]);
xlabel(['Delta E. LAB,', ' Avg: ', num2str(mean(mean(dE_LAB(:))))]);
colorbar;


%% Now calculate cone responses from each image 

%% Create scene from original image first
% create scene from file
[s1, ~] = sceneFromFile(scene, 'rgb');
s1 = sceneSet(s1, 'fov', 2);

ieAddObject(s1);
sceneWindow;


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
cMosaic.integrationTime = 0.05;

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

if max(imd1(:)) ~= 0
    % need this if statement when dealing with uniformly black square
    imd1 = imd1 / (max(imd1(:)));
end

imd1 = lms2srgb(imd1);

% get white point
% whiteRGB = reshape(imd1,[],3);
% wbal = mean(whiteRGB, 1);
wbal = [0.9725 0.8864 0.1830];

% Adjust channels for white point
if wbal_on 
    for i = 1:3
        imd1(:,:,i) = imd1(:,:,i) ./ wbal(i);
    end
    imd1(imd1>1) = 1;
    imd1(imd1<0) = 0; 
end
% 
% imd1 = imgaussfilt(imd1, sigma);


imshow(imd1);

fig_title = sprintf('Reconstructed Original Image ([%.1f %.1f %.1f])',col1(1), col1(2), col1(3));
title(fig_title);


%% Repeat process for modified scene image

%% Create scene from modified image
[s2, ~] = sceneFromFile(scene_mod, 'rgb');
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
cMosaic.integrationTime = 0.05;

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


if (max(imd2(:))) ~= 0
    imd2 = (imd2) / (max(imd2(:)));
end

imd2 = lms2srgb(imd2); 

if wbal_on 
    for i = 1:3
        imd2(:,:,i) = imd2(:,:,i) ./ wbal(i);
    end
    imd2(imd2>1) = 1;
    imd2(imd2<0) = 0; 
end

%
% imd2 = imgaussfilt(imd2, sigma);

imshow(imd2);
fig_title = sprintf('Reconstructed Modified Image ([%.1f %.1f %.1f])',col2(1), col2(2), col2(3));
title(fig_title);


%% Find delta E between reconstructed images 
disp('Getting Delta E between the images');

% get delta E in CIELAB
whiteRGB_1= [1, 1, 1]; 

whiteRGB_2 = wbal; 

img1_XYZ = rgb2xyz(imd1);
img2_XYZ = rgb2xyz(imd2);

% get dE with white point = [1,1,1]
whiteXYZ_1 = rgb2xyz(whiteRGB_1);
[dE, lab1, lab2] = deltaLab(img1_XYZ,img2_XYZ, whiteXYZ_1);
dE_LAB = im2double(dE);

% get dE with white point = wbal
whiteXYZ_2 = rgb2xyz(whiteRGB_2);
[dE2, lab3, lab4] = deltaLab(img1_XYZ,img2_XYZ, whiteXYZ_2);
dE_LAB2 = im2double(dE2);

mean_dE2 = mean(mean(dE_LAB2(:)));


% % get Delta E in S-CIELAB
% fov = sceneGet(s1, 'horizontal fov'); % in degree
% sampPerDeg = sceneGet(s1, 'cols') / fov;
% img1_XYZ = rgb2xyz(imd1);
% img2_XYZ = rgb2xyz(imd2);
% dE = scielab(sampPerDeg, img1_XYZ, img2_XYZ, whiteXYZ, 'xyz');
% dE_SLAB = im2double(dE);

fig_title = sprintf('Deltat E Image([%.1f %.1f %.1f],[%.1f %.1f %.1f])',col1(1), col1(2), col1(3), col2(1), col2(2), col2(3));
figure('name', fig_title)
subplot(2, 2, 1);
imshow(imd1);
xlabel('Original');
subplot(2, 2, 2);
imshow(imd2);
xlabel('Modified');
subplot(2, 2, 3);
imshow(dE_LAB, [min(dE_LAB(:)), max(dE_LAB(:))]);
mean_dE = mean(mean(dE_LAB(:)));
xlabel(['Delta E. LAB,', ' Avg: ', num2str(mean_dE)]);
colorbar;
disp('')
disp('Mean dE between Cone Responses (WP = [1,1,1]): ');
disp(mean_dE);
disp('Mean dE between Cone Responses (WP = wbal): ');
disp(mean_dE2);

% subplot(2, 2, 4);
% imshow(dE_SLAB, [min(dE_SLAB(:)), max(dE_SLAB(:))]);
% xlabel(['Delta E. S-LAB.', 'Max: ', num2str(max(dE_SLAB(:)))]);
% colorbar;

%% Ratios of dE for cone responses to dE S-CIELAB
dE_ratio1 = mean_dE / mean_dE_S;
disp('dE ratio (WP = [1 1 1]');
disp(dE_ratio1);
dE_ratio2 = mean_dE2 / mean_dE_S;
disp('dE ratio (WP = wbal');
disp(dE_ratio2);


%%
