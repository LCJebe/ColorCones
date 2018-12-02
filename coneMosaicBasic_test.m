% Basic introduction to the cone mosaic object.
%
% Description:
%    Show how to create a cone mosaic object and compute cone
%    isomerizatoins across a set of small eye movements. This lets you look
%    at the result in the coneMosaic window.
%
% Notes:
%    * [Note: DHB - Either we should show how to access fields of the
%      mosaic programatically here, or we should point to a different
%      tutorial that does so.]
%

%% Initialize and clear
%addpath(genpath('/Users/Joe1/isetbio/'));
clc;
clear;
close all;
ieInit;

%% Build a simple scene and oi (retinal image) for computing

% First the scene
s = sceneCreate('rings rays'); % or zone plate
s = sceneSet(s, 'fov', 2);

%% create our own scene form a colorful spoke target

disp('Setting up Scene...');
% create spoke
col1 = [1, 1, 0];
col2 = [0, 0, 1];
spoke_size = 512;
spoke = createSpoke(col1, col2, 10, spoke_size);

% create scene from file
[s, ~] = sceneFromFile(spoke, 'rgb');
s = sceneSet(s, 'fov', 2);

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
% cMosaic.integrationTime = 0.01;

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

%% Get image from cone absorptions
disp('Reconstructing Image from Cone Absorptions...')
search_radii = [2, 3, 5]; % in pixels
[im] = imageFromConeAbs(cMosaic, search_radii);

figure(3)
imd = im2double(im);
imd = (imd - min(imd(:))) / (max(imd(:)) - min(imd(:)));
%imshow(imd);

l_abs = im(:,:,1);
m_abs = im(:,:,2);
s_abs = im(:,:,3);

% Weight the different channels, see what happens there 
% Show different channels individually 

%% Get Delta E between recnostructed and original image
disp('Getting Delta E between the images');
[n, n, m] = size(imd);
scale_down = n / spoke_size;
spoke_scaled = imresize(spoke,scale_down);

% get delta E in CIELAB
whiteRGB = [1, 1, 1] * 0.5;
dE_LAB = getDeltaE_LAB(spoke_scaled, imd, whiteRGB);

% get Delta E in S-CIELAB
fov = sceneGet(s, 'horizontal fov'); % in degree
sampPerDeg = sceneGet(s, 'cols') / fov;
dE_SLAB = getDeltaE_SLAB(spoke_scaled, imd, whiteRGB, sampPerDeg);


figure('name', 'Delta E Image')
subplot(2, 2, 1);
imshow(spoke_scaled);
xlabel('Original');
subplot(2, 2, 2);
imshow(imd);
xlabel('Recons.');
subplot(2, 2, 3);
imshow(dE_LAB, [min(dE_LAB(:)), max(dE_LAB(:))]);
xlabel(['Delta E. LAB.', 'Max: ', num2str(max(dE_LAB(:)))]);
subplot(2, 2, 4);
imshow(dE_SLAB, [min(dE_SLAB(:)), max(dE_SLAB(:))]);
xlabel(['Delta E. S-LAB.', 'Max: ', num2str(max(dE_SLAB(:)))]);


% [de00, errComponents] = deltaE2000(original_LAB,recon_LAB,[1 1 1]);

%% White balancing using Gray World

% https://www.mathworks.com/help/images/comparison-of-auto-white-balance-algorithms.html


illuminant_gw1 = illumpca(imd, 3.5);

% this is the illumination returned by white balancing the black/white
% zoneplate. Better whould be to use the Macbeth color chart...
illum = [0.817034824735336,0.573905664424434,0.055555229378056];

imd_adapt = chromadapt(imd, illum, 'ColorSpace', 'linear-rgb');

imd_adapt_sRGB = lin2rgb(imd_adapt);
imd_adapt_sRGB = imd_adapt_sRGB;
imd_adapt_sRGB = scaleToOne(imd_adapt_sRGB, [0, 0]);

% show the unbalanced and balanced image side by side
figure('name', 'White balanced image using PCA auto white balance');
subplot(1, 2, 1);
imshow(imd);
subplot(1, 2, 2);
imshow(imd_adapt_sRGB);

%% compare with plots

% show channels separately
[r, g, b] = deal(zeros(1, 1, 3));
r(1) = 1;
g(2) = 1;
b(3) = 1;

figure()
subplot(3, 4, 1);
imshow(imd.*r);
xlabel('Raw red (L)');
subplot(3, 4, 2);
imshow(imd.*g);
xlabel('Raw green (M)');
subplot(3, 4, 3);
imshow(imd.*b);
xlabel('Raw blue (S)');
subplot(3, 4, 4);
imshow(imd);
xlabel('Raw Image');


percentSaturation = [10, 10];

subplot(3, 4, 5);
imshow(imd_adapt_sRGB.*r);
xlabel('Balanced red (L)');
subplot(3, 4, 6);
imshow(imd_adapt_sRGB.*g);
xlabel('Balanced green (M)');
subplot(3, 4, 7);
imshow(imd_adapt_sRGB.*b);
xlabel('Balanced blue (S)');
subplot(3, 4, 8);
imshow(imd_adapt_sRGB);
xlabel('Balanced Image');

% scale adapted image to 0-1
imd_adapt_sRGB_scaled = scaleToOne(imd_adapt_sRGB, percentSaturation);

subplot(3, 4, 9);
imshow(imd_adapt_sRGB_scaled.*r);
xlabel('Contrast Enhanced red (L)');
subplot(3, 4, 10);
imshow(imd_adapt_sRGB_scaled.*g);
xlabel('Contrast Enhanced green (M)');
subplot(3, 4, 11);
imshow(imd_adapt_sRGB_scaled.*b);
xlabel('Contrast Enhanced blue (S)');
subplot(3, 4, 12);
imshow(imd_adapt_sRGB_scaled);
xlabel('Contrast Enhanced blue (S)');

%% finally, get S-CIELAB difference for the processed image too

% get delta E in CIELAB
dE_LAB = getDeltaE_LAB(spoke_scaled, imd_adapt_sRGB_scaled, whiteRGB);

% get Delta E in S-CIELAB
dE_SLAB = getDeltaE_SLAB(spoke_scaled, imd_adapt_sRGB_scaled, whiteRGB, sampPerDeg);

figure('name', 'Delta E S-CIELAB Image')
subplot(2, 2, 1);
imshow(spoke_scaled);
xlabel('Original');
subplot(2, 2, 2);
imshow(imd_adapt_sRGB_scaled);
xlabel('Recons.');
subplot(2, 2, 3);
imshow(dE_LAB, [min(dE_LAB(:)), max(dE_LAB(:))]);
xlabel(['Delta E. LAB.', 'Max: ', num2str(max(dE_LAB(:)))]);
subplot(2, 2, 4);
imshow(dE_SLAB, [min(dE_SLAB(:)), max(dE_SLAB(:))]);
xlabel(['Delta E. S-LAB.', 'Max: ', num2str(max(dE_SLAB(:)))]);

%%

function img = scaleToOne(img, percentSaturation)
    bot = percentSaturation(1) / 100;
    top = percentSaturation(2) / 100;
    
    idx_bot  = floor(max(1, length(img(:)) * bot));
    idx_top  = ceil(min(length(img(:)) * (1-top), length(img(:))));
    
    img_sorted = sort(img(:));
    ref_min = img_sorted(idx_bot);
    ref_max = img_sorted(idx_top);
    
    
    img = (img - ref_min) / (ref_max - ref_min);
    
    img(img>1) = 1;
    img(img<0) = 0;
end

function dE_d = getDeltaE_LAB(img1RGB, img2RGB, whiteRGB)
    img1_XYZ = rgb2xyz(img1RGB);
    img2_XYZ = rgb2xyz(img2RGB);
    whiteXYZ = rgb2xyz(whiteRGB);

    %% Calculate delta E between original & reconstructed images

    [dE, lab1, lab2] = deltaLab(img1_XYZ,img2_XYZ, whiteXYZ);
    % dE_XYZ = ieLAB2XYZ(dE, whitepointXYZ);

    dE_d = im2double(dE);
end

function dE_d = getDeltaE_SLAB(img1RGB, img2RGB, whiteRGB, sampPerDeg)
    img1_XYZ = rgb2xyz(img1RGB);
    img2_XYZ = rgb2xyz(img2RGB);
    whiteXYZ = rgb2xyz(whiteRGB);
    
    % get S-CIELAB Error image
    dE = scielab(sampPerDeg, img1_XYZ, img2_XYZ, whiteXYZ, 'xyz');
    dE_d = im2double(dE);
end