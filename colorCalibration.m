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
s = sceneCreate('macbeth'); % or zone plate
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
cMosaic.integrationTime = 0.2;

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

%% plot dots in the middle of IMD

[m, n, ~] = size(imd);
imd_test = imd;

captures = cell(24);
avgs = zeros(24, 3);
radius = 20;

count = 1;
for i = 0:3
    for j = 0:5
        r = 71 + i * 51;
        c = 22 + j * 51;
        imd_test(r-2:r+2, c-2:c+2, 1) = 1;
        captures{count} = imd(r-radius:r+radius, c-radius:c+radius, :);
        avgs(count, :) = mean(reshape(captures{count}, [], 3));
        count = count + 1;
    end
end

% get targets
hexTargets = {'#735244', '#c29682', '#627a9d', '#576c43', '#8580b1', '#67bdaa', ...
              '#d67e2c', '#505ba6', '#c15a63', '#5e3c6c', '#9dbc40', '#e0a32e', ...
              '#383d96', '#469449', '#af363c', '#e7c71f', '#bb5695', '#0885a1', ...
              '#f3f3f2', '#c8c8c8', '#a0a0a0', '#7a7a79', '#555555', '#343434'};
          
rgbTargets = zeros(24, 3);
for i = 1:24
    rgbTargets(i, :) = hex2rgb(hexTargets{i});
end


%% linear regression for colors
% select colors to adjust. 1:24 for all
adjust = 13:18;


illum_R = avgs(adjust, 1) \ rgbTargets(adjust, 1);
illum_G = avgs(adjust, 2) \ rgbTargets(adjust, 2);
illum_B = avgs(adjust, 3) \ rgbTargets(adjust, 3);

illum = reshape([illum_R, illum_G, illum_B], 3, 1 );
inv_illum = [1; 1; 1] ./ illum;
inv_illum = inv_illum / norm(inv_illum);

% adjust image
imd_adapt =  imd;
imd_adapt(:, :, 1) = imd(:, :, 1) * illum_R;
imd_adapt(:, :, 2) = imd(:, :, 2) * illum_G;
imd_adapt(:, :, 3) = imd(:, :, 3) * illum_B;


%imd_adapt = scaleToOne(imd_adapt, [3, 3]);
imd_adapt(imd_adapt>1) = 1;
imd_adapt(imd_adapt<0) = 0;


% show image
imshow(imd_test);
figure()
imshow(imd_adapt);


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

