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
addpath(genpath('/Users/Joe1/isetbio/'));
clc;
clear;
close all;
ieInit;

%% Build a simple scene and oi (retinal image) for computing

% First the scene
s = sceneCreate('rings rays'); % or zone plate
s = sceneSet(s, 'fov', 1);

ieAddObject(s);
sceneWindow;


% Then the oi
oi = oiCreate;
oi = oiCompute(oi, s);

ieAddObject(oi);

%% Build a default cone mosaic and compute isomerizatoins

% Create the coneMosaic object
cMosaic = coneMosaic;

% Set size to show about half the scene. Speeds things up.
cMosaic.setSizeToFOV(0.5 * sceneGet(s, 'fov'));

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
[im] = imageFromConeAbs(cMosaic);

figure(3)
imd = im2double(im);
imd = (imd - min(imd(:))) / (max(imd(:)) - min(imd(:)));
imshow(imd);

l_abs = im(:,:,1);
m_abs = im(:,:,2);
s_abs = im(:,:,3);

% Weight the different channels, see what happens there 
% Show different channels individually 

%% White balancing using Gray World

% https://www.mathworks.com/help/images/comparison-of-auto-white-balance-algorithms.html

illuminant_gw1 = illumgray(imd, 0);

B_gw1 = chromadapt(imd, illuminant_gw1, 'ColorSpace', 'linear-rgb');

B_gw1_sRGB = lin2rgb(B_gw1);

imshow(B_gw1_sRGB)
title('White balanced image using Gray World with percentiles=[0 0]')
