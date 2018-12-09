%% Computation of S-CIELAB Values given an input image
clear all
close all

%% read in image
blur = 'gauss'; % 'gauss' or 'lpf' % LPF doesn't work yet. WIP

img = createSpoke([1, 1, 0], [0, 0, 1], 64, 256);

if strcmp(blur, 'gauss')
    img_blur = imgaussfilt(img, 1);
elseif strcmp(blur, 'lpf')
    k = 1;
    % do fft on three channels individually
    img_fft = fftshift(fft2(img));
    
    % remove high frequencies
    img_fft_filt = zeros(size(img));
    center = size(img) / 2 + 0.5;
    
    % coordinates to crop out middle low frequency part
    m = min(size(img));
    y_min = round(center(1) - k*m);
    y_max = round(center(1) + k*m);
    x_min = round(center(2) - k*m);
    x_max = round(center(2) + k*m);
    
    % crop
    img_fft_filt(y_min:y_max, x_min:x_max) = img_fft(y_min:y_max, x_min:x_max);
    
    % retransform to spacial domain
    img_blur = real(ifft2(fftshift(img_fft_filt)));
end

img = im2double(img);
img_blur = im2double(img_blur);

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
imgXYZ = rgb2xyz(img);
img_blurXYZ = rgb2xyz(img_blur);
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

errorImage = scielab(sampPerDeg, imgXYZ, img_blurXYZ, whitepointXYZ, 'xyz');

%% analyze S-CIELAB image (it's not really an error image, just a delta E
% image with respect to a reference white image

%% show images
figure();
subplot(1, 3, 1);
imshow(img);
subplot(1, 3, 2);
imshow(img_blur);
subplot(1, 3, 3);
imshow(errorImage, [0, max(errorImage(:))]);
