%% Initialize and clear
%addpath(genpath('/Users/Joe1/isetbio/'));
clc;
clear;
close all;
ieInit;

%% Build a simple scene and oi (retinal image) for computing

% First the scene
fov = 4;
[s1, cMosaic1] = createMacbethScene(fov);
[s2, cMosaic2] = createMacbethScene(fov);

%% Get image from cone absorptions
WHITE_BALANCE = false;
rec1 = reconstructImage(cMosaic1, WHITE_BALANCE);
rec2 = reconstructImage(cMosaic2, WHITE_BALANCE);

%% Get the error image
mseRGB = getMSE(rec1, rec2);

%% plot centers in macbeth reconstructed image, just to check that they are right
r0 = 71;
c0 = 22;
r_step = 51;
c_step = 51;
margin = 20;

squares1 = evaluateMacbeth(rec1, r0, c0, r_step, c_step, margin, true);
squares2 = evaluateMacbeth(rec2, r0, c0, r_step, c_step, margin, true);

%% get the MSE and SNR in all three channels for every one of the 24 squares
mse_macbeth_rgb = zeros(24, 3);
snr_macbeth_rgb = zeros(24, 3);
for i = 1:24
    mse_macbeth_rgb(i, :) = getMSE(squares1{i}, squares2{i});
    snr_macbeth_rgb(i, :) = getConstantSNR(squares1{i});
end

%% create a labelled macbeth
plotMacbethLabels(rec1, r0, c0, r_step, c_step);

%% plot the MSE for each color
figure('units','normalized','outerposition',[0.1 0.1 0.7 0.5]);
stem(mse_macbeth_rgb(:, 1), 'r', 'filled');
hold all
stem(mse_macbeth_rgb(:, 2), 'g', 'filled');
stem(mse_macbeth_rgb(:, 3), 'b', 'filled');
xticks(1:24);
grid;
title('MSE');

%% plot the SNR for each color

figure('units','normalized','outerposition',[0.1 0.1 0.7 0.5]);
stem(snr_macbeth_rgb(:, 1), 'r', 'filled');
hold all
stem(snr_macbeth_rgb(:, 2), 'g', 'filled');
stem(snr_macbeth_rgb(:, 3), 'b', 'filled');
xticks(1:24);
ylim([0, 140]);
grid;
title('SNR');

%% get the pixelwise CIELAB delta E for the two reconstructed images
if WHITE_BALANCE
    whiteRGB = [1, 1, 1];
else
    whiteRGB = [0.9671, 0.8751, 0.1795];
end

% optional: manually adjust white point
whiteRGB = [1, 1, 1];

dE1 = getDeltaE_CIELAB(rec1, rec2, whiteRGB);

% get average error for each of the 24 colors
errorSquares = evaluateMacbeth(dE1, r0, c0, r_step, c_step, margin, false);
avgError = zeros(24, 1);
for i = 1:24
    avgError(i) = mean(errorSquares{i}(:));
end

%% plot the delta E error with stem plot
figure('units','normalized','outerposition',[0.1 0.1 0.7 0.5]);
stem(avgError, 'k', 'filled');
xticks(1:24);
ylim([0, 30]);
grid;
title(['Average Delta E (CIELAB), WhitePoint ', mat2str(whiteRGB)]);


% %% save stuff
% rootdir = 'NoiseExperiments/';
% ex_name = 'Balanced_Reconstruction';
% exdir = strcat(rootdir, ex_name, '/');
% mkdir(exdir);
% 
% 
% saveas(figure(3), strcat(exdir, 'rec1.jpg'));
% saveas(figure(4), strcat(exdir, 'rec2.jpg'));
% saveas(figure(5), strcat(exdir, 'rec1box.jpg'));
% saveas(figure(6), strcat(exdir, 'rec2box.jpg'));
% 
% saveas(figure(8), strcat(exdir, 'mse.jpg'));
% saveas(figure(8), strcat(exdir, 'mse.fig'));
% saveas(figure(9), strcat(exdir, 'snr.jpg'));
% saveas(figure(9), strcat(exdir, 'snr.fig'));
% saveas(figure(10), strcat(exdir, 'deltaE.jpg'));
% saveas(figure(10), strcat(exdir, 'deltaE.fig'));


%% Find delta E between reconstructed images (in RGB)
function dE = getDeltaE_CIELAB(img1, img2, whiteRGB)
    img1_XYZ = rgb2xyz(img1);
    img2_XYZ = rgb2xyz(img2);
    whiteXYZ = rgb2xyz(whiteRGB);
    [dE, ~, ~] = deltaLab(img1_XYZ,img2_XYZ, whiteXYZ);
    dE = im2double(dE);
end


%% label the macbeth chart
function plotMacbethLabels(img_in, r0, c0, r_step, c_step)
    [m, n, c] = size(img_in);
    % get scale factor so the squares stay centered
    scale = m / 297;
    r0 = round(r0 * scale);
    c0 = round(c0 * scale);
    r_step = round(r_step * scale);
    c_step = round(c_step * scale);

    counter = 1;
    for i = 0:3
        for j = 0:5
            r = r0 + i * r_step;
            c = c0 + j * c_step;
            
            img_in = insertText(img_in, [c-5, r-5], num2str(counter));
            counter = counter + 1;
        end
    end
    figure()
    imshow(img_in);
end

%% label the macbeth chart as
function labelChart()
    chart = imread('NoiseExperiments/colorchart.png');
    chart = imresize(chart, 8);
    
    [m, n, c] = size(chart);
    r_step = round(m / 4);
    c_step = round(n / 6);
    
    
    counter = 1;
    for i = 0:3
        for j = 0:5
            r = round(r_step/2) + i * r_step;
            c = round(c_step/2) + j * c_step;
            
            chart = insertText(chart, [c-5, r-5], num2str(counter));
            counter = counter + 1;
        end
    end
    figure();
    imshow(chart);
end

%% SNR (mean / std) for each channel, assuming a CONSTANT IMAGE
function snrRGB = getConstantSNR(img)
    R = img(:, :, 1);
    G = img(:, :, 2);
    B = img(:, :, 3);
    
    snrR = mean(R(:)) / std(R(:));
    snrG = mean(G(:)) / std(G(:));
    snrB = mean(B(:)) / std(B(:));
    
    snrRGB = [snrR, snrG, snrB];
end

%% channel wise MSE between two images
function mseRGB = getMSE(img1, img2)
    error = abs(img1-img2);
    R = error(:, :, 1);
    G = error(:, :, 2);
    B = error(:, :, 3);
    mseR = mean(R(:).^2);
    mseG = mean(G(:).^2);
    mseB = mean(B(:).^2);
    mseRGB = [mseR, mseG, mseB];
end

%% plot dots in the middle of image
function squares = evaluateMacbeth(imd, r0, c0, r_step, c_step, margin, PLOT)

    imd_test = imd;
    
    squares = cell(24);
    
    if ndims(imd) == 3
        [m, n, num_channels] = size(imd);
    else
        [m, n] = size(imd);
        num_channels = 1;
    end
    
    % get scale factor so the squares stay centered
    scale = m / 297;
    r0 = round(r0 * scale);
    c0 = round(c0 * scale);
    r_step = round(r_step * scale);
    c_step = round(c_step * scale);
    margin = round(margin * scale);
    
    avgs = zeros(24, num_channels);

    count = 1;
    margin2 = margin - 2;
    for i = 0:3
        for j = 0:5
            r = r0 + i * r_step;
            c = c0 + j * c_step;
            
            % add frame
            imd_test(r-margin:r+margin, c-margin:c+margin, :) = ...
                imd_test(r-margin:r+margin, c-margin:c+margin, :) + 10;
            
            imd_test(r-margin2:r+margin2, c-margin2:c+margin2, :) = ...
                imd_test(r-margin2:r+margin2, c-margin2:c+margin2, :) - 10;
            
            % clip to max of 1
            imd_test(imd_test>1) = 1;
            
            squares{count} = imd(r-margin:r+margin, c-margin:c+margin, :);
            avgs(count, :) = mean(reshape(squares{count}, [], num_channels));
            count = count + 1;
        end
    end

    if PLOT
        figure();
        imshow(imd_test);
    end
end

%% reconstruction function
function rgb = reconstructImage(cMosaic, WHITE_BALANCE)
    disp('Reconstructing Image from Cone Absorptions...')
    
    search_radii = [2, 3, 5]; % in pixels
    [im] = imageFromConeAbs(cMosaic, search_radii);
    imd = im2double(im);
    imd = imd / max(imd(:));

    rgb = lms2srgb(imd);

    if WHITE_BALANCE
        r0 = 71;
        c0 = 22;
        r_step = 51;
        c_step = 51;
        margin = 20;
        squares = evaluateMacbeth(rgb, r0, c0, r_step, c_step, margin, false);
        
        % square 19 is reference white
        whiteRGB = squares{19};
        whiteRGB = reshape(whiteRGB, [], 3);
        avgWhiteRGB = mean(whiteRGB, 1)
        factorRGB  = [1, 1, 1] ./ avgWhiteRGB;
        
        % apply
        for i = 1:3
            rgb(:, :, i) = rgb(:, :, i) * factorRGB(i);
        end
        
        % clip
        rgb(rgb>1) = 1;
        rgb(rgb<0) = 0;
    end

    figure();
    imshow(rgb);
end

%% create scenes and absorptions
function [s, cMosaic] = createMacbethScene(fov)
    s = sceneCreate('macbeth'); % or zone plate
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
    cMosaic.integrationTime = 1000.00;

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
