clear all
close all
%% define Experiments to be run
n_exps = 10;
exps = cell(n_exps);
for i = 1:n_exps
    exps{i} = struct;
    exps{i}.n_rings = 10;
    exps{i}.spoke_size = 512;
end

% first exp: yellow white vs. white yellow spoke
exps{1}.name = "Red";
exps{1}.col1 = [1, 0, 0];
exps{1}.col2 = getCol2(exps{1}.col1);

% second exp: black white vs. white black spoke
exps{2}.name = "Yellow";
exps{2}.col1 = [1, 1, 0];
exps{2}.col2 = getCol2(exps{2}.col1);

% third experiment: red, blue vs. red, greenish blue
exps{3}.name = "Green";
exps{3}.col1 = [0, 1, 0];
exps{3}.col2 = getCol2(exps{3}.col1);

% fourth experiment: green yellow vs. green red
exps{4}.name = "Greenblue";
exps{4}.col1 = [0, 1, 1];
exps{4}.col2 = getCol2(exps{4}.col1);

% fifth experiment: black-blue opposited
exps{5}.name = "Blue";
exps{5}.col1 = [0, 0, 1];
exps{5}.col2 = getCol2(exps{5}.col1);

% sixth experiment: yellow-blue opposite
exps{6}.name = "Magenta";
exps{6}.col1 = [1, 0, 1];
exps{6}.col2 = getCol2(exps{6}.col1);

% seventh experiment: Red-GreenBlue-Opposite
exps{7}.name = "White";
exps{7}.col1 = [1, 1, 1];
exps{7}.col2 = getCol2(exps{7}.col1);



%% define white balance (estimated with whitePointCheck script)
whiteRGB_1 = [0.9671, 0.8751, 0.1795];

%% run select an experiment number and run this experiment
for ex = 7:7
    exp = exps{ex};

    for i = 1:4
        
        col1 = exp.col1;
        col2 = exp.col2(i, :);
        
        
        disp(strcat("Comparing ", num2str(col1), " and ", num2str(col2)));
        
        % create the two scenes
        fov = 2;
        [s1, cMosaic1, spoke1] = createSpokeScene(col1, col2, exp.spoke_size, exp.n_rings, fov);
        [s2, cMosaic2, spoke2] = createSpokeScene(col2, col1, exp.spoke_size, exp.n_rings, fov);

        rec1 = reconstructImage(cMosaic1, [1, 1, 1]);
        rec2 = reconstructImage(cMosaic2, [1, 1, 1]);

        % get CIELAB difference for the two reconstructed images
        dE_white1 = getDeltaE_CIELAB(rec1, rec2, whiteRGB_1);

        % Now we have three different delta E. Still missing: Scene S-CIELAB!!!
        fov = sceneGet(s1, 'horizontal fov'); % in degree
        sampPerDeg = sceneGet(s1, 'cols') / fov;
        dE_SLAB = getDeltaE_SLAB(spoke1, spoke2, whiteRGB_1, sampPerDeg);

        %% Plot the two spokes and the four reconstructed images

        f1 = figure('units','normalized','outerposition',[0 0 1 1]);
        subplot(2, 2, 1);
        imshow(spoke1);
        title('Input Spokes');
        subplot(2, 2, 2);
        imshow(spoke2);
        subplot(2, 2, 3);
        imshow(rec1);
        title('Reconstruction')
        subplot(2, 2, 4);
        imshow(rec2);

        %% new figure: plot the four differences

        % find out maximum first
        maximum = max([max(dE_SLAB(:)), max(dE_white1(:))]);

        
        f2 = figure('units','normalized','outerposition',[0 0 1 1]);
        subplot(1, 2, 1);
        imshow(dE_SLAB, [0, maximum]);
        title('S-CIELAB (Scenes)');
        xlabel(['Mean: ', num2str(mean(dE_SLAB(:)))]);
        colorbar;
        subplot(1, 2, 2);
        imshow(dE_white1, [0, maximum]);
        title(['CIELAB (Cones)']);
        xlabel(['Mean: ', num2str(mean(dE_white1(:)))]);
        colorbar;


        %% Saves
        % create directory
        savedir = 'NewSpokeExperiments/';
        expdir = strcat(savedir, exp.name, '/');
        mkdir(expdir);
        saveas(f1, strcat(expdir, 'images_', num2str(i), '.jpg'));
        saveas(f1, strcat(expdir, 'images_', num2str(i), '.fig'));
        saveas(f2, strcat(expdir, 'erorrs_', num2str(i), '.jpg'));
        saveas(f2, strcat(expdir, 'errors_', num2str(i), '.fig'));
        
        close all
    end
end

%% check
col1 = [1, 0, 0];
getCol2(col1)

%% get the 4 experiments depending on the first color
function col2 = getCol2(col1)
    col2 = repmat(col1, 4, 1);
    if col1(1) == 0 col2(1, 1) = 0.05; else col2(1, 1) = 0.95; end
    if col1(2) == 0 col2(2, 2) = 0.05; else col2(2, 2) = 0.95; end
    if col1(3) == 0 col2(3, 3) = 0.05; else col2(3, 3) = 0.95; end
    
    col2(4, :) = [col2(1, 1), col2(2, 2), col2(3, 3)];
end

%% scales each color channel by a factor
function img = whiteBalance(img, whiteRGB)
    for i = 1:3
        img(:, :, i) = img(:, :, i) / whiteRGB(i);
    end
    
    img(img>1) = 1;
    img(img<0) = 0;
end

%% Find delta E between reconstructed images (in RGB)
function dE = getDeltaE_CIELAB(img1, img2, whiteRGB)
    img1_XYZ = rgb2xyz(img1);
    img2_XYZ = rgb2xyz(img2);
    whiteXYZ = rgb2xyz(whiteRGB);
    [dE, ~, ~] = deltaLab(img1_XYZ,img2_XYZ, whiteXYZ);
    dE = im2double(dE);
end


%% reconstruction function
function rgb = reconstructImage(cMosaic, whiteRGB)
    disp('Reconstructing Image from Cone Absorptions...')
    
    search_radii = [2, 3, 5]; % in pixels
    [im] = imageFromConeAbs(cMosaic, search_radii);
    imd = im2double(im);
    imd = imd / max(imd(:));

    rgb = lms2srgb(imd);

    % apply white balance
    for i = 1:3
        rgb(:, :, i) = rgb(:, :, i) / whiteRGB(i);
    end
        
    % clip
    rgb(rgb>1) = 1;
    rgb(rgb<0) = 0;

    figure();
    imshow(rgb);
end


%% create scenes and absorptions
function [s, cMosaic, spoke] = createSpokeScene(col1, col2, spoke_size, n_rings, fov)

    spoke = createSpoke(col1, col2, n_rings, spoke_size);

    % create scene from file
    [s, ~] = sceneFromFile(spoke, 'rgb');
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
    cMosaic.integrationTime = 1000.0;

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


function dE_d = getDeltaE_SLAB(img1RGB, img2RGB, whiteRGB, sampPerDeg)
    img1_XYZ = rgb2xyz(img1RGB);
    img2_XYZ = rgb2xyz(img2RGB);
    whiteXYZ = rgb2xyz(whiteRGB);
    
    % get S-CIELAB Error image
    dE = scielab(sampPerDeg, img1_XYZ, img2_XYZ, whiteXYZ, 'xyz');
    dE_d = im2double(dE);
end
