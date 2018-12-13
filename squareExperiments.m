clear all
close all

%% define Experiments to be run
n_exps = 28;  % Warning: currently takes about 15 min to run all 28 experiments.

exps = cell(n_exps);
for i = 1:n_exps
    exps{i} = struct;
    exps{i}.size = 512;
end

% Experiments 1-4: Red square 
col1 = [1,0,0];
col2 = getCol2(col1);
j = 1;
for i = 1:4
    exps{j}.col1 = col1;
    exps{j}.col2 = col2(i,:);
    j = j+1;
end
exps{1}.name = 'Red-dR';
exps{2}.name = 'Red-dG';
exps{3}.name = 'Red-dB';
exps{4}.name = 'Red-dAll';

% Exp. 5-8: Yellow Square
col1 = [1,1,0];
col2 = getCol2(col1);
for i = 1:4
    exps{j}.col1 = col1;
    exps{j}.col2 = col2(i,:);
    j = j+1;
end
exps{5}.name = 'Yellow-dR';
exps{6}.name = 'Yellow-dG';
exps{7}.name = 'Yellow-dB';
exps{8}.name = 'Yellow-dAll';

% Exp 9-12: Green 
col1 = [0,1,0];
col2 = getCol2(col1);
for i = 1:4
    exps{j}.col1 = col1;
    exps{j}.col2 = col2(i,:);
    j = j+1;
end
exps{9}.name = 'Green-dR';
exps{10}.name = 'Green-dG';
exps{11}.name = 'Green-dB';
exps{12}.name = 'Green-dAll';

% Exp 13-16: Green-blue 
col1 = [0,1,1];
col2 = getCol2(col1);
for i = 1:4
    exps{j}.col1 = col1;
    exps{j}.col2 = col2(i,:);
    j = j+1;
end
exps{13}.name = 'GreenBlue-dR';
exps{14}.name = 'GreenBlue-dG';
exps{15}.name = 'GreenBlue-dB';
exps{16}.name = 'GreenBlue-dAll';

% Exp 17-20: Blue 
col1 = [0,0,1];
col2 = getCol2(col1);
for i = 1:4
    exps{j}.col1 = col1;
    exps{j}.col2 = col2(i,:);
    j = j+1;
end
exps{17}.name = 'Blue-dR';
exps{18}.name = 'Blue-dG';
exps{19}.name = 'Blue-dB';
exps{20}.name = 'Blue-dAll';

% Exp 21-24: Violet  
col1 = [1,0,1];
col2 = getCol2(col1);
for i = 1:4
    exps{j}.col1 = col1;
    exps{j}.col2 = col2(i,:);
    j=j+1;
end
exps{21}.name = 'Violet-dR';
exps{22}.name = 'Violet-dG';
exps{23}.name = 'Violet-dB';
exps{24}.name = 'Violet-dAll';

% Exp 25-28: White  
col1 = [1,1,1];
col2 = getCol2(col1);
for i = 1:4
    exps{j}.col1 = col1;
    exps{j}.col2 = col2(i,:);
    j=j+1;
end
exps{25}.name = 'White-dR';
exps{26}.name = 'White-dG';
exps{27}.name = 'White-dB';
exps{28}.name = 'White-dAll';


%% define white balance (estimated with whitePointCheck script)
whiteRGB_1 = [0.9671, 0.8751, 0.1795];

%% run select an experiment number and run this experiment
for ex = 1:n_exps
    close all
    exp = exps{ex};

    % create the two scenes
    fov = 2;
    [s1, cMosaic1, square1] = createSquareScene(exp.col1, exp.size, fov);
    [s2, cMosaic2, square2] = createSquareScene(exp.col2, exp.size, fov);

    rec1 = reconstructImage(cMosaic1, [1,1,1]);
    rec2 = reconstructImage(cMosaic2, [1,1,1]);

    rec1_avg = mean(mean(rec1));
    rec2_avg = mean(mean(rec2));
    avg_dE_white1 = getDeltaE_CIELAB(rec1_avg,rec2_avg,whiteRGB_1);
    
    dE_white1 = getDeltaE_CIELAB(rec1, rec2, whiteRGB_1);

    fov = sceneGet(s1, 'horizontal fov'); % in degree
    sampPerDeg = sceneGet(s1, 'cols') / fov;
    dE_SLAB = getDeltaE_SLAB(square1, square2, [1,1,1], sampPerDeg);

    %% Plot the two squares and the two reconstructed images

    f1 = figure('units','normalized','outerposition',[0 0 1 1]);
    subplot(2, 2, 1);
    imshow(square1);
    title('Input Squares');
    subplot(2, 2, 3);
    imshow(square2);
    subplot(2, 2, 2);
    imshow(rec1);
    title('Reconstruction')
    subplot(2, 2, 4);
    imshow(rec2);
    

    %% new figure: plot the two differences

    % find out maximum first
    maximum = max([max(dE_SLAB(:)), max(dE_white1(:))]);

    f2 = figure('units','normalized','outerposition',[0 0 1 1]);
    subplot(1,2,1);
    imshow(dE_SLAB, [0, maximum]);
    title('S-CIELAB (Scenes)');
    xlabel(['Mean: ', num2str(mean(dE_SLAB(:)))]);
    colorbar;
    subplot(1,2,2);
    imshow(dE_white1, [0, maximum]);
    title(['WP: ', mat2str(whiteRGB_1)]);
    xlabel(['Mean: ', num2str(avg_dE_white1)]);
    colorbar;


    %% Saves
    % create directory
    savedir = 'SquareExperiments/';
    expdir = strcat(savedir, exp.name, '/');
    mkdir(expdir);
    saveas(f1, strcat(expdir, 'images.jpg'));
    saveas(f1, strcat(expdir, 'images.fig'));
    saveas(f2, strcat(expdir, 'erorrs.jpg'));
    saveas(f2, strcat(expdir, 'errors.fig'));
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
    
    if max(imd(:)) ~= 0
        imd = imd / max(imd(:));
    end
    
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
function [s, cMosaic, square] = createSquareScene(col, square_size, fov)

    square = createSquare(col, square_size);

    % create scene from file
    [s, ~] = sceneFromFile(square, 'rgb');
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
    cMosaic.integrationTime = 1000;

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

function col2 = getCol2(col1)
    col2 = repmat(col1, 4, 1);
    if col1(1) == 0 col2(1, 1) = 0.05; else col2(1, 1) = 0.95; end
    if col1(2) == 0 col2(2, 2) = 0.05; else col2(2, 2) = 0.95; end
    if col1(3) == 0 col2(3, 3) = 0.05; else col2(3, 3) = 0.95; end
    
    col2(4, :) = [col2(1, 1), col2(2, 2), col2(3, 3)];
end

function dE_d = getDeltaE_SLAB(img1RGB, img2RGB, whiteRGB, sampPerDeg)
    img1_XYZ = rgb2xyz(img1RGB);
    img2_XYZ = rgb2xyz(img2RGB);
    whiteXYZ = rgb2xyz(whiteRGB);
    
    % get S-CIELAB Error image
    dE = scielab(sampPerDeg, img1_XYZ, img2_XYZ, whiteXYZ, 'xyz');
    dE_d = im2double(dE);
end
