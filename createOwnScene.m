clear all; close all;
%% create our own scene form a colorful spoke target

% create spoke
col1 = [1, 1, 0];
col2 = [0, 0, 1];
spoke = createSpoke(col1, col2, 25, 512);

% create scene from file
[s, ~] = sceneFromFile(spoke, 'rgb');

%%

% set some scene parameters
s = sceneSet(s, 'fov', 1);

% add object to scene and show scene window
ieAddObject(s);
sceneWindow;