%% script creates a spoke or zoneplate target
% parameters
numRays = 25;
targetSize = 1024;
zoneScale = 1024;
target = 'spoke';

% define two RGB colors
color1 = [1, 1, 0]*0.7;
color2 = [0, 0, 1]*0.7;

image = zeros(targetSize, targetSize, 3);
c = targetSize / 2 + 0.5;


for i = 1:targetSize
    for j = 1:targetSize
        x = i - c;
        y = j - c;

        if strcmp(target, 'spoke')
            % calculate angle 
            phi = atan2(y, x);

            % fit numRays cosine waves into the 2*pi rad
            interpVal = 0.5 * cos(numRays*phi) + 0.5;
        end
        if strcmp(target, 'zoneplate')
            % calculate distance to center
            r = sqrt(x.^2 + y.^2);

            % get opacity value of Zoneplate
            interpVal = 0.5 * cos(zoneScale / targetSize.^2 * r.^2) + 0.5;
        end

        color = interpVal * color1 + (1-interpVal) * color2;

        image(i, j, :) = color;
    end
end


%% show
%imshow(imgaussfilt(image, 20));
figure();
imshow(image);