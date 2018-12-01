function image = createZoneplate(color1, color2, zoneScale, targetSize)
%% function creates a zoneplate target

image = zeros(targetSize, targetSize, 3);
c = targetSize / 2 + 0.5;

for i = 1:targetSize
    for j = 1:targetSize
        x = i - c;
        y = j - c;

        % calculate distance to center
        r = sqrt(x.^2 + y.^2);

        % get opacity value of Zoneplate
        interpVal = 0.5 * cos(zoneScale / targetSize.^2 * r.^2) + 0.5;
        
        color = interpVal * color1 + (1-interpVal) * color2;

        image(i, j, :) = color;
    end
end


%% show
%imshow(imgaussfilt(image, 20));
figure();
imshow(image);