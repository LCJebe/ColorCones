function image = createSpoke(color1, color2, numRays, targetSize)
% function creates a spoke target

image = zeros(targetSize, targetSize, 3);
c = targetSize / 2 + 0.5;

for i = 1:targetSize
    for j = 1:targetSize
        x = i - c;
        y = j - c;

        % calculate angle 
        phi = atan2(y, x);

        % fit numRays cosine waves into the 2*pi rad
        interpVal = 0.5 * cos(numRays*phi) + 0.5;

        color = interpVal * color1 + (1-interpVal) * color2;

        image(i, j, :) = color;
    end
end

end