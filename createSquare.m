function [ image ] = createSquare( color, targetSize )
%creates a uniform color square

image = zeros(targetSize, targetSize, 3);

for i = 1:targetSize
    for j = 1:targetSize
        image(i,j,:) = color;
    end
end

end

