function [ im ] = imageFromConeAbs( cMosaic , search_radii)
% Takes a cMosaic and creates an image based on the cone mosaic and 
% its cone absorptions. 
% Returns im (RGB image)

% get cone absorptions
cone_abs = cMosaic.absorptions;
single_abs = cone_abs(:,:,1); 

% get mask for L, M, S cones
p = cMosaic.pattern;
L_cones = p==2;
M_cones = p==3;
S_cones = p==4;

% get absorptions for cones separately
L_abs = L_cones .* single_abs;
M_abs = M_cones .* single_abs;
S_abs = S_cones .* single_abs;

% stack to 3D matrix
[x,y] = size(L_abs);

LMS_abs = zeros(x,y,3);
LMS_abs(:,:,1) = L_abs;
LMS_abs(:,:,2) = M_abs;
LMS_abs(:,:,3) = S_abs;

LMS_cones = zeros(x,y,3);
LMS_cones(:,:,1) = L_cones;
LMS_cones(:,:,2) = M_cones;
LMS_cones(:,:,3) = S_cones;

for k = 1:3
    absorptions = LMS_abs(:,:,k);
    cone_layer = LMS_cones(:,:,k);
    th = search_radii(k);
    
    temp_layer = absorptions;
    
    for i = 1:x
        for j = 1:y
            % continue, if we already have the cone absoprtion
            curr = absorptions(i,j);
            if curr ~= 0
                continue
            end
            
            % find nonzero coords in bounding box
            mean_absorption = getResponseEstimate(cone_layer, absorptions, i, j, x, y, th);

            temp_layer(i,j) = mean_absorption;
            
%             Next Neighbor Approach
%             [M,ind] = min(dists);
%             neighbor_coords = nonzero_coords(ind,:); %average 
%             
%             n_row = neighbor_coords(1,1);
%             n_col = neighbor_coords(1,2);
%             temp_layer(i,j) = layer(n_row,n_col); 
        end
    end
    
    
    LMS_abs(:,:,k) = temp_layer;
end



im = zeros(x,y,3);
im(:,:,1) = LMS_abs(:,:,1);
im(:,:,2) = LMS_abs(:,:,2);
im(:,:,3) = LMS_abs(:,:,3);

end


function mean_absorption = getResponseEstimate(mask, absorptions, i, j, x, y, th)
    % i, j is the pixel we're looking at
    % x, y is the image size
    % th is the boundig box radius
    i_min = max(1, i-th);
    i_max = min(x, i+th);
    j_min = max(1, j-th);
    j_max = min(y, j+th);
    
    % get nonzero coordinates within boundig box
    maskBox = mask(i_min:i_max, j_min:j_max);
    absorptionsBox = absorptions(i_min:i_max, j_min:j_max);
    [row, col] = find(maskBox);
    nonzero_coords = [row col];
    nonzero_idx = sub2ind(size(maskBox), row, col);
    
    % only retain distances within certain L2 radius
    center = [i-i_min, j-j_min]+1;
    diff = center - nonzero_coords;
    dists = sqrt(sum(diff.^2,2));
    
    small_idx = nonzero_idx(dists < th);
    if ~isempty(small_idx)
        mean_absorption = mean(absorptionsBox(small_idx));
    else % recurse with a larger search radius
        mean_absorption = getResponseEstimate(mask, absorptions, i, j, x, y, th+1);
    end
    
end

