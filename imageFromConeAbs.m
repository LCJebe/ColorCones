function [ im ] = imageFromConeAbs( cMosaic )
% Takes a cMosaic and creates an image based on the cone mosaic and 
% its cone absorptions. 
% Returns im (RGB image)

cone_abs = cMosaic.absorptions;
single_abs = cone_abs(:,:,1); 

p = cMosaic.pattern;
L_cones = p==2;
M_cones = p==3;
S_cones = p==4;

L_abs = L_cones .* single_abs;
M_abs = M_cones .* single_abs;
S_abs = S_cones .* single_abs;

[x,y] = size(L_abs);

LMS_abs = zeros(x,y,3);
LMS_abs(:,:,1) = L_abs;
LMS_abs(:,:,2) = M_abs;
LMS_abs(:,:,3) = S_abs;

LMS_cones = zeros(x,y,3);
LMS_cones(:,:,1) = L_cones;
LMS_cones(:,:,2) = M_cones;
LMS_cones(:,:,3) = S_cones;

thresholds = [4;5;10;];

for k = 1:3
    layer = LMS_abs(:,:,k);
    cone_layer = LMS_cones(:,:,k);
    
    temp_layer = layer;
    [row,col] = find(cone_layer); 
    
    
    nonzero_coords = [row col];
    
    for i = 1:x
        for j = 1:y
            curr = layer(i,j);
            if curr ~= 0
                continue
            end
            coords = [i j];
            diff = coords - nonzero_coords;
            
            dists = sqrt(sum(diff.^2,2));
            
            M = min(dists);
            threshold = thresholds(k,1);
            small_dists = dists<threshold;
            n_small_dists = nnz(small_dists);
            sum_vals = 0;
            for l = 1:numel(small_dists)
                if small_dists(l) == 0
                    continue
                end 
                ind = l;
                neighbor_coords = nonzero_coords(ind,:);
                n_row = neighbor_coords(1,1);
                n_col = neighbor_coords(1,2);
                sum_vals = sum_vals + layer(n_row,n_col);
            end
            
            sum_vals = sum_vals / n_small_dists;
            temp_layer(i,j) = sum_vals;
            
            
            
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

