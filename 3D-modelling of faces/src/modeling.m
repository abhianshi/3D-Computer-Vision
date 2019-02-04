function modeling(im)
[P1,P2,H,Hr,pnts] = extcaoforooshcalib(im);
% % Read the first image from the image set.
% I = imread(im);
% % Change to grayScale if image is RGB
% if size(I,3)==3
%     grayImg = rgb2gray(I);
% end
% % Detect and Extract Features for start Image
% points = detectSURFFeatures(grayImg);
% xy = points(:).Location;
% xy = xy';

% save('face5.mat','pnts');
XY=load('face5.mat');
xy = XY.pnts;
[n_rows,n_columns] = size(xy);
% xy_ = [xy; ones(1,n_columns)];
% pnts1 = hnormalise(Hr*xy_);
pnts1 = hnormalise(Hr*xy);

z = [-1 0 0; 0 1 0; 0 0 1];
pnts2= hnormalise(Hr*z*hnormalise(H*pnts1));

p3d = [];

for i = 1:n_columns
 A = triangulation(P1,P2,pnts1(:,i),pnts2(:,i));
 [U,D,V] = svd(A);
 x_w = V(:,end);
 x_w = x_w./x_w(4);
 X_wn = [x_w(1); x_w(2); x_w(3)];
 p3d = [p3d X_wn];
end


figure; 
scatter3(p3d(1,:), p3d(2,:), p3d(3,:)); 
axis vis3d;

end

function [A] = triangulation(P1,P2,x1,x2)
A = [   x1(1)*P1(3,:) - P1(1,:);
        x1(2)*P1(3,:) - P1(2,:);
        x2(1)*P2(3,:) - P2(1,:);
        x2(2)*P2(3,:) - P2(2,:) ];
end
