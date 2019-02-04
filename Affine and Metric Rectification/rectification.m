function H_a = assignment2(I,O1,O2)
% Part-1
% Read the image
I = imread(I);
imshow(I)

% Get 4 points such that they form 2 pairs of parallel lines
[xi,yi] = getline;

m1 = [xi(1) yi(1) 1];
m2 = [xi(2) yi(2) 1];
m3 = [xi(3) yi(3) 1];
m4 = [xi(4) yi(4) 1];

% Obtain the position of vanishing line
p_infinity1 = cross(cross(m1,m2),cross(m3,m4));
p_infinity2 = cross(cross(m1,m4),cross(m2,m3));
l_inf = cross(p_infinity1,p_infinity2);

% Affine Rectification
H_p = [ 1 0 0; 0 1 0; l_inf(1)/l_inf(3) l_inf(2)/l_inf(3) 1]
temp = maketform('projective', transpose(H_p));
I_a = imtransform(I,temp);

% Show and Save the image
myout = imshow(I_a);
saveas(myout,O1);


% Part-2
% Read the image
O1 = imread(O1);
imshow(O1);

% Get 8 points such that they form 2 pairs of orthogonal lines
[x,y] = getline;

l1 = cross([x(1) y(1) 1],[x(2) y(2) 1]);
m1 = cross([x(3) y(3) 1],[x(4) y(4) 1]);
l2 = cross([x(5) y(5) 1],[x(6) y(6) 1]);
m2 = cross([x(7) y(7) 1],[x(8) y(8) 1]);

% Solve Ax_ = B
A = [l1(1)*m1(1) l1(1)*m1(2) + l1(2)*m1(1); l2(1)*m2(1) l2(1)*m2(2) + l2(2)*m2(1)];
B = [-l1(2)*m1(2); -l2(2)*m2(2)];
x_ = linsolve(A,B);
S = [x_(1) x_(2); x_(2) 1];

% Get the value of K
% KK' = S
% K = U * sqrt(D) * U'
% V = U'
[U,D,V] = svd(S);
K = U * sqrt(D) * V;

% Form matrix H_a 
H_a = [K(1,1) K(1,2) 0; K(2,1) K(2,2) 0; 0 0 1];
if H_a(1,1) < 0
    H_a(1,1) = -H_a(1,1);

elseif H_a(2,2) < 0
    H_a(2,2) = -H_a(2,2);
end

% Metric Rectification
H = H_p * H_a;
temp = maketform('projective',transpose(H));
I_e = imtransform(I, temp);

% Show and Save the image
myout = imshow(I_e);
saveas(myout,O2);
return;