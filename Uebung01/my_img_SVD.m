% function [] = MyFun()
%
%   inputs:
%       imgPath in string,		path to gray scale image
%		k in {1; rank(img)},	rank of approximation
%
%   outputs:
%       
% 
%   MyFun(..) ...


function [] = my_img_SVD(imgPath, k)
    fprintf('-------- my_img_SVD() --------\n');
    A = imread(imgPath);
	A = im2double(A);
    
	% create ||..||_2-optimal rank k approximation of A using SVD
	[U, S, V] = svd(A);
	B =   U(:          , (1 : 1 : k))...
		* S((1 : 1 : k), (1 : 1 : k))...
		* V(:          , (1 : 1 : k))';
	
	% plot original and approximation
	figure('WindowStyle', 'docked');
    subplot(1, 2, 1);
	imshow(A, 'InitialMagnification', 'fit');
    title('Original (A)');
	
	subplot(1, 2, 2);
    imshow(B, 'InitialMagnification', 'fit');
    title(sprintf('SVD (k = %3i), ||A - B||_2=%.3e',...
		k, S(k+1, k+1))...
	);
	
    fprintf('------------------------------\n\n');
end