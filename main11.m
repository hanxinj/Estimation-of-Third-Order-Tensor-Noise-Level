% function [lambda,nlevel] = main11(Sigma,pic_name)
 addpath(genpath('./train/'));
addpath(genpath('./LIVE/'));
addpath(genpath('./test/'));
addpath(genpath('./tid2008_images/'));
Sigma = 20;

a = imread('I01.BMP');
% a = pic_name;
if size(a,3) > 1
    a = circ(a);
%  a = a(:,:,1);
end
na = double(a)+Sigma*randn(size(a));
[nlevel,d] = NoiseLevel0518(na,7);
lambda = d(1:8)
nlevel
% end