function [a] = circ(A)
[i j k] = size(A);
a = zeros(i,j*k);
a(1:i,j+1:2*j) = A(:,:,2);
a(1:i,2*j+1:3*j) = A(:,:,3);
a(1:i,1:j) = A(:,:,1);
end
