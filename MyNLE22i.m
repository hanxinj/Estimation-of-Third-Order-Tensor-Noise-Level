function [nlevel,d] = MyNLE22i(img,w,KK,decim,conf,itr)
% nlevel=[];
if( ~exist('itr', 'var') )
    itr = 3;
end

if( ~exist('conf', 'var') )
    conf = 1-1E-6;
end

if( ~exist('decim', 'var') )
    decim = 0;
end

if( ~exist('patchsize', 'var') )
    patchsize = w;
end
kh = [-1/2,0,1/2];
imgh = imfilter(img,kh,'replicate');
imgh = imgh(:,2:size(imgh,2)-1,:);
imgh = imgh .* imgh;

kv = kh';
imgv = imfilter(img,kv,'replicate');
imgv = imgv(2:size(imgv,1)-1,:,:);
imgv = imgv .* imgv;

Dh = my_convmtx2(kh,patchsize,patchsize);
Dv = my_convmtx2(kv,patchsize,patchsize);
DD = Dh'*Dh+Dv'*Dv;
r = rank(DD);

Dtr = trace(DD);
tau0 = gaminv(conf,double(r)/2, 2 * Dtr / double(r));

n=patchsize*patchsize;

liu = randperm(200,1);

for cha=1:size(img,3)
    X = im2col(img(liu:liu+25,liu:liu+25,cha),[patchsize patchsize]);    %%以子块滑动的方式将图像块分为m*n的子矩阵，沿列的方向装换成B的列
    Xh = im2col(imgh(liu:liu+25,liu:liu+23,cha),[patchsize patchsize-2]);
    Xv = im2col(imgv(liu:liu+23,liu:liu+25,cha),[patchsize-2 patchsize]);
    Xtr = sum(vertcat(Xh,Xv));
    
    if( decim > 0 )
        XtrX = vertcat(Xtr,X);%%垂直串联数组
        XtrX = sortrows(XtrX')';  %%列方向升序排列     
        p = floor(size(XtrX,2)/(decim+1));%%不大于X的最大整数
        p = [1:p] * (decim+1);
        Xtr = XtrX(1,p);
        X = XtrX(2:size(XtrX,1),p);
    end
    
    %%%%% noise level estimation %%%%%
    if( size(X,2) < size(X,1) )
        sigw = 0;
    else
        cov = X*X'/(size(X,2)-1);
        d = flipud(eig(cov));%%flipud函数实现矩阵的上下翻转。如：flipud(QQ)实现了矩阵QQ的上下翻转。
                             %%注意：区别于fliplr函数，该函数是实现矩阵的左右翻转。 
        sigw = d(n);
    end 
    
    for i = 2:itr        
        %%%% weak texture selectioin %%%%%
        tau = sigw * tau0;
        p = (Xtr<tau);
        Xtr = Xtr(:,p);
        X = X(:,p);
        if( size(X,2) < size(X,1) )            
            break;
        end
%         X = bsxfun(@minus,X,mean(X,2));
        cov = X*X'/(size(X,2)-1);
        d = flipud(eig(cov));
        
        CompCnt = NoiseEstimationd(d);%主维数  2
        Pd = n-CompCnt;%冗余维数    47 
        
        alpha = 0.375;
        i1=floor(Pd); % 最后一个，向下取整
        i2 = Pd-KK;   % 47-3 = 44
        d1 = icdf('norm',(Pd-alpha+1-i1)/(Pd-2*alpha+1));
        d2 = icdf('norm',(Pd-alpha+1-i2)/(Pd-2*alpha+1));
        
        sig1 = d(i1);
        sig2 = d(CompCnt+i2);  % 2+44 = 46 %% = n-kk
        sigw = (d1*sig2-d2*sig1)/(d1-d2);
        if sqrt(abs(sigw)) > 35            
            KK = 3;
            i2 = Pd-KK;
            d2 = icdf('norm',(Pd-alpha+1-i2)/(Pd-2*alpha+1));
            sig2 = d(CompCnt+i2);
            sigw = (d1*sig2-d2*sig1)/(d1-d2);
        end
        
        
    end
    nlevel(cha)=sqrt(abs(sigw));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function T = my_convmtx2(H, m, n)
s = size(H);
T = zeros((m-s(1)+1) * (n-s(2)+1), m*n);

k = 1;
for i=1:(m-s(1)+1)
    for j=1:(n-s(2)+1)
        
        for p=1:s(1)
            T(k,(i-1+p-1)*n+(j-1)+1:(i-1+p-1)*n+(j-1)+1+s(2)-1) = H(p,:);
        end
        
        k = k + 1;
    end
end

function CompCnt =NoiseEstimationd(eigvalue)
% p_out = image2cols(im_noisy, PatchSize, 3);

% mu = mean(p_out,2);
% sigma=(p_out-repmat(mu,[1,size(p_out,2)])) ...
%         *(p_out-repmat(mu,[1,size(p_out,2)]))'/(size(p_out,2));
% eigvalue = (sort((eig(sigma)),'ascend')); 
for CompCnt = size(eigvalue,1):-1:1
    Mean = mean(eigvalue(1:CompCnt));    
    if(sum(eigvalue(1:CompCnt)>Mean) == sum(eigvalue(1:CompCnt)<Mean))
        break
    end   
end
% estsigma = sqrt(Mean);
