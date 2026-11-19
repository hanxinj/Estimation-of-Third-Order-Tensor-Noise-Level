function X = main0726(Sig)
t=[];X=[];
% Sig = 5;
Sigma = double(Sig);
K = 19;
% 训练集合中的y值
y =Sigma*Sigma*ones(1,K+1) + 0.02*randn(1,K+1);
% y =Sigma*Sigma*ones(1,K);
xlswrite('train_y_5.xls',y)
num = randperm(60,1);
for i = num:num+K
    imagename1 = strcat(num2str(i),'.jpg');
    t = main11(Sigma,imagename1);
    X = [X,t];
end
XX = X';
% 训练集合中的x值
xlswrite('train_x_5.xls',XX)

% 测试集合中的x值
for j = 1:12
    imagename = strcat(num2str(j),'.bmp');
test_X=[];
for i = 1:3
    t = main11(Sigma,imagename);
    test_X = [test_X,t];
end
test_set = test_X';
xlsname = strcat('test',num2str(j),'_x_5.xls');
xlswrite(xlsname,test_set)
end
end

function lambda = main11(Sigma,pic_name)
a = imread(pic_name);
if size(a,3) > 1
    a = a(:,:,1);
end
na = double(a)+Sigma*randn(size(a));
[nlevel,d] =  NoiseLevel0518(na,7); 
lambda = d(49-7:49);
end
