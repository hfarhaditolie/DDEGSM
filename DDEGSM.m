function [score] = DDEGSM(Im1,Im2)
% Im1=im2double(rgb2ycbcr(Im1));
% Im2=im2double(rgb2ycbcr(Im2));

K1=Im1(:,:,2);
K2=Im1(:,:,3);
colorfulness_r = log(std2(K1)^2/(abs(mean2(K1))^0.2))*log((std2(K2)^2)/(abs(mean2(K2))^0.2));
[DEM_r,IC_r,LoGEM_r] = ref_Edge(Im1(:,:,1));
K1=Im2(:,:,2);
K2=Im2(:,:,3);
colorfulness_d = log(std2(K1)^2/(abs(mean2(K1))^0.2))*log((std2(K2)^2)/(abs(mean2(K2))^0.2));
RCS=(colorfulness_r)/(colorfulness_d);       
filterSize= [9,9];
kerx =[0 1 9 -1 0; 0 7 9 7 0; 9 9 0 -9 -9; 0 -7 -9 -7 0; 0 1 -9 -1 0];
% kerx=[2 2 2 2 2; 1 1 1 1 1; 0 0 0 0 0;-1 -1 -1 -1 -1;-2 -2 -2 -2 -2];

dirNum =12;
Im2Edge = 0;
dist = 180 / dirNum;
for ii = 0 : (dirNum-1)
    ker = (imrotate(kerx, ii*dist, 'bilinear', 'crop'));
    temp = (imfilter(Im2(:,:,1), ker,'replicate'));
    Im2Edge = (Im2Edge  + temp.^2);
end
DEM_d = sqrt(Im2Edge);
%Edge Strength Weight
W1 = fspecial('gaussian',filterSize, 2);
W2 = fspecial('gaussian',filterSize,2.1);
DOG_W = (W1 - W2);
W_d = abs(imfilter(Im2(:,:,1), DOG_W, 'replicate'));
ESW = max(IC_r, W_d);
% ESW=imfilter(ESW,fspecial('gaussian',[11 11],1.8));

% f1=floor(f1');
% m=size(f1,1);
% for i=1:m
%     ESW(f1(i,2),f1(i,1))=(ESW(f1(i,2),f1(i,1)))^.4;
% end
%%LOGEM Calculation

sigma=.85;
LoGKer = fspecial('log',filterSize,sigma);
LoGEM_d = (imfilter(Im2(:,:,1),LoGKer,'replicate'));
multiplyLog = LoGEM_r .* LoGEM_d;

%%ESM Calculation

T1=0.001;
LoGESM = (2*multiplyLog + T1) ./(LoGEM_r.^2 + LoGEM_d.^2 + T1);
T2 =2;
DESM = (2*DEM_r .* DEM_d + T2) ./(DEM_r.^2 + DEM_d.^2 + T2);

%%Pooling Strategy
ES = DESM.*LoGESM;
ES_ESW = ES.*ESW;
scoreY=(sum(ES_ESW(:))/sum(ESW(:)));

% scoreY=median(scores);
% return;
%%Distortion Detection
diff=abs(Im2(:,:,1)-Im1(:,:,1));
[h,w]=size(Im2(:,:,1));

d= (histcounts(diff,256));
% d2=d./(h*w);

% clss=predict(model,sqrt(d2));
% d= histcounts(diff,256)*(h*w);

threshold = compute_contrast(d,diff)/256;
pixel_not_eq=Im2(:,:,1)==Im1(:,:,1);
pixel_not_eq=size(find(pixel_not_eq==0),1)/(h*w) *100;
if pixel_not_eq <0.002
        scoreCSC=quality_evaluation_module(Im1(:,:,3),Im2(:,:,3));
        score = .2*RCS+.3*scoreCSC+.5*scoreY;
        disp('yes-csc')
%     score=scoreY;
elseif threshold>0.1953125
    disp('yes-cc')
    ker =ones(11,11);
    ker(6,6)=-1;
%     ker=ker./121;
    differenceDist= (imfilter(Im2(:,:,1),ker,'replicate'));
    differenceRef= (imfilter(Im1(:,:,1),ker,'replicate'));
    minMeans = min([mean2(differenceDist), mean2(differenceRef)]);
    maxMeans = max([mean2(differenceDist), mean2(differenceRef)]);
    B1=(minMeans)/(maxMeans);
    score=(threshold)*scoreY+(1-threshold)*B1;
else
    score = scoreY;
end
end
function score=quality_evaluation_module(ref_im,dist_im)
filterSize= [9 9];
kerx =[0 1 9 -1 0; 0 7 9 7 0; 9 9 0 -9 -9; 0 -7 -9 -7 0; 0 1 -9 -1 0];
dirNum =12;
Im2Edge = 0;
dist = 180 / dirNum;
for ii = 0 : (dirNum-1)
    ker = imrotate(kerx, ii*dist, 'bilinear', 'crop');
    temp = (imfilter(dist_im, ker,'replicate'));
    Im2Edge = (Im2Edge  + temp.^2);
end
DEM_d = sqrt(Im2Edge);

LoGKer = fspecial('log',filterSize,0.85);
LoGEM_d = (imfilter(dist_im,LoGKer,'replicate'));

[DEM_r, IC_r, LoGEM_r]= ref_Edge(ref_im);

T1=0.001;
multiplyLog = LoGEM_r .* LoGEM_d;
LoGESM = (2*multiplyLog + T1) ./(LoGEM_r.^2 + LoGEM_d.^2 + T1);
T2 =2;
DESM = (2*DEM_r .* DEM_d + T2) ./(DEM_r.^2 + DEM_d.^2 + T2);

W1 = fspecial('gaussian',filterSize, 2);
W2 = fspecial('gaussian',filterSize,2.1);
DOG_W = (W1 - W2);
W_d = abs(imfilter(dist_im, DOG_W, 'replicate'));
ESW = max(IC_r,W_d);

ES = DESM.*LoGESM;
ES_ESW = ES.*ESW;
score = (sum(ES_ESW(:))/sum(ESW(:)));
end
function contrast = compute_contrast(d,im)
contrast = 0;
[h,w]=size(im);
for i=1:256
    contrast=contrast+i*d(i);
end
contrast=contrast/(h*w);
end
function [Im1Edge1, Im1GausWeight, Im1EdgeLog]= ref_Edge(Im1)
kerx =[0 1 9 -1 0; 0 7 9 7 0; 9 9 0 -9 -9; 0 -7 -9 -7 0; 0 1 -9 -1 0];
filterSize= [9 9];

dirNum =12;
Im1Edge1 = 0;
dist = 180 / dirNum;
for ii = 0 : (dirNum-1)
    ker = imrotate(kerx, ii*dist, 'bilinear', 'crop');
    Im1Edge1 =  (Im1Edge1 + (imfilter(Im1, ker, 'replicate')) .^ 2);
end
Im1Edge1 = sqrt(Im1Edge1);
W1 = fspecial('gaussian',filterSize, 2);
W2 = fspecial('gaussian',filterSize,2.1);
DOG_W = (W1 - W2);
Im1GausWeight = abs(imfilter(Im1, DOG_W, 'replicate'));
logKer = fspecial('log',filterSize,0.85);
Im1EdgeLog = (imfilter(Im1,logKer,'replicate'));
end
