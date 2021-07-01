function [score] = DDEGSM(Im1,Im2)
filterSize= [9,9];
%Color Space Conversion
Im1=im2double(rgb2ycbcr(Im1));
Im2=im2double(rgb2ycbcr(Im2));

%Computing Directional Edge Map (DEM),LoG Edge Map (LoGEM), and Image Contour (IC) of the reference image
[DEM_r,IC_r,LoGEM_r] = ref_Edge(Im1(:,:,1));

%Computing Directional Edge Map (DEM) of the distorted image
kerx =[0 1 9 -1 0; 0 7 9 7 0; 9 9 0 -9 -9; 0 -7 -9 -7 0; 0 1 -9 -1 0];
dirNum =12;
Im2Edge = 0;
dist = 180 / dirNum;
for ii = 0 : (dirNum-1)
    ker = (imrotate(kerx, ii*dist, 'bilinear', 'crop'));
    temp = (imfilter(Im2(:,:,1), ker,'replicate'));
    Im2Edge = (Im2Edge  + temp.^2);
end
DEM_d = sqrt(Im2Edge);

%Computing IC of the distorted image
W1 = fspecial('gaussian',filterSize, 2);
W2 = fspecial('gaussian',filterSize,2.1);
DOG_W = (W1 - W2);
IC_d = abs(imfilter(Im2(:,:,1), DOG_W, 'replicate'));
ESW = max(IC_r, IC_d);

%Computing LoGEM of the distorted image
sigma=.85;
LoGKer = fspecial('log',filterSize,sigma);
LoGEM_d = (imfilter(Im2(:,:,1),LoGKer,'replicate'));

%Computing the Edge Similarity Maps
T1=0.001;
LoGESM = (2*LoGEM_r .* LoGEM_d + T1) ./(LoGEM_r.^2 + LoGEM_d.^2 + T1);
T2 =2;
DESM = (2*DEM_r .* DEM_d + T2) ./(DEM_r.^2 + DEM_d.^2 + T2);

%Pooling Strategy
ES = DESM.*LoGESM;
ES_ESW = ES.*ESW;
scoreY=(sum(ES_ESW(:))/sum(ESW(:)));

%Computing NWF and tr
DiffIm=abs(Im2(:,:,1)-Im1(:,:,1));
img_histogram= (histcounts(DiffIm,256));
tr = NWF(img_histogram,DiffIm)/256;

%Computing RNNEP
pixel_not_eq=Im2(:,:,1)==Im1(:,:,1);
[h,w]=size(Im2(:,:,1));
RNNEP=size(find(pixel_not_eq==0),1)/(h*w) *100;

%CSC Distortion Detection
if RNNEP <0.002
    K1=Im1(:,:,2);
    K2=Im1(:,:,3);
    colorfulness_r = log(std2(K1)^2/(abs(mean2(K1))^0.2))*log((std2(K2)^2)/(abs(mean2(K2))^0.2));
    K1=Im2(:,:,2);
    K2=Im2(:,:,3);
    colorfulness_d = log(std2(K1)^2/(abs(mean2(K1))^0.2))*log((std2(K2)^2)/(abs(mean2(K2))^0.2));
    
    
    RCS= (colorfulness_r)/(colorfulness_d);
    score_Cr=quality_evaluation_module(Im1(:,:,3),Im2(:,:,3));
    score_CSC = .2*RCS+.3*score_Cr+.5*scoreY;
    score=score_CSC;
%CC Distortion Detection
elseif tr>0.1953125
    ker =ones(11,11);
    ker(6,6)=-1; 
    LocalMean_d= mean2(imfilter(Im2(:,:,1),ker,'replicate'));
    LocalMean_r= mean2(imfilter(Im1(:,:,1),ker,'replicate'));
    minMeans = min([LocalMean_d, LocalMean_r]);
    maxMeans = max([LocalMean_d, LocalMean_r]);
    ro=(minMeans)/(maxMeans);
    
    score_CC=(tr)*scoreY+(1-tr)*ro;
    score=score_CC;
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
function contrast = NWF(d,im)
contrast = 0;
[h,w]=size(im);
for i=1:256
    contrast=contrast+i*d(i);
end
contrast=contrast/(h*w);
end
function [DEM, IC, LoGEM]= ref_Edge(Im1)
filterSize= [9 9];

%Computing DEM
kerx =[0 1 9 -1 0; 0 7 9 7 0; 9 9 0 -9 -9; 0 -7 -9 -7 0; 0 1 -9 -1 0];
dirNum =12;
DEM = 0;
dist = 180 / dirNum;
for ii = 0 : (dirNum-1)
    ker = imrotate(kerx, ii*dist, 'bilinear', 'crop');
    DEM =  (DEM + (imfilter(Im1, ker, 'replicate')) .^ 2);
end
DEM = sqrt(DEM);

%Computing IC
W1 = fspecial('gaussian',filterSize, 2);
W2 = fspecial('gaussian',filterSize,2.1);
DOG_W = (W1 - W2);
IC = abs(imfilter(Im1, DOG_W, 'replicate'));

%Computing LoGEM
logKer = fspecial('log',filterSize,0.85);
LoGEM = (imfilter(Im1,logKer,'replicate'));
end
