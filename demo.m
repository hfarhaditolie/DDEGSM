
%==========================================================================
%Paper: ( H. Farhadi Tolie and M. R. Faraji, "Screen Content Image Quality Assessment Using Distortion-based Directional Edge
%and Gradient Similarity Maps", Manuscript submitted to Elsevier’s journal of Signal Processing: Image Communication, 2021)

%E-mail Addresses: (h.farhadi@iasbs.ac.ir, h.farhaditolie@gmail.com), (m.faraji@iasbs.ac.ir, mohammadreza.faraji@aggiemail.usu.edu)
%==========================================================================

clear;
clc;

im1 = imread('im1.bmp');
im2 = imread('im2.bmp');

score = DDEGSM(im1,im2);
score
