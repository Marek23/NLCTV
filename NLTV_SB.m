%% clear data on start
clear; close all; clc

%% load image
fp = imread('C:\MAREK\MAGISTERKA\Obrazy\flaga2.png');
fp = im2double(fp);

%% split of R, G and B
fp1 = fp(:,:,1); fp2 = fp(:,:,2); fp3 = fp(:,:,3);

%% finding the mask of image. Mask in image is defined as [R G B] =[0 1 0]
[nx, ny] = size(fp1);
lambdap = ones(size(fp2));
for x=1:nx
    for y=1:ny
        if(fp2(x,y) == 1 && fp3(x,y) ==0 && fp1(x,y)==0)
            lambdap(x,y) = 0;
        end
    end
end
lambda   = lambdap(:);
or       = find(lambdap >0);  %piksels of original image
m        = find(lambdap < 1); %piksels of mask 
[mx, my] = find(lambdap < 1); %piksels of mask in x,y coordinates
%% Initial values of image for algorithm
fp1(m) = fp(1,1); fp2(m) = fp(1,1); fp3(m) = fp(1,1);
f1     = fp1(:);  f2     = fp2(:);  f3     = fp3(:);
%% Initial image and mask
figure; imshow(fp);
figure; imshow(lambdap)

%% Params
SW = 2; %Search window size  
PS = 1; %Patch window size
%% Other params
iter_num = 20;
beta     = 70;
fp_num   = size(fp1,1)*size(fp1,2);%pixels amount
SWp_num = (2*SW +1)^2;%searchwindow pixels amount

%% Initialization of arguments for algorithm
uk1    = f1;  uk2    = f2;  uk3    = f3;
ukrsh1 = fp1; ukrsh2 = fp2; ukrsh3 = fp3;
bk1    = zeros(SWp_num,fp_num); bk2 = zeros(SWp_num,fp_num); bk3 = zeros(SWp_num,fp_num);
dk1    = zeros(SWp_num,fp_num); dk2 = zeros(SWp_num,fp_num); dk3 = zeros(SWp_num,fp_num);

aff_matrix1 = weight(fp1,lambdap,SW,PS);
aff_matrix2 = weight(fp2,lambdap,SW,PS);
aff_matrix3 = weight(fp3,lambdap,SW,PS);
%% Start of algorithm iterations
for krok = 1 : iter_num
    
    krok 
    %% finding b_k values
    bk1 = solveB(bk1,uk1,dk1,ukrsh1,aff_matrix1,fp_num,SWp_num,SW);
    bk2 = solveB(bk2,uk2,dk2,ukrsh2,aff_matrix2,fp_num,SWp_num,SW);
    bk3 = solveB(bk3,uk3,dk3,ukrsh3,aff_matrix3,fp_num,SWp_num,SW);
    
    %% finding u_k values
    uk1 = solveU(f1,ukrsh1,lambda,aff_matrix1,bk1,dk1,fp_num,SWp_num,beta,SW);
    uk2 = solveU(f2,ukrsh2,lambda,aff_matrix2,bk2,dk2,fp_num,SWp_num,beta,SW);
    uk3 = solveU(f3,ukrsh3,lambda,aff_matrix3,bk3,dk3,fp_num,SWp_num,beta,SW);
    
    %% actualization reshaped u
    ukrsh1 = uk1; 
    ukrsh2 = uk2; 
    ukrsh3 = uk3;
    
    %% finding BETA value
    [BETA1, BETA2, BETA3] = solveBETA(dk1,dk2,dk3,beta);

    %% finding d_k values
    dk1 = solveD(uk1,bk1,aff_matrix1,ukrsh1,fp_num,SWp_num,BETA1,beta,SW);
    dk2 = solveD(uk2,bk2,aff_matrix2,ukrsh2,fp_num,SWp_num,BETA2,beta,SW);
    dk3 = solveD(uk3,bk3,aff_matrix3,ukrsh3,fp_num,SWp_num,BETA3,beta,SW);
    
    %% reshape for uk
    ukrsh1 = reshape(uk1,nx,ny);
    ukrsh2 = reshape(uk2,nx,ny);
    ukrsh3 = reshape(uk3,nx,ny);
    
    %% Actualization of pixel values from mask in image fp
    fp1(m) = ukrsh1(m);
    fp2(m) = ukrsh2(m);
    fp3(m) = ukrsh3(m);

    %% Visualization of image at every iteration
    fwyn(:,:,1) = fp1;
    fwyn(:,:,2) = fp2;
    fwyn(:,:,3) = fp3;
    figure
    imshow(fwyn);
    
    %% Weight actualization
    aff_matrix1 = weightR(fp1,aff_matrix1,lambdap,mx,my,SW,PS);
    aff_matrix2 = weightR(fp2,aff_matrix2,lambdap,mx,my,SW,PS);
    aff_matrix3 = weightR(fp3,aff_matrix3,lambdap,mx,my,SW,PS);
end

%% Visualization of final image
fwyn(:,:,1) = fp1;
fwyn(:,:,2) = fp2;
fwyn(:,:,3) = fp3;
figure
imshow(fwyn);
