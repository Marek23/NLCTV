function [output] = weight(im,mask,SW,f)

sp_num = size(im,1)*size(im,2);
[nx, ny]=size(im);
fpom = f;

SWp_num = (2*SW +1)^2;
sigma = 0.01; %param r from page 857 in article
%obraz z z ³ódk¹ ma 0.4
%dla bungee 0.06!
output=zeros(SWp_num,sp_num);

imSWf   = padarray(im,[f+SW, f+SW],'symmetric');
imshow(imSWf);
maskSWf = padarray(mask,[f+SW f+SW],'symmetric');
%% making kernel
sigmaKennel  = 5; 
kernel = fspecial('gaussian',2*f+1,sigmaKennel);

%%
ITER=1;
for l=1:ny
    Y = l + f + SW;
    
    for k=1:nx
        
        X = k + f + SW;
        
        clear Fim
        clear chil
        clear kernel
        
        if maskSWf(X,Y) == 2
            f = 17;
            kernel = fspecial('gaussian',2*f+1,sigmaKennel);
        else
            f = fpom;
            kernel = fspecial('gaussian',2*f+1,sigmaKennel);
        end
        Fim = imSWf(X-f:X+f , Y-f:Y+f);
        chil = maskSWf(X-f:X+f , Y-f:Y+f);
        
        iter = 1;
        clear d
        clear fim
        for y=Y-SW:Y+SW
            
            for x=X-SW:X+SW
                
                fim = imSWf(x-f:x+f , y-f:y+f);
                
                d = -sqrt(sum(sum(chil.*(kernel.*(Fim-fim).*(Fim-fim)))));
                
                %w = exp(d/sigma^2);
                w = exp(d*10);
                
                output(iter,ITER) = w;
                iter = iter+1;
            end
        end
%         output(:,ITER) = output(:,ITER)/sum(output(:,ITER));
        ITER = ITER+1;
    end
end
% sumOutput = sum(output(:));
% output = output/sumOutput;
end
