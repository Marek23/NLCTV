function [output] = weight(im,mask,SW,f)

sp_num = size(im,1)*size(im,2);
[nx, ny]=size(im);

SWp_num = (2*SW +1)^2;
sigma = 0.0001; %param r from page 857 in article
%obraz z z ³ódk¹ ma 0.4
%dla bungee 0.06!
output=zeros(SWp_num,sp_num);

imSWf = padarray(im,[f+SW, f+SW],'symmetric');
maskSWf = padarray(mask,[f+SW f+SW],'symmetric');
%% 
% Used kernel
kernel = make_kernel(f);
kernel = kernel / sum(sum(kernel));
%kernel = fspecial('gaussian',2*f+1,1);


%%
ITER=1;
for l=1:ny
    Y = l + f + SW;
    
    for k=1:nx
        
        X = k + f + SW;
        
        Fim = imSWf(X-f:X+f , Y-f:Y+f);
        chi = maskSWf(X-f:X+f , Y-f:Y+f);
        
        iter = 1;
        for y=Y-SW:Y+SW
            
            for x=X-SW:X+SW
                
                fim = imSWf(x-f:x+f , y-f:y+f);
                
                d = sum(sum(chi.*(kernel.*(Fim-fim).*(Fim-fim))));
                
                w = exp(-d/sigma^2);
                
                
                output(iter,ITER) = w;
                iter = iter+1;
            end
        end
        ITER = ITER+1;
    end
end
output = output/max(output(:));
end
