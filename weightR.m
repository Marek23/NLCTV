function [output] = weightR(im,weight,mask,mx,my,SW,f)

%% poczatkowa inicjalizacja outputu
output = weight;
[nx, ny]=size(im);

sigma = 0.0001;



imSWf = padarray(im,[f+SW, f+SW],'symmetric');
maskSWf = padarray(mask,[f+SW f+SW],'symmetric');

%% Used kernel
kernel = make_kernel(f);
kernel = kernel / sum(sum(kernel));

%%
for pkt=1:length(mx)
    
    YP = my(pkt);
    XP = mx(pkt);
    ITER = (YP-1)*nx + XP;
    Y = YP+f+SW;
    X = XP+f+SW;
    
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
end
end
