function [output] = solveU(f,ukrsh,lambda,aff_matrix,bk,dk,fpx_num,SWp_num,beta,SW)
%% u
lampar = 0.01;
output = f;
[nx, ny] = size(ukrsh);

for iter_u = 1:2
    
    for  i = 1: length(lambda)

        if(lambda(i) ~= 0)
            output(i) = f(i);
        else
            sumw = sum(aff_matrix(:,i));
            
            fDen = (lampar+beta*sumw);
            
            [pom, XP, YP] = searchWindow(i,ukrsh,SW);
            aff_tmp = aff_matrix(:,i).*pom;
            sumwu = sum(aff_tmp);
            
            sumdb = 0;
            for j=1:SWp_num
                
                yp = floor((j-1)/(SW*2+1))+1;
                xp = rem(j-1,(SW*2+1))+1;
                
                %% pocz¹tek uk³adu SW
                PXP = XP-2*SW;
                PYP = YP-2*SW;
                
                hx = PXP-1+xp;
                hy = PYP-1+yp;
                
                yr = 2*SW+1 +1 - yp;
                xr = 2*SW+1 +1 - xp;
                hp = xr+(2*SW+1)*(yr-1);
                
                jl = hx+(nx)*(hy-1);
                
                sumdb = sumdb + sqrt(aff_matrix(j,i))*(dk(j,i) - dk(hp,jl) - bk(j,i) + bk(hp,jl));
                
            end
            
            output(i) = (sumwu*beta + lampar*f(i) - beta*sumdb)/(fDen+eps);
        end

    end
end
end

