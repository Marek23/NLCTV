function [output] = solveU(f,ukrsh,lambda,aff_matrix,bk,dk,fpx_num,SWpx_num,SW,beta)
%% u
lampar = 1;
output = f;
MfDen = zeros(size(f));
Msumwu = zeros(size(f));
Msumdb = zeros(size(f));

for iter_u = 1:2
    
    for  i = 1: length(lambda)
        
        pom = lambda(i);
        sumw = sum(aff_matrix(:,i));

        fDen = 1/(lampar*lambda(i)+beta*sumw);

        
        pom = searchWindow(i,ukrsh,SW);
        aff_tmp = aff_matrix(:,i).*pom;
        sumwu = sum(aff_tmp);
        
        sumdb = 0;
        for j=1:SWpx_num
            
            sumdb = sumdb + sqrt(aff_matrix(j,i))*(2*dk(j,i) -...
                2*bk(j,i));
            
        end
        MfDen(i) = fDen;
        Msumwu(i) = sumwu;
        Msumdb(i) =  sumdb;
        output(i) = fDen*(sumwu*beta + lampar*lambda(i)*f(i) - beta/2*sumdb);

        
    end
end
end

