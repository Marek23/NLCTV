function [output] = solveD(uk,bk,aff_matrix,ukrsh,fpx_num,SWpx_num,BETA,beta,SW)
%% d
output = zeros(size(bk));

for i=1:fpx_num
    
    max_value_tmp =0;
    ukp = searchWindow(i,ukrsh,SW);
    
    for j=1:SWpx_num
        
%         max_value_tmp = max_value_tmp + (sqrt(aff_matrix(j,i))*(ukp(j)...
%             - uk(i))+ bk(j,i))^2;

        max_value_tmp = max_value_tmp + (sqrt(aff_matrix(j,i))*(ukp(j)...
            - uk(i)) + bk(j,i))^2;
    end
    
    for j=1:SWpx_num
        
        sumf2 = sqrt(max_value_tmp);
        max_value = sumf2 - BETA;
%         max_value = sumf2 - beta;
        max_value = max(0,max_value);
        
        
        sumf1 = sqrt(aff_matrix(j,i))*(ukp(j) - uk(i)) + bk(j,i);
        
        output(j,i) = sumf1/(sumf2+eps)*max_value;
%         if sumf2 ~= 0
%             output(j,i) = sumf1/sumf2*max_value;
%         end
        
    end
end
end

