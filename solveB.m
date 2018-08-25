function [output] = solveB(bk,uk,dk,ukrsh,aff_matrix,fpx_num,SWpx_num,SearchWindow)
    %% b
    output = bk;
    for i=1:fpx_num
        ukp = searchWindow(i,ukrsh,SearchWindow);
        for j=1:SWpx_num
            output(j,i) = bk(j,i) + sqrt(aff_matrix(j,i))*(ukp(j) - uk(i)) - dk(j,i);
        end
    end
end

