function [B1, B2, B3] = solveBETA(dk1,dk2,dk3,beta)
eps =0.0000001;
[SWpx_num, fpx_num ] = size(dk1);
BETA1 = zeros(1,fpx_num);
BETA2 = zeros(1,fpx_num);
BETA3 = zeros(1,fpx_num);
for j=1:fpx_num
    for i=1:SWpx_num
        BETA1(1,j) = BETA1(1,j) + dk1(i,j)^2;
        BETA2(1,j) = BETA2(1,j) + dk2(i,j)^2;
        BETA3(1,j) = BETA3(1,j) + dk3(i,j)^2;
    end
    BETA1(1,j) = sqrt(BETA1(1,j));
    BETA2(1,j) = sqrt(BETA2(1,j));
    BETA3(1,j) = sqrt(BETA3(1,j));
end

SBETA1 = sum(BETA1);
SBETA2 = sum(BETA2);
SBETA3 = sum(BETA3);

MBETA = sqrt(SBETA1^2 + SBETA2^2 + SBETA3^2) * beta;
B1 = SBETA1/(MBETA+eps);
B2 = SBETA2/(MBETA+eps);
B3 = SBETA3/(MBETA+eps);

end

