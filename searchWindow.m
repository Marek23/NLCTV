function [output] = searchWindow(i,im,SW)

[nx, ny] = size(im);
yp = floor((i-1)/nx)+1+SW;
ody = yp-SW;
doy = yp+SW;

xp = rem(i-1,nx)+1+SW;
odx = xp-SW;
dox = xp+SW;

imSW = padarray(im,[SW, SW],'symmetric');
pom = imSW(odx:dox,ody:doy);
output = pom(:);

end