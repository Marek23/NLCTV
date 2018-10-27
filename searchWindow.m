function [output, XP, YP, yp, xp] = searchWindow(i,im,SW)

[nx, ny] = size(im);
yp = floor((i-1)/nx)+1+SW;
ody = yp-SW;
doy = yp+SW;

YP = yp;
yp = floor((i-1)/nx)+1+SW;

xp = rem(i-1,nx)+1+SW;
odx = xp-SW;
dox = xp+SW;

XP = xp;
xp = rem(i-1,nx)+1+SW;

imSW = padarray(im,[SW, SW],'symmetric');
pom = imSW(odx:dox,ody:doy);
output = pom(:);

end