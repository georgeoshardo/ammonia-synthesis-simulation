exceldata='pipetemp.xlsx';
varydata=xlsread(exceldata);

X=varydata(:,1);
Y=varydata(:,2);
Z=varydata(:,3);
plot3(X,Y,Z,'.-')

tri = delaunay(X,Y);
plot(X,Y,'.')


[r,c] = size(tri);
disp(r)

h = trisurf(tri, X, Y, Z);
axis vis3d
shading interp
grid minor
ylabel('Temperature (K)');
xlabel('Pipe numer');
zlabel('Volume (m^3)');

lighting phong
material dull 


colorbar EastOutside
title('Optimising reactor volume');
colormap(jet)
box on