J=imcomplement(largestBlob);
imshow(J)
saveas(gcf,'area.jpg')
props=regionprops(largestBlob,'Area','Circularity','Eccentricity','Extent','MajorAxisLength','MinorAxisLength','Perimeter','Centroid','Orientation');

%Plot major axis length
x=props.Centroid(1);
y=props.Centroid(2);

deltax=props.MajorAxisLength*cosd(props.Orientation);
deltay=props.MajorAxisLength*sind(props.Orientation);

xVals=[x-deltax/2 x+deltax/2];
yVals=[y-deltay/2 y+deltay/2];
hold on
line(xVals,yVals)

%plot ellipse
phi = linspace(0,2*pi,50);
cosphi = cos(phi);
sinphi = sin(phi);


xbar = props.Centroid(1);
ybar = props.Centroid(2);

a = props.MajorAxisLength/2;
B = props.MinorAxisLength/2;

theta = pi*props.Orientation/180;
R = [ cos(theta)   sin(theta)
     -sin(theta)   cos(theta)];

xy = [a*cosphi; B*sinphi];
xy = R*xy;

x = xy(1,:) + xbar;
y = xy(2,:) + ybar;

plot(x,y,'r','LineWidth',2);

%plot convex hull perimeter
CH=bwconvhull(largestBlob); %convex hull.
CHP=bwperim(CH,8);
CHP=imcomplement(CHP);
imshow(CHP)

%plot cell perimeter
CP = bwperim(J,8);
CP=imcomplement(CP);
imshow(CP)