%perptool4

im=double(imgcf);
hroi=impoly('Closed',false)
polypoints=hroi.getPosition;
delete(hroi);
width=8

%Get wide and narrow polygons
[xy1wide,xy2wide]=create2poly(polypoints,8,im)
[xy1narrow,xy2narrow]=create2poly(polypoints,2,im)

%Create left and right polygons
xyL=[xy1wide; flipud(xy1narrow)]
xyR=[xy2wide; flipud(xy2narrow)]

hroi=impoly(gca,xyL,'Closed',true);
hroi2=impoly(gca,xyR,'Closed',true);
setColor(hroi,'red');

        set(hroi,'Tag','impoly')
        set(hroi2,'Tag','impoly')