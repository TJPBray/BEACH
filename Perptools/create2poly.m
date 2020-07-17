%Creates polygon
function [xy1,xy2] = createpoly(polypoints,width,im)

%Preallocate coords
x1=[];
y1=[];
x2=[];
y2=[];

%Define flaring factor
fl=1.2;
%Define spacing factor (determines how far apart the lines are)
sp=15;

for l=2:(0.5*numel(polypoints)-3) %Range of l determines which drawn lines have perplines (The penultimate line doesn't have any perplines - these are provided by the last line)
%Get coords and length for each line
coords(:,1)=polypoints(l:l+1,1);
coords(:,2)=polypoints(l:l+1,2);

linelength=sqrt(((coords(1,1)-coords(2,1))^2+(coords(1,2)-(coords(2,2))^2)));

%Determine number of intersections needed
line=improfile(im,coords(:,1),coords(:,2));
numintersections=round(numel(line)/sp);
numintersections(numintersections<1)=1; %Minimum is one intersection at end of line

%Generate line for each intersection
for k=1:numintersections
xcoord(k,1)=coords(1,1) + ((coords(2,1)-coords(1,1))*(k/numintersections));
ycoord(k,1)=coords(1,2) + ((coords(2,2)-coords(1,2))*(k/numintersections));

grad(k,1)=(coords(2,2)-coords(1,2))/(coords(2,1)-coords(1,1));
perpgrad(k,1)=-1/grad(k,1);

xlength(k,1)=(width^2/(1+perpgrad(k,1)^2))^0.5; %Gives length in x and y direction of perplines
ylength(k,1)=perpgrad(k,1)*xlength(k,1);

%Record x and y coords of all perplines
x1(numel(x1)+1,1)=xcoord(k,1)+xlength(k,1);
y1(numel(y1)+1,1)=ycoord(k,1)+ylength(k,1);
x2(numel(x2)+1,1)=xcoord(k,1)-xlength(k,1);
y2(numel(y2)+1,1)=ycoord(k,1)-ylength(k,1);

x=[x1(numel(x1),1) x2(numel(x2),1)]; %Chooses current (last) x1 value
y=[y1(numel(y1),1) y2(numel(y2),1)];
end

end

%Generate start and end lines
%For first line
grad1=(polypoints(2,2)-polypoints(1,2))/(polypoints(2,1)-polypoints(1,1));
perpgrad1=-1/grad1;

xleng=((fl*width)^2/(1+perpgrad1^2))^0.5; %Gives length in x and y direction of perplines %fl is flaring factor
yleng=perpgrad1*xleng;

fx1=polypoints(2,1)+xleng;
fy1=polypoints(2,2)+yleng;

fx2=polypoints(2,1)-xleng;
fy2=polypoints(2,2)-yleng;

fx=[fx1 fx2];
fy=[fy1 fy2];

%For last line
n=numel(polypoints)*0.5;

gradl=(polypoints(n,2)-polypoints(n-1,2))/(polypoints(n,1)-polypoints(n-1,1));
perpgradl=-1/gradl;

xleng=((fl*width)^2/(1+perpgradl^2))^0.5; %Gives length in x and y direction of perplines
yleng=perpgradl*xleng;

lx1=polypoints(n-1,1)+xleng;
ly1=polypoints(n-1,2)+yleng;

lx2=polypoints(n-1,1)-xleng;
ly2=polypoints(n-1,2)-yleng;

lx=[lx1 lx2];
ly=[ly1 ly2];

%Combine coordinate data - creates x and y columns for each row of
%perpoints
%xy1 is one side of the polygon ('left')
xy1(1,1)=fx1;
xy1(1,2)=fy1;
xy1(2:numel(x1)+1,1)=x1;
xy1(2:numel(y1)+1,2)=y1;
xy1(numel(x1)+2,1)=lx1;
xy1(numel(y1)+2,2)=ly1;

xy2(1,1)=fx2;
xy2(1,2)=fy2;
xy2(2:numel(x2)+1,1)=x2;
xy2(2:numel(y2)+1,2)=y2;
xy2(numel(x2)+2,1)=lx2;
xy2(numel(y2)+2,2)=ly2;

%xy2 is the other side ('right')
xy2=flipud(xy2); %inverts xy2