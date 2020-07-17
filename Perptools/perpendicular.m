

%Get coords and length of defined line
coords=roi.getPosition
line=improfile(im,coords(:,1),coords(:,2))
linelength=((coords(1,1)-coords(2,1))^2+(coords(1,2)-(coords(2,2))^2))^0.5

%Determine number of intersections needed
numintersections=abs(round(linelength/5))

%Generate line for each intersection
for k=1:numintersections-1
xcoord(k,1)=coords(1,1) + (coords(2,1)-coords(1,1))*(k/numintersections)
ycoord(k,1)=coords(1,2) + (coords(2,2)-coords(1,2))*(k/numintersections)

grad(k,1)=(coords(2,2)-coords(1,2))/(coords(2,1)-coords(1,1))
perpgrad(k,1)=-1/grad(k,1)

xlength(k,1)=(length^2/(1+perpgrad(k,1)^2))^0.5 %Gives length in x and y direction of added lines
ylength(k,1)=perpgrad(k,1)*xlength(k,1)

x1(k,1)=xcoord(k,1)+xlength(k,1)
y1(k,1)=ycoord(k,1)+ylength(k,1)

x2(k,1)=xcoord(k,1)-xlength(k,1)
y2(k,1)=ycoord(k,1)-ylength(k,1)

x=[x1(k,1) x2(k,1)]
y=[y1(k,1) y2(k,1)]

p(k,1)=imline(gca,x,y) %gca is the current axis handle
end

