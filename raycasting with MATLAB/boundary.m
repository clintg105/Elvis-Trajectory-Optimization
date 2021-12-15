function [hull]=boundary(xt,yt,xp,yp,xm,ym)

reachable=zeros(10000,3); a=1;
for i=1:size(xt,1)
    for j=2:1:size(xt,2)
        if xt(i,j)~=0
            reachable(a,1)=xt(i,j);
            reachable(a,2)=yt(i,j);
            a=a+1;
        end
    end
end
for i=1:size(xp,1)
    for j=2:2:size(xp,2)
        if xp(i,j)~=0
            reachable(a,1)=xp(i,j);
            reachable(a,2)=yp(i,j);
            reachable(a+1,1)=xm(i,j);
            reachable(a+1,2)=ym(i,j);
            a=a+2;
        end
    end
end
array=reachable(1:a-1,:);

X=array(:,1);
Y=array(:,2);
k=convhull(X,Y);

n=length(k);
hull=zeros(n,2);
for i=1:n
    hull(i,:)=[X(k(i)),Y(k(i))];
end