function[xp,xm,yp,ym]=reflections(i,b,v,xp,xm,yp,ym,ttime)

p1=2*i-1; p2=p1+1;
v1=v(i); v2=v(i+1);
DF=size(xp,1)-1; thc=asin(v1/v2);
if i==1
    % xpCrit is the positive x-coordinate of the trajectory which travels
    % at the critical angle to reach the next boundary.  
    xpCrit=((b(i+1)-b(i))/sin((pi/2)-thc))*cos((pi/2)-thc);
    % xpMax is the positive x-coordinate of the trajectory which travels at
    % the critical angle to reach the next boundary and then travels
    % horizontally for time remaining
    xpMax=xpCrit+v2*(ttime-((b(i+1)-b(i))/sin((pi/2)-thc))/v1);
    % xmCrit and xmMax are identical to xpCrit and xpMax except for
    % negative x-coordinates
    xmCrit=-((b(i+1)-b(i))/sin((pi/2)-thc))*cos((pi/2)-thc);
    xmMax=xmCrit-v2*(ttime-((b(i+1)-b(i))/sin((pi/2)-thc))/v1);
    % xp and yp contain the x and y coordinates of reflections with
    % positive x-coordinates.  Coordinates are listed in pairs.  xm and ym
    % contain the same information for reflections with negative
    % x-coordinates.
    xp(DF+1,p1:p2)=[xpCrit,xpMax]; yp(DF+1,p1:p2)=[b(i+1),b(i+1)];
    xm(DF+1,p1:p2)=[xmCrit,xmMax]; ym(DF+1,p1:p2)=[b(i+1),b(i+1)];
    xp(1:DF,p1)=transpose(linspace(xpCrit,xpMax,DF)); yp(1:DF,p1)=b(i+1)*ones(DF,1);
    xm(1:DF,p1)=transpose(linspace(xmCrit,xmMax,DF)); ym(1:DF,p1)=b(i+1)*ones(DF,1);
    for j=1:DF
        dx=v1*cos((pi/2)-thc)*(ttime-((b(i+1)-b(i))/sin((pi/2)-thc))/v1-(xp(j,p1)-xpCrit)/v2);
        dy=v1*sin((pi/2)-thc)*(ttime-((b(i+1)-b(i))/sin((pi/2)-thc))/v1-(xp(j,p1)-xpCrit)/v2);
        xp(j,p2)=xp(j,p1)+dx; yp(j,p2)=yp(j,p1)-dy;
        xm(j,p2)=xm(j,p1)-dx; ym(j,p2)=ym(j,p1)-dy;
    end
else
    ltc=(b(i+1)-b(i))/sin((pi/2)-thc); 
    ttc=ltc/v1;  xc=ltc*cos((pi/2)-thc);
    composite=zeros(i-1,3);
    for k=1:i-1
        composite(k,1)=(pi/2)-asin((v(k)/v(k+1))*sin((pi/2)-thc))
        composite(k,2)=(b(k+1)-b(k))/sin(composite(k,1));
        composite(k,3)=composite(k,2)/v(k);
        ttc=ttc+composite(k,3);
        xc=xc+(composite(k,2)*cos(composite(k,1)));
    end
    xpCrit=+xc; xpMax=xpCrit+(ttime-ttc)*v(i+1);
    xp(DF+1,p1:p2)=[xpCrit,xpMax]; yp(DF+1,p1:p2)=[b(i+1),b(i+1)];
    xp(1:DF,p1)=transpose(linspace(xpCrit,xpMax,DF)); yp(1:DF,p1)=b(i+1)*ones(DF,1);
    xmCrit=-xc; xmMax=xmCrit-(ttime-ttc)*v(i+1);
    xm(DF+1,p1:p2)=[xmCrit,xmMax]; ym(DF+1,p1:p2)=[b(i+1),b(i+1)];
    xm(1:DF,p1)=transpose(linspace(xmCrit,xmMax,DF)); ym(1:DF,p1)=b(i+1)*ones(DF,1);
    for j=1:DF
        dx=v1*cos((pi/2)-thc)*(ttime-ttc-(xp(j,p1)-xpCrit)/v2);
        dy=v1*sin((pi/2)-thc)*(ttime-ttc-(xp(j,p1)-xpCrit)/v2);
        xp(j,p2)=xp(j,p1)+dx; yp(j,p2)=yp(j,p1)-dy;
        xm(j,p2)=xm(j,p1)-dx; ym(j,p2)=ym(j,p1)-dy;
    end    


%     ltCrit=(b(i+1)-b(i))/sin((pi/2)-thc);
%     ttCrit=ltCrit/v(i); xCrit=ltCrit*cos((pi/2)-thc);
%     % Composite contains (1) angle (2) distance (3) time
%     composite=zeros(i-1,3);
%     for k=1:i-1
%         composite(k,1)=(pi/2)-asin((v(k)/v(k+1))*sin((pi/2)-thc));
%         composite(k,2)=(b(k+1)-b(k))/sin(composite(k,1));
%         composite(k,3)=composite(k,2)/v(k);
%         ttCrit=ttCrit+composite(k,3);
%         xCrit=xCrit+(composite(k,2)*cos(composite(k,1)));
%     end
%     distance=(ttime-ttCrit)*v(i+1);
%     xmCrit=-xCrit; xmMax=xmCrit-distance;
%     xm(1:DF,p1)=transpose(linspace(xmCrit,xmMax,DF)); ym(1:DF,p1)=b(i+1)*ones(DF,1);
%     xpCrit=xCrit; xpMax=xpCrit+distance;
%     xp(1:DF,p1)=transpose(linspace(xpCrit,xpMax,DF)); yp(1:DF,p1)=b(i+1)*ones(DF,1);
%     for k=1:DF
%         timeleft=ttime-ttCrit-(abs(xm(k,p1)-xmCrit)/v(i+1));
%         xm(k,p2)=xm(k,p1)-timeleft*v(i)*cos(thc);
%         ym(k,p2)=ym(k,p1)-timeleft*v(i)*sin(thc);
%         xp(k,p2)=xp(k,p1)+timeleft*v(i)*cos(thc);
%         yp(k,p2)=yp(k,p1)-timeleft*v(i)*sin(thc);
%     end
end