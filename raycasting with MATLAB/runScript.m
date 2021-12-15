clear
hold on

% v is an array containing the velocities in the different regions
v=[2;3;8;4];

% b is an array containing the y-coordinates of the boundaries separating
% different regions.  B is the number of different regions
b=[0;4;10;12]; B=length(b);

% DF1 is the discretization factor for the number of trajectories, e.g. if
% DF1=360 there will be ~360 trajectories emanating from the starting point
DF1=360; 

% DF2 is the discretization factor for the number of reflections at
% boundaries where the velocity increases (STL)
DF2=30;

% ttime is the total time for which each trajectory is calculated
ttime=9;

% angleGen is a function which takes DF1 as an input and outputs an array
% of angles which are used as the basis for individual trajectories.  A
% function is used to make sure that the trajectories are symmetric about
% both axes
angles=angleGen(DF1,B); DFT=size(angles,1); % angles measured from x axis

% time is an array containing the time left for each discretized angle upon
% reaching a boundary (e.g. time(j,2) is the time left after the trajectory
% followed by angle j reaches the second boundary)
time=[ttime*ones(DFT,1),zeros(DFT,B-1)];

% xt and yt contain the x-coordinate and y-coordinate of each trajectory 
% at their point(s) of intersection with each boundary (if they exist).  
% The first coordinates are always (0,0).
% is always 0
xt=zeros(DFT,B+1);        yt=zeros(DFT,B+1);

% xp contains the x-coordinates of the reflections at each boundary which
% are positive (if they exist).  The array is organized with columns in
% pairs (e.g. entries in column 1 and column 2 correspond to the same
% reflected ray with the entry in the 1st column containing the
% x-coordinate on the boundary and the 2nd column containing the
% x-coordinate of the reflection)
% yp contains the y-coordinates corresponding to the x-coordinates in xp
xp=zeros(DF2+1,2*(B-1));    yp=zeros(DF2+1,2*(B-1));

% xm contains the x-coordinates of the reflections at each boundary which
% are negative (if they exist).  The array organization is identical to xp.
% ym contains the y-coordinates corresponding to the y-coordinates in yp
xm=zeros(DF2+1,2*(B-1));   ym=zeros(DF2+1,2*(B-1));

% This loop calculates the x and y coordinates of different trajectories
% and reflections (if they occur).  All calculations are done by functions.
for i=1:B-1 % loop over boundaries - replace with while (nonzeroTimes) loop
    if v(i+1)<v(i)
        % boundaryLTS is used to calculate the x-coordinates and
        % y-coordinates of trajectories which move from a higher velocity
        % to a lower velocity by crossing a boundary
        [xt,yt,time,angles]=boundaryLTS(i,v(i),v(i+1),b(i),b(i+1),xt,yt,time,angles);
    elseif v(i+1)>v(i)
        % reflections is used to calculate the coordinates of reflections
        % which occur when moving from a lower velocity to a higher
        % velocity when crossing a boundary
        [xp,xm,yp,ym]=reflections(i,b,v,xp,xm,yp,ym,ttime);
        % boundaryLTS is used to calculate the x-coordinates and
        % y-coordinates of trajectories which move from a lower velocity
        % to a higher velocity by crossing a boundary
        [xt,yt,time,angles]=boundarySTL(i,v(i),v(i+1),b(i),b(i+1),xt,yt,time,angles);
    end
end
[hull]=boundary(xt,yt,xp,yp,xm,yp);

% This loop generates the image output
for j=1:DFT
    if j<size(xp,1)
        for i=1:2:B
            if xp(j,i)~=0 && j<size(xp,1)
                plot([xp(j,i) xp(j,i+1)],[yp(j,i) yp(j,i+1)],'r');
                plot([xm(j,i) xm(j,i+1)],[ym(j,i) ym(j,i+1)],'r');
            end
        end
    elseif j==size(xp,1)
        for i=1:2:B
            if xp(j,i)~=0 && rem(i,2)==0
                plot([xp(j,i) xp(j,i+1)],[yp(j,i) yp(j,i+1)],'b');
                plot([xm(j,i) xm(j,i+1)],[ym(j,i) ym(j,i+1)],'b');
            elseif xp(j,i)~=0
                plot([xp(j,i) xp(j,i+1)],[yp(j,i) yp(j,i+1)],'g');
                plot([xm(j,i) xm(j,i+1)],[ym(j,i) ym(j,i+1)],'g');
            end
        end
    end
    for i=1:B
        if xt(j,i+1)~=0 && rem(i,2)==0
            plot([xt(j,i) xt(j,i+1)],[yt(j,i) yt(j,i+1)],'g');
        elseif xt(j,i+1)~=0 
            plot([xt(j,i) xt(j,i+1)],[yt(j,i) yt(j,i+1)],'b');
        end
    end
end
% for i=1:size(hull,1)-1
%     plot([hull(i,1) hull(i+1,1)],[hull(i,2) hull(i+1,2)],'r');
% end
% plot([hull(size(hull,1),1) hull(1,1)],[hull(size(hull,1),2) hull(1,2)],'r');
    
    

axis equal

    

