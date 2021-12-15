function [xt,yt,time,angles]=boundaryLTS(i,v1,v2,b1,b2,xt,yt,time,angles)

ttime=time(1,1);
DF=size(xt,1);

% minTheta is the minimum angle to hit the boundary
minTheta=0;
% maxTheta is the maximum angle to hit the boundary
maxTheta=pi;

for j=1:DF
    if angles(j,i)>minTheta && angles(j,i)<maxTheta
        % Calculate the angle of refraction and time remaining after
        % hitting the next boundary for each trajectory.
        angles(j,i+1)=(pi/2)-asin((v2/v1)*sin((pi/2)-angles(j,i))); % base on F types
        time(j,i+1)=time(j,i)-((b2-b1)/sin(angles(j,i)))/v1; %
        if time(j,i+1)>0
            % If the next boundary will be hit then the next x and y
            % coordinates are set to where contact occurs
            xt(j,i+1)=xt(j,i)+((b2-b1)/sin(angles(j,i)))*cos(angles(j,i)); % xt(j,i)+time(j,i+1)*rvel(1)
            yt(j,i+1)=b2;
        elseif time(j,i+1)<0 && time(j,i)>0
            % If the next boundary will not be hit but the trajectory will
            % pass the current boundary the next x and y coordinates are
            % set to where the trajectory will stop
            xt(j,i+1)=xt(j,i)+time(j,i)*v1*cos(angles(j,i));
            yt(j,i+1)=yt(j,i)+time(j,i)*v1*sin(angles(j,i));
        else
            % Otherwise the end of the trajectory has been reached
            xt(j,i+1)=0;
            yt(j,i+1)=0;
        end
    elseif i==1
        % All trajectories should have 1 non-zero coordinate, this makes
        % sure that happens for trajectories which do not hit the boundary
        xt(j,i+1)=ttime*v1*cos(angles(j,i));
        yt(j,i+1)=ttime*v1*sin(angles(j,i));
    else
        xt(j,i+1)=0;
        yt(j,i+1)=0;
    end
    if i==size(time,2)-1 && time(j,i+1)>0
        % This makes sure that trajectories in the final region are
        % calculated correctly
        xt(j,i+2)=xt(j,i+1)+(time(j,i+1)*v2*cos(angles(j,i+1)));
        yt(j,i+2)=yt(j,i+1)+(time(j,i+1)*v2*sin(angles(j,i+1)));
    end
end

        
