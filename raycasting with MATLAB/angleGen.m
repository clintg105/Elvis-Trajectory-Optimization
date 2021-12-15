function[angles]=angleGen(DF,B)

% angleGen is a function which takes DF1 as an input and outputs an array
% of angles which are used as the basis for individual trajectories.  To
% make sure that the trajectories are symmetric about both axes the angles
% are generated for a single quadrant and then expanded to all quadrants
DF=floor(DF/4);
angles1=[transpose(linspace(0,pi/2,DF)),zeros(DF,B-1)];
angles2=zeros(DF,B); angles3=zeros(DF,B); angles4=zeros(DF,B);
for i=1:DF
    angles2(i,1)=angles1(i,1)+pi/2;
    angles3(i,1)=angles1(i,1)+pi;
    angles4(i,1)=angles1(i,1)+3*pi/2;
end
angles=[angles1;angles2;angles3;angles4];