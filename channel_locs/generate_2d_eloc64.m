clear all;
close all;
clc;

A=imread('eloc2d.png');
[centers, radii, metric] = imfindcircles(A,[5 60]);
viscircles(centers,radii)
axis equal;
hold on;

N=max(size(centers));

indices=zeros(1,N);
for i=1:N
    a=gtext(num2str(i)); 
    
    mpos=a.Position(1:2);
    [~,ind]=min(sqrt((centers(:,1)-mpos(1)).^2+(centers(:,2)-mpos(2)).^2));
    indices(i)=ind;
    
    plot(a.Position(1), a.Position(2),'r.')

end

eloc642d=centers(indices,:);

figure;
hold on
plot(eloc642d(:,1),eloc642d(:,2),'.');
for i=1:N
    text(eloc642d(i,1),eloc642d(i,2),num2str(i));
end
    