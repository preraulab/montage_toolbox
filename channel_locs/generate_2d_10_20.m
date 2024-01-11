clear all;
close all;
clc;

A=imread('eloc_10_20.png');
[centers, radii, metric] = imfindcircles(A,[10 30]);
inds=(radii>=18);
centers=centers(inds,:);
radii=radii(inds);
viscircles(centers,radii)
axis equal;
hold on;

N=21;

chan_labels={'Fp1','Fp2','F7','F3','Fz','F4','F8',...
    'A1','T3','C3','Cz','C4','T4','A2',...
    'T5','P3','Pz','P4','T6',...
    'O1','O2'};
    

indices=zeros(1,N);
for i=1:N
    a=gtext(num2str(i)); 
    
    mpos=a.Position(1:2);
    [~,ind]=min(sqrt((centers(:,1)-mpos(1)).^2+(centers(:,2)-mpos(2)).^2));
    indices(i)=ind;
    
    plot(a.Position(1), a.Position(2),'r.')

end

eloc10_20=centers(indices,:);

figure;
hold on
plot(eloc10_20(:,1),eloc10_20(:,2),'.');
for i=1:N
    text(eloc10_20(i,1),eloc10_20(i,2),chan_labels{i});
end

save eloc_10_20.mat centers chan_labels
    