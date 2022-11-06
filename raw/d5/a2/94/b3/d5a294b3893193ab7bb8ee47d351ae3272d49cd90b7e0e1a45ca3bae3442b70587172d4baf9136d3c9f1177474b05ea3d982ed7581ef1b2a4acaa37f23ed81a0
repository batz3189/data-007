%PATCHGENERATOR1
%   Imports a patch antenna structure (ground plane plus patch)
%   Creates antenna feed(s) using Matlab mouse input
%
%   The following parameters need to be specified prior to 
%   calculations:
%
%   Patch height                                h
%   Number of rectangles of the feeding strip:  Number
%   Binary filename to import                   'filename'
%
%   Format:
%   t(4,:)=0 triangles of the metal ground plane
%   t(4,:)=1 triangles of the probe feed
%   t(4,:)=2 triangles of the patch
%   t(4,:)=3 triangles of the upper boundary of dielectric (not seen)
%
%   Copyright 2002 AEMM. Revision 2002/03/22 
%   Chapter 10

clear all
warning off

Color1=[1.00 1.00 1.00];    %dielectric
Color2=[0.45 0.45 0.45];    %ground plane
Color3=[0.70 0.70 0.70];    %patch
Color4=[0.90 0.90 0.90];    %feed

h=0.15;                     %patch height
Number=2;                       

load('eshape')              %import file           
p(3,:)=-h;

Plate=find(t(4,:)==1);                  %from PDE toolbox
Patch=find((t(4,:)==2)|(t(4,:)==3));    %from PDE toolbox
%   Indexes 1, 2 and 3 can be interchanged if another PDE file is used

t(4,:)=3;
t(4,Patch)=2;
save plate p t;

%Identify junction edge(s)
%Use RETURN key to stop the process
clear figure
hold off
viewer plate; view(0,90); hold on
FeedingTriangle=[];
TRI=t(1:3,:)';
while ~isempty(t)
    [xi,yi]=ginput(1);
    TriangleNumber = tsearch(p(1,:),p(2,:),TRI,xi,yi);
    n=t(1:3,TriangleNumber);
    FeedingTriangle= [FeedingTriangle TriangleNumber];
    x= p(1,n);
    y= p(2,n);
    if isempty(xi|yi) break; end
    fill(x,y,Color4)
    clear xi yi
end
clear figure

%Create structure patch-ground plane
tbase=t; pbase=p;
tbase(4,:)=0;
p(3,:)=p(3,:)+h;

T=[tbase t+length(pbase)];
T(4,:)=[tbase(4,:) t(4,:)];
P=[pbase p];
p=P; t=T;

%Create the probe feed
FeedingTriangle=[FeedingTriangle FeedingTriangle+length(tbase)];
for n=1:length(FeedingTriangle)/4
    FT=[FeedingTriangle(2*n-1) FeedingTriangle(2*n)];
    N=t(1:3,FT(1));
    M=t(1:3,FT(2));
    a=1-all([N-M(1) N-M(2) N-M(3)]);
    Edge_B=M(find(a)); 
    Edge_T  =[Edge_B'+length(pbase)];

    Edge_MM=Edge_B;
    for k=1:Number-1
        p(:,length(p)+1)=k/Number*(p(:,Edge_T(1))-p(:,Edge_B(1)))+p(:,Edge_B(1));
        p(:,length(p)+1)=k/Number*(p(:,Edge_T(2))-p(:,Edge_B(2)))+p(:,Edge_B(2));
        Edge_M=[length(p)-1,length(p)];
        tFeed1(:,k)  =[Edge_MM(1);Edge_MM(2);Edge_M(2);1];
        tFeed2(:,k)  =[Edge_MM(1);Edge_M(1);Edge_M(2);1];
        Edge_MM=Edge_M;
    end
        
    tFeed3  =[Edge_M(1);Edge_M(2);Edge_T(2);1];
    tFeed4  =[Edge_M(1);Edge_T(1);Edge_T(2);1];
    t=[t tFeed1 tFeed2 tFeed3 tFeed4];
end

p=p/10; h=h/10;        %scale properly!!!

save patch p t h FeedingTriangle
hold off
clear figure
viewer patch
    
    




