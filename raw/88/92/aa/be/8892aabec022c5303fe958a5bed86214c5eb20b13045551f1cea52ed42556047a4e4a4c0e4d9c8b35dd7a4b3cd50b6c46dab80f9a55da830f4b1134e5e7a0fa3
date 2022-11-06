%PATCHGENERATOR
%   Creates a patch antenna structure using the method of images 
%   in the ground plane; assumes the uniform mesh.
%   Mutliple (parasitic) patches and multiple probe feeds are allowed. 
%   Uses Matlab mouse input (ginput); see Section 10.4  
%
%   The following parameters need to be specified prior to 
%   calculations:
%   
%   Patch height                                h
%   Plate length (along the x-axis)             L
%   Plate width (along the y-axis)              W
%   Discretization parameter (length)           Nx
%   Discretization parameter (width)            Ny      
%   Number of rectangles of the feeding strip:  Number
%
%   Format:
%   t(4,:)=0 triangles of the metal ground plane
%   t(4,:)=1 triangles of the probe feed
%   t(4,:)=2 triangles of the patch
%   t(4,:)=3 triangles of the upper boundary of dielectric (not seen)
%
%   Note: Matlab function FILL may not work properly in some cases 
%   (usually, when marking the feeding triangles). Pay no attention to it. 
%   The final result (after mouse clicking and pressing Enter) will be 
%   obtained in the correct form
%
%   Copyright 2002 AEMM. Revision 2002/03/22 
%   Chapter 10

clear all
warning off

Color1=[1.00 1.00 1.00];    %dielectric
Color2=[0.45 0.45 0.45];    %ground plane
Color3=[0.70 0.70 0.70];    %patch
Color4=[0.90 0.90 0.90];    %feed

%Separation distance between patch and ground plane
h=0.01;

%Identify finite ground plane
L=0.10;    %Plate length (along the x-axis)
W=0.05;    %Plate width (along the y-axis)
Nx=17;     %Discretization parameter (length)
Ny=9;      %Discretization parameter (width)

%Number of rectangles of the feeding strip:
Number=3;

%Set vertexes of the plate
epsilon=1e-6; 
M=1;
for i=1:Nx+1
    for j=1:Ny+1
        X(M)=-L/2+(i-1)/Nx*L;
        Y(M)=-W/2+(j-1)/Ny*W-epsilon*X(M);
        M=M+1;
    end
end

%Use Delaunay triangulation
TRI = delaunay(X,Y,0); 
t=TRI'; t(4,:)=0;
p=[X; Y; -h*ones(1,length(X))];
save plate p t 


%Identify triangles of the patch
%Use RETURN key to stop the process
PatchNumber=[];
viewer plate; view(0,90); hold on
m=0;
while ~isempty(t)
    m=m+1;
    [xi,yi]=ginput(1);
    TriangleNumber = tsearch(X,Y,TRI,xi,yi);
    n=t(1:3,TriangleNumber);
    PatchNumber= [PatchNumber TriangleNumber];
    x= p(1,n);
    y= p(2,n);
    if isempty(xi|yi) break; end
    fill(x,y,Color4)
    clear xi yi
end
t(4,:)=3;
t(4,PatchNumber)=2;
save plate p t

%Identify junction edge(s)
%Use RETURN key to stop the process
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
    clear xi yi;
end

%Create structure patch-ground plane
tbase=t; pbase=p;
tbase(4,:)=0;
p(3,:)=p(3,:)+h;

T=[tbase t+length(pbase)];
T(4,:)=[tbase(4,:) t(4,:)];
P=[pbase p];
p=P; t=T;

%Create the probe feed (code is too long - should be simplified!)
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

save patch p t h FeedingTriangle
hold off
clear figure
viewer patch
    
    
    




