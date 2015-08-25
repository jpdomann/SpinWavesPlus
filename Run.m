clear
clc
clf
addpath('functions')

%% Load Run Settings

%% Create lattice
n_points = 20;
X = linspace(0,1,n_points);
Y = zeros(size(X));
Z = Y;
Y = X;

%Grid of lattice points
[X,Y,Z] = meshgrid(X,Y,Z);
%Vectorize
X = unique([X(:) Y(:) Z(:)],'rows');
%Size of lattice
nX = numel(unique((X(:,1))));
nY = numel(unique((X(:,2))));
nZ = numel(unique((X(:,3))));

plot3(X(:,1),X(:,2),X(:,3),'.')
xlabel('X')
ylabel('Y')
zlabel('Z')
grid on

%% Initialize spins
spin0 = [0 0 1]; 
spin0 = repmat(spin0,size(X,1),1);

%Plot initial spin setup
figure(1)
clf
grid on
vecScale = 1;
atomScale = 0;
nFaces = 6;
tic
spinPlot( X, spin0, vecScale,atomScale,nFaces)
toc
light('Position',[0,0,1])


%% Run Program

% Call main function
SpinWavesPlus(options)

