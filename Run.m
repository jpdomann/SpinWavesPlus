clear
clc
clf
addpath('functions')

%% Global Variables
global gamma mu0 alpha Ms Aex tmax
gamma = 1e9;            % - gyromagnetic ratio
mu0 = 4*pi*10^(-7);     % - permeability of free space
alpha = 0.05;            %unitless - Gilbert Damping
Ms = 4.8e5;             %A/m - saturation magnetization
Aex = 1.05e-11;         %J/m - exchange strength

tmax = 100e-9;

%% Load Run Settings
%Number of lattice sites
nX = 1;
nY = 1;
nZ = 1;
%Lattice Dimensions
Lx = 10e-9;
Ly = 0;
Lz = 0;

%Initial Spin Orientation
S0 = [0 0 1];

%Spin Plot Settings
vecScale = Lx/nX;   %size of spin vectors
atomScale = .25;    %size of atoms (set =0 to turn off atom display)
nFaces = 20;        %number of sides on the atom (sphere(nFaces))

%Main Routine Options
H0_mag = 5e6;          %Uniform Magnetic Bias Field magnitude
H0_dir = [-1 0 0];    %Uniform Magnetic Bias Field direction

%Numerical options
ODE_options=odeset('RelTol',1*10^(-6),'AbsTol',1*10^(-6),...
    'InitialStep',[],...
    'InitialSlope',0,...
    'Stats','off');
    
%% Create lattice

X = linspace(0,Lx,nX);
Y = linspace(0,Ly,nY);
Z = linspace(0,Lz,nZ);

%Grid of lattice points
[X,Y,Z] = meshgrid(X,Y,Z);
%Vectorize
X = unique([X(:) Y(:) Z(:)],'rows');

clear Y Z

% plot3(X(:,1),X(:,2),X(:,3),'.')
% xlabel('X')
% ylabel('Y')
% zlabel('Z')
% grid on

%% Initialize spins / plot

%initial spin arrangement (uniform)
spin = repmat(S0,size(X,1),1);

%Plot initial spin setup
figure(1)
clf
grid on
[h_spin,h_X,h_light] = spinPlot( X, spin, vecScale,atomScale,nFaces);

%% Program Options

%Update magnetic field direction for each lattice site
options.H0_dir = repmat(H0_dir,size(X,1),1);
options.H0_mag = repmat(H0_mag,size(X,1),1);
options.ODE_options = ODE_options;
options.matrixSize = size(X);

%% Run Program: Call LLG+Mechanics function
[t,x,m] = LLG_Mechanics(X,spin,options);

%% Animate

%Plot initial spin setup
figure(1)
clf
grid on
[h_spin,h_X,h_light] = spinPlot( X, spin, vecScale,atomScale,nFaces);
xlim([-Lx/nX,Lx+Lx/nX])
a = vecScale;
ylim([-a,a])
zlim([-a,a])
view(3)

%initialize variables
tipStorage = {};
tipStorage{size(x,1)} = [];
col = [];
cols = colormap(parula(numel(t)));

for i = 1:size(X,1)
    h_tip(i) = matlab.graphics.chart.primitive.Surface;
    h_tip(i).FaceColor = 'none';
    h_tip(i).EdgeColor = 'interp';
    h_tip(i).LineWidth = 2;
    h_tip(i).Marker = 'none';
    h_tip(i).Parent = gca;
end

%turn off warning that appears during first iteration
s = warning('off','MATLAB:gui:array:InvalidArrayShape');

%Animate!
for i = 1:numel(t)
    
    %update atom position
    for j = 1:size(X,1)
        h_X(j).XData = h_X(j).XData + X(j,1) - x(j,1,i);
        h_X(j).YData = h_X(j).YData + X(j,2) - x(j,2,i);
        h_X(j).ZData = h_X(j).ZData + X(j,3) - x(j,3,i);
    end
    
    %update spin position
    h_spin.XData = x(:,1,i);
    h_spin.YData = x(:,2,i);
    h_spin.ZData = x(:,3,i);
    
    %update magnetization vector
    m1 = m.*vecScale;
    h_spin.UData = m1(:,1,i);
    h_spin.VData = m1(:,2,i);
    h_spin.WData = m1(:,3,i);
                
    %plot tip position
    tip = x(:,:,i) + m1(:,:,i);           
    for j = 1:size(X,1)
        tipStorage{j} = [tipStorage{j}; tip(j,:)];
        tips = tipStorage{j};
        col = t(1:i)/max(t);
        %h_tip = surface([tips(:,1) tips(:,1)],[tips(:,2) tips(:,2)],[tips(:,3) tips(:,3)],[col col]);        
        h_tip(j).XData = [tips(:,1) tips(:,1)];
        h_tip(j).YData = [tips(:,2) tips(:,2)];
        h_tip(j).ZData = [tips(:,3) tips(:,3)];
        h_tip(j).CData = [col col];
    end

    
    pause(.001)
end

%% Done!
disp('fin')
