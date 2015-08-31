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

tmax = 50e-9;          %run time

%% Load Run Settings
%Lattice Spacing 
a = 3e-10;

%Number of lattice sites
nX = 25;
nY = 1;
nZ = 1;

%Lattice Dimensions
Lx = a*nX;
Ly = 0;
Lz = 0;

%Initial Magnetization Orientation
M0 = [0 0 1];

%Spin Plot Settings
vecScale = Lx/nX;   %size of spin vectors
atomScale = .25;    %size of atoms (set =0 to turn off atom display)
nFaces = 20;        %number of sides on the atom (sphere(nFaces))

%Main Routine Options

%Bias Field
H0_mag = 0e6;         %Uniform Magnetic Bias Field magnitude
H0_dir = [-1 0 0];    %Uniform Magnetic Bias Field direction

%Numerical options
ODE_options=odeset('RelTol',1*10^(-4),'AbsTol',1*10^(-4),...
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
Spin = repmat(M0,size(X,1),1);

%Plot initial spin setup
figure(1)
clf
grid on
[h_spin,h_X,h_light] = spinPlot( X, Spin, vecScale,atomScale,nFaces);

%% Program Options

%Update magnetic field direction for each lattice site
options.H0_dir = repmat(H0_dir,size(X,1),1);
options.H0_mag = repmat(H0_mag,size(X,1),1);
options.ODE_options = ODE_options;
options.matrixSize = size(X);

%% Boundary Conditions
%Time vector
nt_bc = 1e3;
t_bc = linspace(0, tmax, nt_bc);

%circular frequency
n_cycles = 3;
omega = n_cycles*2*pi/(tmax);

%Boundary Conditions
%Control location and spin of first element at x = [0 0 0]
BCs(1).t = t_bc;
mtemp = [zeros(nt_bc,1), sin(omega.*t_bc)', cos(omega.*t_bc)'];
xtemp = [zeros(nt_bc,1), zeros(nt_bc,1), zeros(nt_bc,1)];

%stop rotating spin after one cycle
condition = t_bc > tmax/n_cycles;
mtemp(condition,:) = repmat(M0,sum(condition),1);

BCs(1).m = mtemp;
BCs(1).x  = xtemp;

%Control location and spin of last element at x = [Lx 0 0]
% BCs(2).t = t_bc;
% BCs(2).m = [zeros(nt_bc,1), sin(omega.*t_bc)', cos(omega.*t_bc)'];
% BCs(2).x = [Lx*ones(nt_bc,1), zeros(nt_bc,1), zeros(nt_bc,1)];

nBCs = numel(BCs);
%The variable ind keeps track of the indices of elements that are subject
%to boundary conditions. Using this in the ode solver will allow explicit
%control of the desired nodes. 
ind = [];
for i = 1:nBCs
    xtemp = BCs(i).x;
    indX = X(:,1) == xtemp(1,1);
    indY = X(:,2) == xtemp(1,2);
    indZ = X(:,3) == xtemp(1,3);
    
    ind = any([ind,all([indX,indY,indZ],2)],2);
    BCs(nBCs+1).ind = ind; 
end

%% Final Options
options.vecScale = vecScale;
options.atomScale = atomScale;
options.nFaces = nFaces;

%% Run Program: Call LLG+Mechanics function
tic
[t,x,m] = LLG_Mechanics(X,Spin,options,BCs);
toc

%% Animate

%Plot initial spin setup
figure(1)
clf
grid on
[h_spin,h_X,h_light] = spinPlot( X, Spin, vecScale,atomScale,nFaces);
xlim([-Lx/nX,Lx+Lx/nX])

ylim([-vecScale,vecScale])
zlim([-vecScale,vecScale])
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
                
%     %plot tip position
%     tip = x(:,:,i) + m1(:,:,i);           
%     for j = 1:size(X,1)
%         tipStorage{j} = [tipStorage{j}; tip(j,:)];
%         tips = tipStorage{j};
%         col = t(1:i)/max(t);
%         %h_tip = surface([tips(:,1) tips(:,1)],[tips(:,2) tips(:,2)],[tips(:,3) tips(:,3)],[col col]);        
%         h_tip(j).XData = [tips(:,1) tips(:,1)];
%         h_tip(j).YData = [tips(:,2) tips(:,2)];
%         h_tip(j).ZData = [tips(:,3) tips(:,3)];
%         h_tip(j).CData = [col col];
%     end

    
    pause(.001)
end

%% Done!
disp('fin')
