clear
clc
clf
addpath('functions')

%% Global Variables
global gamma
gamma = 1e9;

%% Load Run Settings
%Number of lattice sites
nX = 5;
nY = 1;
nZ = 1;
%Lattice Dimensions
Lx = 10e-9;
Ly = 0;
Lz = 0;

%Initial Spin Orientation
S0 = [0 0 1];

%Spin Plot Settings
vecScale = Lx/10;   %size of spin vectors
atomScale = .25;    %size of atoms (set =0 to turn off atom display)
nFaces = 25;        %number of sides on the atom (sphere(nFaces))

%Main Routine Options
H0_mag = 10;          %Uniform Magnetic Bias Field magnitude
H0_dir = [1 0 0];    %Uniform Magnetic Bias Field direction

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

%% Run Program

%Update magnetic field direction for each lattice site
options.H0_dir = repmat(H0_dir,size(X,1),1);
options.H0_mag = repmat(H0_mag,size(X,1),1);
options.ODE_options = ODE_options;
options.matrixSize = size(X);

%% Call LLG function
[t,x,s] = LLG_Mechanics(X,spin,options);

%% Animate
tipStorage = {};
tipStorage{size(x,1)} = [];
col = [];
cols = colormap(parula(numel(t)));
for i = 1:numel(t)
%     %update atom position
%     h_X.XData = h_X.XData + X(:,1) - x(:,1,i);
%     h_X.YData = h_X.YData + X(:,2) - x(:,2,i);
%     h_X.ZData = h_X.ZData + X(:,3) - x(:,3,i);
%     
%     %update spin position
%     h_spin.XData = x(:,1,i);
%     h_spin.YData = x(:,2,i);
%     h_spin.ZData = x(:,3,i);
%     
%     %update spin vector
%     h_spin.UData = s(:,1,i);
%     h_spin.VData = s(:,2,i);
%     h_spin.WData = s(:,3,i);
        
    figure(1)
    clf
    grid on
    [~,~,~,tip] = spinPlot( x(:,:,i), s(:,:,i), vecScale,atomScale,nFaces);
    hold on
    
    %plot tip position
%     tips = [tips;tip];
%     col = [col; repmat(cols(i,:),size(tip,1),1)];
        

    %     plot3(tips(:,1),tips(:,2),tips(:,3),'.')
    for j = 1:size(tip,1)
        tipStorage{j} = [tipStorage{j}; tip(j,:)];
        tips = tipStorage{j};
        col = t(1:i)/max(t);
        surface([tips(:,1) tips(:,1)],...
            [tips(:,2) tips(:,2)],...
            [tips(:,3) tips(:,3)],...
            'Cdata',[col col],...
            'FaceColor','none',...
            'EdgeColor','interp',...
            'LineWidth',2,...
            'Marker','none')
    end
    
    xlim([0,Lx])
    a = 1e-9;
    ylim([-a,a])
    zlim([-a,a])
    view(3)
    pause(.001)
end

disp('fin')
