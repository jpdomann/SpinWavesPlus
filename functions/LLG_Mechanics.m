function [t,x,m] = LLG_Mechanics( X,Spin,options,BCs )
%LLG_MECHANICS is used to analyzed magnon-phonon interaction (coupled
%spin and acoustic waves)
%   This program couples the landau-lifshitz-gilbert equation and the
%   acoustic wave equation. Coupling is accomplished through 
%

%% Create Anonymous ODE Function to be Solved
% Monotonically increasing time 
global tmax
t0 = 0; 
t1 = tmax;
tspan = [t0, t1];
tspan = BCs(1).t;

%Anonymous ODE function
dFdt =@(t,state) LLG_Mechanics_equation(t,state,X,Spin,t0,t1,BCs,options);

%pdepe Initial-boundary value problems for parabolic-elliptic PDEs in 1-D

%% Solve Differential Equation

%Initial spins and lattice positions except for locations with an imposed
%BC
indBC = BCs(end).ind;
Xsolve = X(~indBC,:);
SpinSolve = Spin(~indBC,:);

initial = [Xsolve(:);SpinSolve(:)];

[t,final] = ode23(dFdt,tspan,initial,options.ODE_options);

%% Reconstruct Output State
%input is a column vector, reconfigue into 2D matrices.
x = zeros([size(X),numel(t)]);
m = zeros([size(Spin),numel(t)]);
deleteLast = false;

for i = 1:numel(t)
    xState = final(i,1:numel(final(2,:))/2);
    sState = final(i,1+numel(final(2,:))/2:end);
    
    matrixSize = options.matrixSize;
    nBCs = numel(BCs)-1;
    matrixSize(1) = matrixSize(1) - nBCs;
    
    xtemp = reconstructMatrix(xState,matrixSize);
    mtemp = reconstructMatrix(sState,matrixSize);
  
    %Interpolate BCs to current time point and add to current state
    for j = 1:numel(BCs)-1
        xBCsInterp = interp1(BCs(j).t,BCs(j).x,t(i),'pchip');
        mBCsInterp = interp1(BCs(j).t,BCs(j).m,t(i),'pchip');
        xtemp = reconstructBCs(xtemp,xBCsInterp,BCs(end).ind,j,deleteLast);
        mtemp = reconstructBCs(mtemp,mBCsInterp,BCs(end).ind,j,deleteLast);        
    end
    
    %Incorporate the temporary solutions into the final storage matrix
    x(:,:,i) = xtemp;
    m(:,:,i) = mtemp;    
end


end

