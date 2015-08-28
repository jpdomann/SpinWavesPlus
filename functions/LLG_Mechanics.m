function [t,x,s] = LLG_Mechanics( X,Spin,options )
%LLG_MECHANICS is used to analyzed magnon-phonon interaction (coupled
%spin and acoustic waves)
%   This program couples the landau-lifshitz-gilbert equation and the
%   acoustic wave equation. Coupling is accomplished through 
%

%% Create Anonymous ODE Function to be Solved
% Monotonically increasing time 
t0 = 0; 
t1 = 1e-9;
tspan = [t0, t1];

%Anonymous ODE function
dFdt =@(t,state) LLG_Mechanics_equation(t,state,X,Spin,t0,t1,options);

%% Solve Differential Equation

%Initial spins and lattice positions
initial = [X(:);Spin(:)];

[t,final] = ode23(dFdt,tspan,initial,options.ODE_options);

%% Reconstruct Output State
%input is a column vector, reconfigue into 2D matrices.
x = zeros([size(X),numel(t)]);
s = zeros([size(Spin),numel(t)]);

for i = 1:numel(t)
    xState = final(i,1:numel(final(2,:))/2);
    sState = final(i,1+numel(final(2,:))/2:end);
    
    x(:,:,i) = reconstructMatrix(xState,options.matrixSize);
    s(:,:,i) = reconstructMatrix(sState,options.matrixSize);
end


end

