function [dFdt, x,spin ] = LLG_Mechanics_equation(t,state,X,Spin,t0,t1,options )
%WAVESWAVESWAVES is used to analyzed magnon-phonon interaction (coupled
%spin and acoustic waves)
%   This program couples the landau-lifshitz-gilbert equation and the
%   acoustic wave equation. Coupling is accomplished through 
%

global gamma mu0 alpha Ms Aex

%% Reconstruct Input state
%input is a column vector, reconfigue into 2D matrices.
xState = state(1:numel(state)/2);
sState = state(1+numel(state)/2:end);

x = reconstructMatrix(xState,options.matrixSize);
m = reconstructMatrix(sState,options.matrixSize);

%% Parse out options
parseStructure(options)

%% Effective magnetic field 
%Effective magnetic field at each lattice site
H_eff_options. H0_dir = options.H0_dir;
H_eff_options. H0_mag = options.H0_mag;

[ H_eff ] = H_effective( X,x,H_eff_options);


%% LLG Equation
dMdt = -mu0*abs(gamma) * cross(m,H_eff) - alpha*mu0*abs(gamma) * cross(m,cross(m,H_eff));
dMdt = dMdt(:);

%% Mechanics
dxdt = zeros(size(dMdt));

%% Output
dFdt = [dxdt;dMdt];

end

