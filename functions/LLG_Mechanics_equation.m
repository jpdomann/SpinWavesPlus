function [dFdt, x,spin,t ] = LLG_Mechanics_equation(t,state,X,Spin,t0,t1,BCs,options )
%WAVESWAVESWAVES is used to analyzed magnon-phonon interaction (coupled
%spin and acoustic waves)
%   This program couples the landau-lifshitz-gilbert equation and the
%   acoustic wave equation. Coupling is accomplished through 
%

global tmax gamma mu0 alpha Ms Aex
persistent t_old

if isempty(t_old);
    t_old = t;
end
%% Reconstruct Input state
%input is a column vector, reconfigue into 2D matrices.
xState = state(1:numel(state)/2);
sState = state(1+numel(state)/2:end);

matrixSize = options.matrixSize;
nBCs = numel(BCs)-1;
matrixSize(1) = matrixSize(1) - nBCs;
x = reconstructMatrix(xState,matrixSize);
m = reconstructMatrix(sState,matrixSize);

%% Interpolate BCs to current time point and add to current state
for i = 1:numel(BCs)-1
    xBCsInterp = interp1(BCs(i).t,BCs(i).x,t,'pchip');
    mBCsInterp = interp1(BCs(i).t,BCs(i).m,t,'pchip');
    x = reconstructBCs(x,xBCsInterp,BCs(end).ind,i);
    m = reconstructBCs(m,mBCsInterp,BCs(end).ind,i);    
end


%% Parse out options
parseStructure(options)

%% Effective magnetic field 
%Effective magnetic field at each lattice site
H_eff_options.H0_dir = options.H0_dir;
H_eff_options.H0_mag = options.H0_mag;
H_eff_options.Options = options;
H_eff_options.t = t;
H_eff_options.tspan = BCs(1).t;

[ H_eff ] = H_effective( X,x,m,H_eff_options);


%% LLG Equation
%Vector form
dMdt = -mu0*abs(gamma)*( cross(m,H_eff) + alpha * cross(m,cross(m,H_eff)) );

%Remove rows where BCs are applied
indBCs = BCs(end).ind;
dMdt(indBCs,:) = [];

dMdt = dMdt(:);

%% Component form
% dm1dt = -mu0*abs(gamma)*(...
%     LeviCivita([1 2 3]).*m(:,2).*H_eff(:,3) + ...
%     LeviCivita([1 3 2]).*m(:,3).*H_eff(:,2) + ...
%     alpha .* m(:,1).*(m(:,1).*H_eff(:,1) + m(:,2).*H_eff(:,2) + m(:,3).*H_eff(:,3) )- ...
%     alpha .* (m(:,1).*m(:,1) + m(:,2).*m(:,2) + m(:,3).*m(:,3)).*H_eff(:,1)  );
% 
% dm2dt = -mu0*abs(gamma)*(...
%     LeviCivita([2 3 1]).*m(:,3).*H_eff(:,1) + ...
%     LeviCivita([2 1 3]).*m(:,1).*H_eff(:,3) + ...
%     alpha .* m(:,2).*( m(:,1).*H_eff(:,1) + m(:,2).*H_eff(:,2) + m(:,3).*H_eff(:,3) )- ...
%     alpha .* (m(:,1).*m(:,1) + m(:,2).*m(:,2) + m(:,3).*m(:,3)).*H_eff(:,2)  );
% 
% dm3dt = -mu0*abs(gamma)*(...
%     LeviCivita([3 1 2]).*m(:,1).*H_eff(:,2) + ...
%     LeviCivita([3 2 1]).*m(:,2).*H_eff(:,1) + ...
%     alpha .* m(:,3).*(m(:,1).*H_eff(:,1) + m(:,2).*H_eff(:,2) + m(:,3).*H_eff(:,3) )- ...
%     alpha .* (m(:,1).*m(:,1) + m(:,2).*m(:,2) + m(:,3).*m(:,3)).*H_eff(:,3)  );
% dmdt = [dm1dt; dm2dt; dm3dt];

%% Mechanics
dxdt = zeros(size(dMdt));

%% Output
dFdt = [dxdt;dMdt];

%update time for next loop
t_old = t;
disp(tmax - t);

end

