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
mState = state(1+numel(state)/2:end);

matrixSize = options.matrixSize;
nBCs = numel(BCs)-1;
matrixSize(1) = matrixSize(1) - nBCs;
x = reconstructMatrix(xState,matrixSize);
m = reconstructMatrix(mState,matrixSize);

%% Interpolate BCs to current time point and add to current state
for i = 1:numel(BCs)-1
    xBCsInterp = interp1(BCs(i).t,BCs(i).x,t,'pchip');
    mBCsInterp = interp1(BCs(i).t,BCs(i).m,t,'pchip');
    x = reconstructBCs(x,xBCsInterp,BCs(end).ind,i);
    m = reconstructBCs(m,mBCsInterp,BCs(end).ind,i);    
end

%normalize m so magnitude is 1
mMag = sum(abs(m),2);
m = bsxfun(@rdivide,m,mMag); 

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


%% Smooth the numerical derivative 
% N = 3;                 % Order of polynomial fit
% F = 5;                % Window length
% [b,g] = sgolay(N,F);   % Calculate S-G coefficients
% 
% HalfWin  = ((F+1)/2) -1;
% 
% for n = (F+1)/2:size(x,1)-(F+1)/2,
%   % Zeroth derivative (smoothing only)
%   SG0(n,1) = dot(g(:,1),m(n - HalfWin:n + HalfWin,1));
%   SG0(n,2) = dot(g(:,1),m(n - HalfWin:n + HalfWin,2));
%   SG0(n,3) = dot(g(:,1),m(n - HalfWin:n + HalfWin,3));
% %   
%   % 1st differential
%   SG1(n,1) = dot(g(:,2),m(n - HalfWin:n + HalfWin,1));
%   SG1(n,2) = dot(g(:,2),m(n - HalfWin:n + HalfWin,2));
%   SG1(n,3) = dot(g(:,2),m(n - HalfWin:n + HalfWin,3));
% 
% %   % 2nd differential
% %   SG2(n) = 2*dot(g(:,3)',y(n - HalfWin:n + HalfWin))';
% end


%% Constant Magnitude M Constraint
%Constraint so that the magnitude of M is constant 
%determine whether vector is increasing or decreasing
resid = dot(m,dMdt,2);
if max(abs(resid)) > 1e-14
    dMdt(:,3) = -dot(m(:,1:2),dMdt(:,1:2),2)./m(:,3);
    resid2 = dot(m,dMdt,2);
end

% check = sum(dMdt,2);
% if any(check>0)
%    disp('Non constant change') 
% end

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

