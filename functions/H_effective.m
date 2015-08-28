function [ H_eff ] = H_effective( X,x,varargin)
%H_EFFECTIVE computes the effective field at all grid / lattice sites 
%   The program generates an energy surface using the desired terms, and
%   then takes the derivative of this surface to compute the effective
%   field
%   inputs: X - refference lattice positions
%           x - current lattice positions

global gamma mu0 alpha Ms Aex

options = varargin{1};

%% Bias Field
H_b_mag = options.H0_mag;
H_b_dir = options.H0_dir;

H_b = bsxfun(@times,H_b_mag,H_b_dir);

%% Exchange Field
H_ex =zeros(size(H_b));


%% Magnetoelastic Field
H_mech = zeros(size(H_b));


%% Total Effective Field
H_eff = H_b + H_ex + H_mech;

end

