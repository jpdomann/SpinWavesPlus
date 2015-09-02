function [ H_eff ] = H_effective( X,x,m,varargin)
%H_EFFECTIVE computes the effective field at all grid / lattice sites 
%   The program generates an energy surface using the desired terms, and
%   then takes the derivative of this surface to compute the effective
%   field
%   inputs: X - refference lattice positions
%           x - current lattice positions

global gamma mu0 alpha Ms Aex

options = varargin{1};

%% Diagnostic plot
%Show Starting Conditions

t = options.t;
tspan = options.tspan;
if min(abs(t-tspan)) < 1e-14;
    figure(2)
    subplot(1,2,1)
%     vecScale = options.Options.vecScale;
%     atomScale = options.Options.atomScale;
%     nFaces = options.Options.nFaces;
%     spinPlot( x, m, vecScale,atomScale,nFaces);
    hold off
    plot(m(:,1),'r')
    hold on
    plot(m(:,2),'g')
    plot(m(:,3)-1,'b')
    plot(sum(sqrt(m.^2),2)-1,'c')
    legend('x','y','z-1', 'total-1','Location','Northeast')
end

%% Bias Field
H_b_mag = options.H0_mag;
H_b_dir = options.H0_dir;

H_b = bsxfun(@times,H_b_mag,H_b_dir);

%% Exchange Field
%Second derivatives (see doc del2 for why the 4 is there)
d2xx = 4*del2(m(:,1),mean(diff(x(:,1))));
d2xy = 4*del2(m(:,1),mean(diff(x(:,2))));
d2xz = 4*del2(m(:,1),mean(diff(x(:,3))));

d2yx = 4*del2(m(:,2),mean(diff(x(:,1))));
d2yy = 4*del2(m(:,2),mean(diff(x(:,2))));
d2yz = 4*del2(m(:,2),mean(diff(x(:,3))));

d2zx = 4*del2(m(:,3),mean(diff(x(:,1))));
d2zy = 4*del2(m(:,3),mean(diff(x(:,2))));
d2zz = 4*del2(m(:,3),mean(diff(x(:,3))));

%condition if lattice has restricted dimensions
if all(diff(x(:,1)) == 0)
    d2xx(:) = 0;
    d2yx(:) = 0;
    d2zx(:) = 0;
end
if all(diff(x(:,2)) == 0)
   d2xy(:) = 0;
   d2yy(:) = 0;
   d2zy(:) = 0;
end
if all(diff(x(:,2)) == 0)
    d2xz(:) = 0;
    d2yz(:) = 0;
    d2zz(:) = 0;
end

d2X = nansum([d2xx d2xy d2xz],2);
d2Y = nansum([d2yx d2yy d2yz],2);
d2Z = nansum([d2zx d2zy d2zz],2);


if min(abs(t-tspan)) < 1e-14;
    figure(2)
    subplot(1,2,2)
    h_old = findobj('Parent',gca,'Type','line');
    delete(h_old)
     h(1) = plot(d2xx,'r'); 
    hold on
    h(2) = plot(d2yx,'g'); 
    h(3) = plot(d2zx,'b'); 
    legend('xx','yx','zx')
    pause(.01)
end

%Compute Exchange Field
H_ex = 2*Aex./(mu0.*Ms)*[d2X(:),d2Y(:),d2Z(:)];

%If any columns return NaN, set to zero instead. This occurs when there is
%not a fully 3D array of atoms (ie there isn't a gradient in one or more of
%the directions)
H_ex(isnan(H_ex))=0;

% H_ex =zeros(size(H_b));


%% Magnetoelastic Field
H_mech = zeros(size(H_b));


%% Total Effective Field
H_eff = H_b + H_ex + H_mech;

end

