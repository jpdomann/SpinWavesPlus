function [h_vec] = spinPlot( x, m, vecScale, atomScale,nFaces )
%SPINPLOT Summary of this function goes here
%   Detailed explanation goes here

numPoints = size(x,1);

%lattice point spacing
spacing = mean(diff(x));
spacing = min(spacing(spacing>0));

if atomScale > 0
    %atom radius
    rs = atomScale*spacing;
    
    [Xs,Ys,Zs] = sphere(nFaces);
    xs = zeros([numPoints,size(Xs)]);
    ys = xs;
    zs = xs;
    
    h = zeros(1,numPoints);
    for i = 1:numPoints
        xs(i,:,:) = rs*Xs+x(i,1);
        ys(i,:,:) = rs*Ys+x(i,2);
        zs(i,:,:) = rs*Zs+x(i,3);
        
        xTemp = reshape(xs(i,:,:),size(Xs));
        yTemp = reshape(ys(i,:,:),size(Xs));
        zTemp = reshape(zs(i,:,:),size(Xs));
        
        h(i) = surfl(xTemp,yTemp,zTemp);
        hold on
    end
    shading interp
    material shiny
    
    for i = 1:numPoints
        h(i).FaceLighting = 'phong';
        h(i).EdgeAlpha = 0.5;
        % h(i).FaceColor = [255 50 50]/255;
        % h(i).FaceAlpha = .6;
    end
end

h_vec = quiver3(x(:,1),x(:,2),x(:,3),m(:,1),m(:,2),m(:,3),vecScale);

grid on
axis equal

end

