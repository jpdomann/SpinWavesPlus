function [h_vec,h_light] = spinPlot( x, m, vecScale, atomScale,nFaces )
%SPINPLOT Summary of this function goes here
%   Detailed explanation goes here

numPoints = size(x,1);

%lattice point spacing
spacing = mean(diff(unique(x)));
spacing = min(spacing(spacing>0));

if atomScale > 0
    %atom radius
    rs = atomScale*spacing;
    
    [Xs,Ys,Zs] = sphere(nFaces);
    xs = zeros([numPoints,size(Xs)]);
    ys = xs;
    zs = xs;
    
    h = matlab.graphics.chart.primitive.Surface;
    h = repmat(h,1,numPoints);
    for i = 1:numPoints
        xs(i,:,:) = rs*Xs+x(i,1);
        ys(i,:,:) = rs*Ys+x(i,2);
        zs(i,:,:) = rs*Zs+x(i,3);
        
        xTemp = reshape(xs(i,:,:),size(Xs));
        yTemp = reshape(ys(i,:,:),size(Xs));
        zTemp = reshape(zs(i,:,:),size(Xs));
        
        h(i) = surf(xTemp,yTemp,zTemp);
        hold on
    end
    shading interp
    %material([ka kd ks n sc]) sets the ambient/diffuse/specular strength, specular exponent, and specular color reflectance of the objects.
    material([1 1 1 7 1]) 
    h_light = light('Position',[-1,-1,1]);
    
    h_ax = gca;
    h_ax.AmbientLightColor = [1 1 1];
    
    for i = 1:numPoints
        h(i).FaceLighting = 'phong';
        h(i).FaceColor = [255 50 50]/255;
        h(i).FaceAlpha = 1.0;
    end
end

h_vec = quiver3(x(:,1),x(:,2),x(:,3),m(:,1),m(:,2),m(:,3),vecScale);

grid on
axis equal

end

