function [h_vec,h_atoms,h_light, tip] = spinPlot(x,m,vecScale,atomScale,nFaces )
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
    
    h_atoms = matlab.graphics.chart.primitive.Surface;
    h_atoms = repmat(h_atoms,1,numPoints);
    for i = 1:numPoints
        xs(i,:,:) = rs*Xs+x(i,1);
        ys(i,:,:) = rs*Ys+x(i,2);
        zs(i,:,:) = rs*Zs+x(i,3);
        
        xTemp = reshape(xs(i,:,:),size(Xs));
        yTemp = reshape(ys(i,:,:),size(Xs));
        zTemp = reshape(zs(i,:,:),size(Xs));
        
        h_atoms(i) = surf(xTemp,yTemp,zTemp);
        hold on
    end
    shading interp
    %material([ka kd ks n sc]) sets the ambient/diffuse/specular strength, specular exponent, and specular color reflectance of the objects.
    material([1 1 1 7 1]) 
    h_light = light('Position',[-1,-1,1]);
    
    h_ax = gca;
    h_ax.AmbientLightColor = [1 1 1];
    
    for i = 1:numPoints
        h_atoms(i).FaceLighting = 'phong';
        h_atoms(i).FaceColor = [255 50 50]/255;
        h_atoms(i).FaceAlpha = 1.0;
    end
end
m1 = m.*vecScale;

h_vec = quiver3(x(:,1),x(:,2),x(:,3),m1(:,1),m1(:,2),m1(:,3),0);

% vecMag = (h_vec.UData).^2 + (h_vec.VData).^2 + (h_vec.WData).^2;
tip = x + m1; 

grid on
axis equal

end

