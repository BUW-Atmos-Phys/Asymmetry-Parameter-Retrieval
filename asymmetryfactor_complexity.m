% scain_s --- scattering intensity distribution ----
% size_s ---- particle size ----- 
% instrument ---- PHIPS-HALO or PHIPS-POL ----- 
%  g ---- asymmetry factor  ---- 
%  Cp ---- complexity ----- 
%  C0 ---- integral of geometrical optics phase function ----- 

function [g, Cp, C0, legcoefs_m] = asymmetryfactor_complexity(scain_s, size_s, instrument)

if nargin < 3
    instrument = 'PHIPS-HALO';
    disp('PHIPS-HALO configuration assumed.')
end

Nsignals = length(size_s);
% coefficient for polynomial-fitting of the asymmetry factor; 
poc = [-5.92702829979099e-05,0.00130522506929561,...
    -0.0108666805589366,0.0409324338201189,0.940294146231213];
maxdeg= 8;
maxdeg= 3000;

deg_range=(1:maxdeg)-1;
legnodes= legpts(maxdeg);
nodes_indeg = rad2deg(acos(legnodes));
gm = zeros(Nsignals,1);
Cp = gm;
C0 = gm;

for k = 1:Nsignals
    P11_m = scain_s(k,:);
    if strcmp(instrument,'PHIPS-HALO')
        angs_det = linspace(18,170,20);
    elseif strcmp(instrument,'PHIPS-POL')
        %angs_det = linspace(6,166,21);
        angs_det = linspace(14,166,20);
    else
        disp('Instrument not defined.')
        return
    end
    % Look for NaNs and adjust angles of detection
    ind = isnan(P11_m); P11_m(ind) = []; angs_det(ind) = [];
    P11_m = flip(P11_m'); angs_det = flip(angs_det');
    method = 'nearest';
    vq1 = interp1(angs_det,P11_m,nodes_indeg,method, 'extrap') ;
    method = 'makima';
    vq2 = interp1(angs_det,P11_m,nodes_indeg,method, 'extrap') ;
    method = 'pchip';
    vq3 = interp1(angs_det,P11_m,nodes_indeg,method, 'extrap') ;
    method = 'linear';
    vq4 = interp1(angs_det,P11_m,nodes_indeg,method, 'extrap') ;
    vq =(vq1 + vq2 + vq4)/3;
    legcoefs_m = legvals2legcoeffs(vq);
    C0(k) = legcoefs_m(1); % scaling factor for normalization or total scattering cross section
    legcoefs_nor=legcoefs_m/legcoefs_m(1);
    g_m = legcoefs_nor(2)/3;
    
    % asymmetry factors; 
    coef = legcoefs_nor./(deg_range'*2+1);
    IoC = 1/sum(abs(coef));
    Cp(k) = IoC;
    gm(k) = g_m; 
    
end
gdif = polyval(poc, log(size_s)); 
% total asymmetry factor;
g = (gm+gdif)/2;

end
