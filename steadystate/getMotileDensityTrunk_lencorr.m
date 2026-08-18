function [rhoW0] = getMotileDensityTrunk_lencorr(param,radtrunk)
% get linear density of motile mitos in the trunk
% ** Includes correction for finite mito length, accounting for mitos that
% fuse before they enter the region where fission is possible, thereby
% reducing the incoming mobile flux at the trunk ***
% Param includes all the kinetic parameters
% radtrunk = radius of trunk

%%
v = param.v; 
kp = param.kp;
kr = param.kr;
ka = param.kfiss; % reactivation rate
ks = param.As/radtrunk^param.gammas;

% fusion probability in trunk
Pu = param.Au/radtrunk^param.gamma;

% delta and alpha ratios in trunk
d0 = v*Pu/ka;
a0 = v*Pu/2/param.kb;

% mitochondrial unit length
len= param.mitolen;

% polynomial coefficients for 5th deg polynomial, highest to lowest
coeff = zeros(1,6);
coeff(1) = Pu^2*a0^2*len*v^2; % x^5
coeff(2) = Pu*v*(a0^2*(2*v+len*kr) - Pu*len*v*(2*a0-d0)); %x^4
coeff(3) = v*(2*a0^2*kr - 4*Pu*a0*v - 4*Pu*a0^2*kp + 2*Pu^2*len*v - 2*Pu*a0*len*kr + Pu*d0*len*(kr+2*ks)); % x^3
coeff(4) = 2*Pu*v^2 - 4*a0*kr*v - 4*a0^2*kp*kr + 8*Pu*a0*kp*v + 2*Pu*len*v*(kr+ks); % x^2
coeff(5) = 2*kr*v + 8*a0*kp*kr - 4*Pu*kp*v; % x^1
coeff(6) = -4*kp*kr;

rhoW0 = roots(coeff);
% find positive real root that is smaller than the uncorrected density
ind = find(real(rhoW0)>0 & abs(imag(rhoW0))<2*eps & rhoW0< 2*kp/v);
if (length(ind)~=1)
    error('No unambiguous solution for mobile density in trunk')
end
rhoW0 = rhoW0(ind);

%% check that the equation is properly solved
% rhosd = ks*rhoW0/(kr + v*Pu*rhoW0);
% m1 = rhoW0/2 + (rhoW0/2 + rhosd)*(1 + d0*rhoW0)/(a0*rhoW0-1)^2;
% check = v*rhoW0/2 + v/2*Pu*rhoW0*len*m1-kp;

end
