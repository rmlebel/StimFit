function [RF,M] = calcflipRF(RF,z,B1,Gbg)
%   Computes the tip angles from the RF pulse across the slice profile
%   
%   Usage: RF = calcflipRF(RF,z,B1,Gbg)
%   Author: RM Lebel
%   Date: 06/2011
%   
%   Input:
%   RF = RF waveform stucture:
%       .path: path to external waveform file
%       .RF: RF waveform (G)
%       .phase: global phase (degrees)
%       .tau: pulse duration (s)
%       .G: slice-select gradient (G/cm)
%       .ref: refocusing fraction (x2, i.e. near unity for excite; zero for refocus)
%       .angle: prescribed nutation angle (degrees)
%       .alpha: spatial distribution of nutation angles (degrees) (EMPTY)
%   z = vector of slice positions (cm)
%   B1 = relative B1+ scaling factor (unitless, near 1)
%   Gbg = Background field gradient (G/cm) 
%   
%   Output:
%   RF = RF waveform stucture:
%       .path: path to external waveform file
%       .RF: RF waveform (G)
%       .phase: global phase (degrees)
%       .tau: pulse duration (s)
%       .G: slice-select gradient (G/cm)
%       .ref: refocusing fraction (x2, i.e. near unity for excite; zero for refocus)
%       .angle: prescribed nutation angle (degrees)
%       .alpha: spatial distribution of nutation angles (degrees) (POPULATED)
%   M = Magnetization [Mx;My;Mz] after RF pulse

%   Check input arguments
if nargin ~= 4
    error('Function requires four inputs');
end

%   Perform Bloch simulation
M0 = zeros(3,length(RF.angle),length(z));M0(3,:,:) = 1;
[M,RF.alpha]  = pulse_sim(M0,z,RF,Gbg,B1);






%   Force small tip angle regime
%   This is one of the major approximations of this method
%SCL = 1e-4;
%SCL = 1;

% %   Define initial magnetization and initialize some variables
% M0 = zeros(3,length(z));M0(3,:) = 1;
% M  = zeros(3,length(z),length(RF.angle));
% RFtemp = RF;
% 
% for i = 1:length(RF.angle)
%     
%     %   Get current waveform and scale for B1
%     RFtemp.RF = SCL * B1 * RF.RF(i,:);
%     
%     %   Perform Bloch simulation
%     M(:,:,i) = pulse_sim(M0,z,RFtemp,Gbg,B1);
%     
%     %   Compute tip angle
%     %   Is this right?!?!?
%     RF.alpha(i,:) = (1/SCL) * real(acos(M(3,:,i)));
% end

