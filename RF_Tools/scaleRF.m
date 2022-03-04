function RF = scaleRF(RF,tau,angle)
%   Scales an RF waveform to produce a desired tip angle
%   
%   Author: RML
%   Date: 06/2011
%   
%   Usage: RF = scaleRF(RF,tau,angle)
%            ----OR----
%          RF = scaleRF(RF)
%   
%   Input (case 1):
%   RF: RF waveform (arbitrary scaling)
%   tau: pulse duration (s)
%   angle: desired flip angle (degrees)
%   
%   Input (case 2):
%   RF = RF waveform stucture:
%       .path: path to external waveform file
%       .RF: RF waveform (arbitrary scaling)
%       .phase: global phase (degrees)
%       .tau: pulse duration (s)
%       .G: slice-select gradient (G/cm)
%       .ref: refocusing fraction (x2, i.e. near unity for excite; zero for refocus)
%       .angle: prescribed nutation angle (degrees) 
%       .alpha: spatial distribution of nutation angles (degrees)
%   
%   
%   Output (case 1):
%   RF: complex RF waveform (G)
%   
%   Output (case 2):
%   RF = RF waveform stucture:
%       .path: path to external waveform file
%       .RF: RF waveform (G)
%       .phase: global phase (degrees)
%       .tau: pulse duration (s)
%       .G: slice-select gradient (G/cm)
%       .ref: refocusing fraction (x2, i.e. near unity for excite; zero for refocus)
%       .angle: prescribed nutation angle (degrees) 
%       .alpha: spatial distribution of nutation angles (degrees)

%   Define some parameters
gamma = 2*pi*42.575e6;  %   Hz/T

if ~isstruct(RF)
    
    warning('scaleRF: untested');
    
    for i = 1:length(angle)
        %   Get current angle
        alphaC = gamma * tau/length(RF(i,:)) * abs(sum(RF(i,:)));
        alphaC = alphaC * 180/pi;   %   Degrees
        
        %   Compute scaling factor
        SF = angle(i)./alphaC;
        
        %   Scale the RF waveform
        RF(i,:) = SF*RF(i,:)*10000;  %   Gauss
    end
    
elseif isstruct(RF)
    
    for i = 1:length(RF.angle)
        %   Get current angle
        alphaC = gamma * RF.tau/length(RF.RF(i,:)) *abs(sum(RF.RF(i,:)));
        alphaC = alphaC * 180/pi;   %   Degrees
        
        %   Compute scaling factor
        SF = RF.angle(i) ./ alphaC;
        
        %   Scale the RF waveform
        RF.RF(i,:) = SF * RF.RF(i,:) * 10000;     %   Gauss
    end
end
