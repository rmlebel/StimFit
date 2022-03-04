function [T2, amp, opt] = ExpFit(S,opt)
%   T2 fitting routine with stimulated echo compensation
%   
%   Author: R. Marc Lebel
%   Date:   06/2011
%   
%   Usage:
%   [T2 B1 amp opt] = StimFit(S,opt)
%   
%   Input:
%   S:   Signal vector to be fit (1 x ETL)
%   opt: Options structure defined by StimFit_optset (optional)
%   
%   Output:
%   T2:  Decay time (s)
%   amp: Relative amplitude (au)
%   opt: Modified options structure


%   Check inputs
if nargin < 1
    error('Function requires at least one input');
end
if nargin < 2 || isempty(opt) || ~isstruct(opt)
    opt = StimFit_optset;
end


%   Check to ensure echo train lengths match
if length(S) ~= opt.etl
    error('Inconsistent echo train length');
end


%   Perform bounded non-linear least squares fitting
[X,~,~,~,info] = lsqnonlin(@diff_sig,opt.lsq.Icomp.X0,...
    opt.lsq.Icomp.XL,opt.lsq.Icomp.XU,opt.lsq.fopt);
T2 = X(1); amp = X(2);
if opt.debug
    fprintf('Fitting info:\n\tIterations: %g\n\tFunction Count:%g\n',...
        info.iterations,info.funcCount);
    fprintf('Fitting info: %s\n',info.message);
end


%   Objecitve function for least squares fitting
function delta = diff_sig(X)
    
    %   Compute candidate signal via EPG algorithm
    te = opt.esp:opt.esp:opt.esp*opt.etl;
    Sfit = X(2).*exp(-te./X(1));
    
    
    %   Compute objective function (signal difference)
    delta = S(:) - Sfit(:);
    
    
    %   Plot signal
    if opt.debug
        figure(42575);
        
        subplot(1,2,1);
        plot(te,S,'o',te,Sfit);
        title(sprintf('T2 = %08.3f ms  A = %04.2f',X(1)*1e3,X(2)));

        set(gca,'XLim',[0 te(end)+opt.esp/2],'YLim',[0 Inf],'LineWidth',1);grid on;
        xlabel('Echo time (s)');ylabel('Mt (au)');
        
        
        subplot(1,2,2);
        resid = 100/S(1)*delta;
        plot(te,resid);
        set(gca,'XLim',[0 te(end)+opt.esp/2],'LineWidth',1);grid on;
        xlabel('Echo time (s)');ylabel('Residual (%)');
        title(sprintf('Residuals (MSE = %07.4f)',sum(delta.^2)));
        drawnow;
    end
end


end
