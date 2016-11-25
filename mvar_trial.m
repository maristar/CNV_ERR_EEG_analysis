function [ARF,RCF,PE,DC,varargout] = mvar_trial(Y, Pmax, Mode);
% Inititialization
[N,M] = size(Y);
if iscell(Y)
        Pmax = min(max(N ,M ),Pmax);
        C    = Y;
end;
if nargin<3,
        % according to [1], and other internal validations this is in many cases the best algorithm 
        Mode=2;
end;
[C(:,1:M),n] = covm(Y,'M');
PE(:,1:M)  = C(:,1:M)./n;
if 0,
elseif Mode==2 %%%%% Partial Correlation Estimation: Vieira-Morf Method [2] with unbiased covariance estimation
	%===== In [1] this algorithm is denoted "Nutall-Strand with unbiased covariance" =====%
        %C(:,1:M) = C(:,1:M)/N;
        F = Y;
        B = Y;
        %PEF = C(:,1:M);
        %PEB = C(:,1:M);
        PEF = PE(:,1:M);
        PEB = PE(:,1:M);
        for K = 1:Pmax,
                [D,n]	= covm(F(K+1:N,:),B(1:N-K,:),'M');
                D = D./n;

		       ARF(:,K*M+(1-M:0)) = D / PEB;	
                ARB(:,K*M+(1-M:0)) = D'/ PEF;	
                
                tmp        = F(K+1:N,:) - B(1:N-K,:)*ARF(:,K*M+(1-M:0)).';
                B(1:N-K,:) = B(1:N-K,:) - F(K+1:N,:)*ARB(:,K*M+(1-M:0)).';
                F(K+1:N,:) = tmp;
                
                for L = 1:K-1,
                        tmp      = ARF(:,L*M+(1-M:0))   - ARF(:,K*M+(1-M:0))*ARB(:,(K-L)*M+(1-M:0));
                        ARB(:,(K-L)*M+(1-M:0)) = ARB(:,(K-L)*M+(1-M:0)) - ARB(:,K*M+(1-M:0))*ARF(:,L*M+(1-M:0));
                        ARF(:,L*M+(1-M:0))   = tmp;
                end;
                
                RCF(:,K*M+(1-M:0)) = ARF(:,K*M+(1-M:0));
                RCB(:,K*M+(1-M:0)) = ARB(:,K*M+(1-M:0));
                
                [PEF,n] = covm(F(K+1:N,:),F(K+1:N,:),'M');
                PEF = PEF./n;

		[PEB,n] = covm(B(1:N-K,:),B(1:N-K,:),'M');
                PEB = PEB./n;

                PE(:,K*M+(1:M)) = PEF;        
        end;
        



end;


if any(ARF(:)==inf),
% Test for matrix division bug. 
% This bug was observed in LNX86-ML5.3, 6.1 and 6.5, but not in SGI-ML6.5, LNX86-ML6.5, Octave 2.1.35-40; Other platforms unknown.
p = 3;
FLAG_MATRIX_DIVISION_ERROR = ~all(all(isnan(repmat(0,p)/repmat(0,p)))) | ~all(all(isnan(repmat(nan,p)/repmat(nan,p))));

if FLAG_MATRIX_DIVISION_ERROR, 
	%fprintf(2,'### Warning MVAR: Bug in Matrix-Division 0/0 and NaN/NaN yields INF instead of NAN.  Workaround is applied.\n');
	warning('MVAR: bug in Matrix-Division 0/0 and NaN/NaN yields INF instead of NAN.  Workaround is applied.');

	%%%%% Workaround 
	ARF(ARF==inf)=NaN;
	RCF(RCF==inf)=NaN;
end;
end;	

%MAR   = zeros(M,M*Pmax);
DC     = zeros(M);
for K  = 1:Pmax,
%       VAR{K+1} = -ARF(:,K*M+(1-M:0))';
        DC  = DC + ARF(:,K*M+(1-M:0))'.^2; %DC meausure [3]
end;
% figure;imagesc(DC)
% input('ss')
% Developers:  What's New (Koders Blog) | Getting Started | Add Projects | Downloads
% Add Search to Your Site | Open Source Zeitgeist | Advertise | About Koders
% 
% RSS Feed
% 
% Copyright © 2006 Koders  -  Searching 512,549,710 lines of code.