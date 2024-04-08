function [uHatSmooth,sigmaSmooth, stdSmooth, idxWindStarts, lambdaSmooth] = bayesian_smoothing(hr, wlen, kstd, varargin)

    defaultShowPlot = false;
    defaultNoiseCorrPar = 1;
    params = inputParser;
    params.CaseSensitive = false;
    addOptional(params,'showplots',defaultShowPlot,@(x) x == true || x == false);
    addOptional(params,'showreconciliationplots',defaultShowPlot,@(x) x == true || x == false);
    addOptional(params, 'noisecorrpar', defaultNoiseCorrPar)
    %addOptional(params,'maxRebound',defaultRebound,@(x) x > 0 && x < 1000);
    parse(params,varargin{:});

    showplots = params.Results.showplots;
    showreconciliationplots = params.Results.showreconciliationplots;
    noiseCorrPar = params.Results.noisecorrpar;

    hr = hr(:).';

    nSample = length(hr);
    kerWeigthPre = normpdf(0: kstd-1, 0, kstd/2); %weigths of the gaussian kernel for the reconciliation phase
    kerWeigthPost = fliplr(kerWeigthPre); %weigths of the gaussian kernel for the reconciliation phase
    
    %kerWeigth = kerWeigth/norm(kerWeigth,1); %normalized weights
    
    
    % Regularization initialization (regularization factor == gamma)
    par.gammamin = 1e-5; %Minimum value of regularization factor
    par.gammamax = 1e7; %Maximum value of regularization factor
    par.max_iter = 100; %Maximum number of iterations to obtain the regularization parameter
    par.tol = 10e-10; %Tolerance
    
    
    %%%% filterPartData
    %noiseCorrPar = [-1.3, 0.42]; %parameters of the AR(2) describing the measurement noise
    sigma2 = [];
    lambda2 = [];
    uHat = nan(ceil(nSample/wlen),wlen+2*kstd);
    covUHat = nan(ceil(nSample/wlen),wlen+2*kstd);
    hrin = hr';
    if(ceil(nSample/wlen) ~= nSample/wlen)
        missing = (ceil(nSample/wlen) - nSample/wlen) * wlen;
        hrin = [hrin; flipud(hrin(end-missing+1:end))];
    end
    
    idxWindStarts = 1 : wlen : nSample;
    
    for w = 1:ceil(nSample/wlen)
        idxWindStart = idxWindStarts(w);
        idxWindEnd = idxWindStart + wlen - 1;
    
        %add overlap to window
        if idxWindStart == 1 % first window
            idxWindEnd = idxWindEnd + kstd;
        elseif idxWindEnd == length(hrin) % last window
            idxWindStart = idxWindStart - kstd;
        else
            idxWindEnd = idxWindEnd + kstd;
            idxWindStart = idxWindStart - kstd;
        end
    
        window = hrin(idxWindStart:idxWindEnd);
        mirrorHead = flipud(window(1:round(wlen/2))); %innesco
        yWind = [mirrorHead; window];
        n_gap = length(find(isnan(yWind)));
        
        if n_gap >= wlen-1
            uHat(idxWindStart,:) = nan;
            sigma2Hat = nan;
            lambda2Hat = nan;
       
        else
                
            
            [idxNan,~] = find(isnan(yWind));     %find missing samples
            N = length(yWind);
            yWind(idxNan) = [];                  %timeseries without nan 
            n = length(yWind);
            
            vt = [1:1:N];
            vt(idxNan)=[];
            
            Gv = eye(N);
            G = Gv(vt',:);
    
            %SMOOTHING
%             a1 = noiseCorrPar(1); %noise correlation parameters
%             a2 = noiseCorrPar(2);
            
            clear A_row A_col
            % Initialization of A (describes the correlation of measurement error)
            A_col = zeros(n,1);
%             A_col(1) = 1;
%             A_col(2) = a1;
%             A_col(3) = a2;
            A_col(1:length(noiseCorrPar)) = noiseCorrPar;
            
            A_row(1) = 1;
            A_row(1,n) = 0;
            A = toeplitz(A_col,A_row);
            
            % White noise for measurement error
            if noiseCorrPar == 1
                A = eye(n);
            end
            
            % matrix F
            m = 2; %2 for double-integrated random walk, to model IG
            
            clear F_row F_col
            F_col = zeros(N,1);
            F_col(1) = 1;
            invF_col = F_col;
            for k=1:m
                F_col = filter([1 -1],1,F_col);
                invF_col = filter(ones(N,1),1,invF_col);
            end
            F_row(1) = 1;
            F_row(1,N) = 0;
            F = toeplitz(F_col,F_row); %Measures the roughness of the estimate
            
            B = inv(A'*A);
            H = B^(-1/2)*G*inv(F);
            [U,D,V] = svd(H);
            diagD = diag(D);
            DtraspD = diagD.^2;
            
            xi = U'*B^(-1/2)*(yWind-yWind(1));
            
            C = filter(invF_col,1,V);
            
            pfgammas = [];
            pfdiff = [];
            pfwrsss = [];
            parfor i = 0:7
                gammamin = 10^i;
                gammamax = 10^(i+1);
%             gammamin = par.gammamin;
%             gammamax = par.gammamax;
                max_iter = par.max_iter;
                toll = par.tol;
                
                itergamma = 1; %number of iterations
                reachedConvergence=0;
                while reachedConvergence==0
                   % update gamma
                    pfgamma = 10^((log10(gammamax)+log10(gammamin))/2);
                    % estimation in the new coordinates
                    eta = diagD./(DtraspD+pfgamma).*xi;
                    eta(n+1:N)=0;
                    % residuals in the new coordinates
                    ro = pfgamma*xi./(DtraspD+pfgamma);
                    % degrees of freedom in the new coordinates
                    qgamma = sum(diagD.^2./(DtraspD+pfgamma));
                    pfwrss = sum(ro.^2);
                    wess = sum((diagD.*xi./(DtraspD+pfgamma)).^2);
                
                    wrss_term = pfwrss/(N-qgamma);
                    wess_term = pfgamma*wess/qgamma;
                    difference = (wrss_term-wess_term); %Consistency criterion #3
                    diff = abs(difference);
                      
    %                 if showplots
    %                     figure(i+11)
    %                     plot(yWind-yWind(1)); hold on
    %                     plot(C*eta)
    %                     hold off
    %                     title([num2str(w) ', diff: ' num2str(diff), ', gamma: ' num2str(pfgamma)])
    %                     pause(0.1)
    %                 end
    
                    if difference<0
                        gammamax=pfgamma; % reduce regularity
                    else
                        gammamin=pfgamma; % increase regularity
                    end %if
                
                    if diff<toll || itergamma==max_iter
                        reachedConvergence=1; 
                        pfgammas = [pfgammas pfgamma];
                        pfdiff = [pfdiff diff];
                        pfwrsss = [pfwrsss pfwrss];
                    end %if
                    itergamma=itergamma+1;
                end
            end
            
%             figure
%             semilogy(0:7, pfdiff)
            
            gamma = pfgammas(pfdiff==min(pfdiff));
            
            wrss = pfwrsss(pfdiff == min(pfdiff));
            eta = (diagD.*xi)./(DtraspD+gamma);
            eta(n+1:N) = 0;
            qgamma = sum(DtraspD./(DtraspD+gamma));
            
            uHatw = C*eta;
    
            sigma2Hat=wrss/(N-qgamma); %estimated noise variance
            lambda2Hat=sigma2Hat/gamma;   %estimated model error variance
            
            covEta = sigma2Hat*ones(n,1)./(DtraspD+gamma);
            covEta(n+1:N,1)=lambda2Hat;
            
            c2 = C.^2;
            covuhatw = c2*covEta;
            
            if idxWindStart == 1 % first window
                uHat(w,:) = [nan(kstd,1); uHatw(round(wlen/2)+1:end)+yWind(1)];  %Estimated signal for first window
                covUHat(w,:) = [nan(kstd,1); covuhatw(round(wlen/2)+1:end)];  %Estimated covariance
            elseif idxWindEnd == length(hrin) % last window
                uHat(w,:) = [uHatw(round(wlen/2)+1:end)+yWind(1); nan(kstd,1)];  %Estimated signal for last window
                covUHat(w,:) = [ covuhatw(round(wlen/2)+1:end); nan(kstd,1)];  %Estimated covariance
            else
                uHat(w,:) = uHatw(round(wlen/2)+1:end)+yWind(1);  %Estimated signal 
                covUHat(w,:) = covuhatw(round(wlen/2)+1:end);  %Estimated covariance
            end
           
            if showplots
                figure(10)
                
                plot(window,'k.-','linewidth',1); hold on;
                plot(uHat(w,:),'r*','LineWidth',1.5); 
                %plot(uHat(w,:)+covUHat(w,:), 'b')
                plot(kstd+1, yWind(round(wlen/2)), 'd')
                hold off;
                legend window uHat windowStart
                title(w)
                pause
            end
    
        end
           
        sigma2 = [sigma2; sigma2Hat];
        lambda2 = [lambda2; lambda2Hat];
    
    end

    %% Filtered signal reconciliation

    %close 11 12
    uHatSmooth = NaN(nSample,1);
    stdSmooth = NaN(nSample,1);
    sigmaSmooth = NaN(nSample,1);
    lamdaSmooth = NaN(nSample,1);

    for w = 1:ceil(nSample/wlen)
        idxWindStart = idxWindStarts(w);
        idxWindEnd = idxWindStart + wlen - 1;
        uHatSmooth(idxWindStart:idxWindEnd) = uHat(w, kstd+1:wlen+kstd); % full actual window
        stdSmooth(idxWindStart:idxWindEnd) = covUHat(w, kstd+1:wlen+kstd);
        sigmaSmooth(idxWindStart:idxWindEnd) = ones(wlen,1)*sigma2(w);
        lambdaSmooth(idxWindStart:idxWindEnd) = ones(wlen,1)*lambda2(w);

        if showreconciliationplots
            % current window
            figure(11)
            subplot(2,1,1)
            plot(uHat(w,kstd:end))
            title(w)
            hold on
            
            subplot(2,1,2)
            plot(-kstd+1:wlen,uHat(w,1:wlen+kstd))
            title(w)
            hold on
        end
    
        
        if idxWindStart == 1 % first window
            uHatSmooth(idxWindEnd-kstd+1:idxWindEnd) = (kerWeigthPre.*uHat(w, wlen+1:wlen+kstd) + kerWeigthPost.*uHat(w+1,1:kstd))./(kerWeigthPost+kerWeigthPre); % last part of current window
            stdSmooth(idxWindEnd-kstd+1:idxWindEnd) = (kerWeigthPre.*covUHat(w, wlen+1:wlen+kstd) + kerWeigthPost.*covUHat(w+1,1:kstd))./(kerWeigthPost+kerWeigthPre); % last part of current window
            if showreconciliationplots         
                figure(11)
                %subplot(2,1,2)
                %plot(hr((w-1)*wlen+1:wlen*(w+1)))
                plot(wlen-kstd:2*wlen-1+kstd, uHat(w+1,:))
                plot(wlen-kstd+1:wlen+kstd, uHatSmooth((w)*wlen+1-kstd:wlen*(w)+kstd),'-*')

                legend uhat windownext smooth
                title(w)
                pause
            end
        elseif idxWindEnd == length(hrin) % last window
%             uHatSmooth(idxWindStart:idxWindStart+kstd-1) = (kerWeigthPre.*uHat(w-1, end-kstd+1:end) + kerWeigthPost.*uHat(w,kstd+1:2*kstd))./(kerWeigthPost+kerWeigthPre); % first part of current window
%             stdSmooth(idxWindStart:idxWindStart+kstd-1) = (kerWeigthPre.*covUHat(w-1, end-kstd+1:end) + kerWeigthPost.*covUHat(w,kstd+1:2*kstd))./(kerWeigthPost+kerWeigthPre); % first part of current window
%             
%             if showreconciliationplots
%                 subplot(2,1,1)
%                 plot(-wlen-kstd+1:kstd, uHat(w-1,:))
%                 plot(-kstd+1:kstd, uHatSmooth((w-1)*wlen+1-kstd:wlen*(w-1)+kstd),'-*')
%                 legend uhat windowpre smooth
%                 title(w)
%                 pause
%             end
        else
            uHatSmooth(idxWindEnd-kstd+1:idxWindEnd) = (kerWeigthPre.*uHat(w, wlen+1:wlen+kstd) + kerWeigthPost.*uHat(w+1,1:kstd))./(kerWeigthPost+kerWeigthPre); % last part of current window
            uHatSmooth(idxWindStart:idxWindStart+kstd-1) = (kerWeigthPre.*uHat(w-1, end-kstd+1:end) + kerWeigthPost.*uHat(w,kstd+1:2*kstd))./(kerWeigthPost+kerWeigthPre); % first part of current window
            stdSmooth(idxWindEnd-kstd+1:idxWindEnd) = (kerWeigthPre.*covUHat(w, wlen+1:wlen+kstd) + kerWeigthPost.*covUHat(w+1,1:kstd))./(kerWeigthPost+kerWeigthPre); % last part of current window
            stdSmooth(idxWindStart:idxWindStart+kstd-1) = (kerWeigthPre.*covUHat(w-1, end-kstd+1:end) + kerWeigthPost.*covUHat(w,kstd+1:2*kstd))./(kerWeigthPost+kerWeigthPre); % first part of current window
            if showreconciliationplots
                figure(11)
                subplot(2,1,2)
                %plot(hr((w-1)*wlen+1:wlen*(w+1)))
                plot(wlen-kstd:2*wlen-1+kstd, uHat(w+1,:))
                plot(wlen-kstd+1:wlen+1, uHatSmooth((w)*wlen-kstd+1:wlen*(w)+1),'-*')
                legend uhat windownext smooth
                title(w)
        
                subplot(2,1,1)
                plot(-wlen-kstd+1:kstd, uHat(w-1,:))
                plot(1:kstd, uHatSmooth((w-1)*wlen+1-0:wlen*(w-1)+kstd),'-*')
                legend uhat windowpre smooth
                title(w)
                pause
            end
        end
        if showreconciliationplots
            figure(11)
            subplot(2,1,1)
            hold off
            subplot(2,1,2)
            hold off
            
        end
    end
    
    if exist('missing', 'var') == 1
        uHatSmooth = uHatSmooth(1:end-missing);
        stdSmooth = stdSmooth(1:end-missing);
        sigmaSmooth = sigmaSmooth(1:end-missing);
    end

end

