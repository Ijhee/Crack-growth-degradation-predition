function[thetaHat,rul]=PF
    clear global;
    global DegraUnit initDisPar TimeUnit time y thres ParamName thetaTrue signiLevel ns ny nt np
    %=== PROBLEM DEFINITION 1 (Required Variables) ============== 
    WorkName='Crack_PF';            % work results are saved by WorkName 
    DegraUnit = 'Crack size (m)';
    TimeUnit = 'Cycles';
    time = [0:100:3500]';
    y = [0.0100 0.0109 0.0101 0.0107 0.0110 0.0123 0.0099 0.0113 0.0132 0.0138 0.0148 0.0156 0.0155 0.0141 0.0169 0.0168]';                 
    thres=0.05;                          % threshold (critical value) 
    dt=20;
    ParamName = ['m '; 'C '; 's '; 'a '];   %[npx1]: parameters' name to be estimated 
    initDisPar=[4.0  0.2;  -23  1.1;  0.001  0;  0.01  0]; %[npx2]: prob. parameters of init./prior dist 
    thetaTrue=[3.8;log(1.5e-10); 0.001; 0.01]; %[npx1]: true values of parameters 
    signiLevel=5;             % significance level for C.I. and P.I.
    ns=5e3;                     % number of particles/samples
    %============================================================ 
    % % % PROGNOSISusing PF 
    ny=length(y); nt=ny; 
    np=size(ParamName,1); 
    for j=1:np                           %% Initial Distribution 
        param(j,:,1)=normrnd(initDisPar(j,1),initDisPar(j,2),1,ns); 
    end
    for k=2:length(time)          %% Update Process or Prognosis 
        % step1. prediction (prior) 
        paramPredi=param(:,:,k-1); 
        nk=(time(k)-time(k-1))/dt; 
    for k0=1:nk; paramPredi(np,:)=MODEL(paramPredi,dt); end 
    if k<=ny                                  % (Update Process) 
        % step2. update (likelihood) 
        mu=paramPredi(np,:); s=paramPredi(np-1,:);
        zeta=sqrt(log(1+(s./mu).^2)); eta=log(mu)-0.5*zeta.^2;
        likel=lognpdf(y(k),eta,zeta); 
        % step3. resampling 
        cdf=cumsum(likel)./sum(likel); 
        for i=1:ns
            u=rand; 
            loca=find(cdf>=u,1); 
            if isempty(loca)
                loca = 1; % or any other valid default value or handling strategy
            end
            param(:,i,k)=paramPredi(:,loca); 
        end

    else % (Prognosis) 
        param(:,:,k)=paramPredi; 
    end 
    end 
     thetaHat=param(1:np-1,:,ny);        %% Final Sampling Results 
     paramRearr=permute(param,[3 2 1]);  %% Degradation Prediction 
     zHat=paramRearr(:,:,np); 
     mu=zHat(ny:end,:); s=paramRearr(ny:end,:,np-1);
     zeta=sqrt(log(1+(s./mu).^2)); eta=log(mu)-0.5*zeta.^2;
     degraPredi=lognrnd(eta,zeta);
     % % % POST-PROCESSING 
     degraTrue=[];
     if~isempty(thetaTrue); k=1; 
        degraTrue0(1)=thetaTrue(np); degraTrue(1)=thetaTrue(np); 
        for k0=2:max(time)/dt+1
            degraTrue0(k0,1)=MODEL([thetaTrue(1:np-1); degraTrue0(k0-1)],dt); 
        loca=find((k0-1)*dt==time,1); 
        if ~isempty(loca); k=k+1; degraTrue(k)=degraTrue0(k0);end 
        end 
     end 
     rul=POST(thetaHat,degraPredi,degraTrue);%% RUL & Result Disp 
     Name = [WorkName ' at ' num2str(time(ny)) '.mat'];
     save(Name); 
end 

 function z1=MODEL(param,dt) 
    global ParamName np 
    for j=1:np; eval([ParamName(j,:) '=param(j,:);']); end 
     %===== PROBLEM DEFINITION 2 (model equation) =============== 
     dsig=75;
     z1=exp(C).*(dsig.*sqrt(pi*a)).^m.*dt+a;
     %=========================================================== 
 end
 
% [POST]: MATLAB Code for RUL Calculation and Results Plot
function rul = POST(thetaHat, degraPredi, degraTrue)
    global DegraUnit TimeUnit time y thres ParamName thetaTrue signiLevel ns ny nt
    np = size(thetaHat, 1);
    perceValue = [50 signiLevel 100 - signiLevel];

    figure(1);                      %% Distribution of Parameters
    for j = 1:np
        subplot(1, np, j);
        [frq, val] = hist(thetaHat(j, :), 30);
        bar(val, frq / ns / (val(2) - val(1)));
        xlabel(ParamName(j, :));
    end

    figure(2);                                %% Degradation Plot
    degraPI = prctile(degraPredi', perceValue)';
    f1(1, :) = plot(time(1:ny), y(:, 1), '.k', 'MarkerSize', 15); hold on;
    f1(2, :) = plot(time(ny:end), degraPI(:, 1), '--r', 'LineWidth', 1.3); % Median line
    f1(3:4, :) = plot(time(ny:end), degraPI(:, 2:3), ':r', 'LineWidth', 1.3); % 90% PI lines
    f2 = plot([0 time(end)], [thres thres], 'g');
    legend([f1(1:3, :); f2], 'Data', 'Median', [num2str(100 - 2 * signiLevel) '% PI'], 'Threshold')
    xlabel(TimeUnit); ylabel(DegraUnit);
    ylim([0 0.06]); % Set y-axis limit
    xlim([0 3500]); % Set x-axis limit
    xticks([1000 2000 3000]); % Set x-axis ticks

    i0 = 0;                                       %% RUL Prediction
    if y(nt(1)) - y(1) < 0; coeff = -1; else coeff = 1; end
    rul = zeros(1, ns);
    for i = 1:ns
        loca = find(degraPredi(:, i) * coeff >= thres * coeff, 1);
        if isempty(loca)
            i0 = i0 + 1;
            disp([num2str(i) 'th not reaching thres']);
        elseif loca == 1
            rul(i - i0) = 0;
        else
            rul(i - i0) = interp1([degraPredi(loca, i) degraPredi(loca - 1, i)], ...
                [time(ny - 1 + loca) time(ny - 2 + loca)], thres) - time(ny);
        end
    end

    rulPrct = prctile(rul, perceValue);
    figure(3);                             %% RUL Results Display
    [frq, val] = hist(rul, 30);
    bar(val, frq / ns / (val(2) - val(1))); hold on;
    xlabel(['RUL' ' (' TimeUnit ')']);
    titleName = ['at ' num2str(time(ny)) ' ' TimeUnit];
    title(titleName)
    fprintf(['\n Percentiles of RUL at %g ' TimeUnit], time(ny))
    fprintf('\n  %gth: %g,  50th (median): %g,  %gth: %g \n', ...
        perceValue(2), rulPrct(2), rulPrct(1), perceValue(3), rulPrct(3))

    if ~isempty(degraTrue)                  %% True Results Plot
        figure(1); % parameters
        for j = 1:np; subplot(1, np, j); hold on;
            plot(thetaTrue(j), 0, 'kp', 'markersize', 18);
        end

        figure(2); % degradation
        sl = 0; if ~isempty(nt); sl = ny - nt(end); end
        f3 = plot(time(sl + 1:sl + length(degraTrue)), degraTrue, 'k');
        legend([f1(1:3, :); f2; f3], 'Data', 'Median', ...
            [num2str(100 - 2 * signiLevel) '% PI'], 'Threshold', 'True')

        figure(3); % RUL
        loca = find(degraTrue * coeff >= thres * coeff, 1);
        rulTrue = interp1([degraTrue(loca) degraTrue(loca - 1)], ...
            [time(sl + loca) time(sl + loca - 1)], thres) - time(ny);
        plot(rulTrue, 0, 'kp', 'markersize', 18);
    end
end