function [thetaHat, rul]=BM(para0,weigh)
    clear global;
    global DegraUnit initDisPar TimeUnit time y thres ParamName thetaTrue signiLevel ns burnIn ny nt np

    % PROBLEM DEFINITION 1 (Required Variables)
    WorkName='Crack_BM';             
    DegraUnit = 'Crack size (m)';
    TimeUnit = 'Cycles';
    time = [0:100:3500]';
    y = [0.0100 0.0109 0.0101 0.0107 0.0110 0.0123 0.0099 0.0113 0.0132 0.0138 0.0148 0.0156 0.0155 0.0141 0.0169 0.0168]';                 
    thres=0.05;                  % threshold (critical value)
    ParamName = ['m'; 'C'; 's'];           % [npx1]: parameters' name to be estimated
    initDisPar=[4.0  0.2;  -23  1.1;  0.001  1e-5];          % [npx2]: prob. parameters of init./prior dist
    thetaTrue = [3.8; log(1.5e-10); 0.001];        % [npx1]: true values of parameters
    signiLevel=5;             % significance level for C.I. and P.I.
    ns=5e3;                     % number of particles/samples
    burnIn=0.2;                 % ratio for burn-in

    % % % PROGNOSIS using BM with MCMC
    ny=length(y); nt=ny;
    np=size(ParamName,1);
    sampl(:,1)=para0;        % Initial Samples
    [~,jPdf0]=MODEL(para0,time(1:ny));
    for i=2:ns/(1-burnIn)               % MCMC Process
        % proposal distribution (uniform)
        para1(:,1)=unifrnd(para0-weigh,para0+weigh);
        % (joint) PDF at new samples
        [~,jPdf1]=MODEL(para1,time(1:ny));
        % Acceptance/rejection criterion
        if rand<min(1,jPdf1/jPdf0)  
            para0=para1; jPdf0=jPdf1;
        end
        sampl(:,i)=para0;
    end
    nBurn=ns/(1-burnIn)-ns;              % Sampling Trace
    figure(4); for j=1:np; subplot(np,1,j); plot(sampl(j,:)); hold on
    plot([nBurn nBurn],[min(sampl(j,:)) max(sampl(j,:))],'r'); end
    thetaHat=sampl(:,nBurn+1:end);       % Final Sampling Results
    for k=1:length(time(ny:end))        % Degradation Prediction
        [zHat(k,:),~]=MODEL(thetaHat,time(ny-1+k));
        degraPredi(k,:)=normrnd(zHat(k,:),thetaHat(end,:));
    end

    % % % POST-PROCESSING
    degraTrue=[];                       % True Degradation
    if ~isempty(thetaTrue); degraTrue=MODEL(thetaTrue,time); end
    rul=POST(thetaHat,degraPredi,degraTrue); % RUL & Result Disp
    Name=[WorkName ' at ' num2str(time(ny)) '.mat']; save(Name);

end

function [z, poste]=MODEL(param, t)
    global y ParamName initDisPar ny np
    for j=1:np; eval([ParamName(j,:) ' =param(j,:);']); end
    %===== PROBLEM DEFINITION 2 (model equation) =====
    a0=0.01; dsig=75;
    z=(t.*exp(C).*(1-m./2).*(dsig*sqrt(pi)).^m ...
    +a0.^(1-m./2)).^(2./(2-m));
    loca=imag(z)~=0; z(loca)=1;
    %===========================================================
    if length(y)~=length(t); poste=[];
    else
        prior=prod(normpdf(param,initDisPar(:,1),initDisPar(:,2)));
        likel=1;
        for k=1:ny
            zeta=sqrt(log(1+(s./z(k)).^2)); 
            eta=log(z(k))-0.5*zeta.^2;
            likel=lognpdf(y(k),eta,zeta).*likel;
        end
        poste=likel.*prior; % Posterior
    end
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
