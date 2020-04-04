if ~exist('pie')||~ishandle(pie.hFigure)
    cd ..
    launch_PIE;
%     cd CustomSimulations
    msgbox('Please change PIE experiment setup and run again');
    return;
end

%% stage uncertainty analysis
vals = linspace(0,25000,11); % Transform from 64-px detector to 1k detector
nTrials = 10;

% simulation parameters
pie.cb(pie.uibReset); % clear error sources
pie.uieMaxIteration.set(200);
pie.uieAccuracy.set(0);
pie.uieAlpha.set(0.5);
pie.uieBeta.set(0.5);
pie.uieSigma.set(0);

RMS = zeros(length(vals), nTrials);
tstart = datestr(now, 31);

for m = 1:nTrials
    for k = 1:length(vals)
        
        fprintf('Trial %d-%d/%d\n', k, m, length(vals));
        tic
                
        val = vals(k);
        pie.setSimParams('xyStageError',val);
        pie.simStackAndReconstruct();
        pie.uipSelectObject.setSelectedIndex(uint8(12));
        pie.cb(pie.uibAnalyze);
%         object = pie.dSelectedObject;
        rmsText = pie.uitRMS.get();
        RMS(k,m) = abs(str2double(rmsText(11:end)));

        fprintf('\n\n===========Simulation complete! ============\n\n');
        
        fprintf('Simulation took %s\n', s2f(toc) );
    end
end

fprintf('Simulation finished at %s\n\n', tstart);
fprintf('Simulation finished at %s\n\n', datestr(now, 31));

%% plot
% figure,h=plot(vals(1:end),mean(RMS,2));xlabel('Stage uncertainty (1\sigma) / nm'),ylabel('RMS object phase difference / rad');
% set(h,'LineWidth',2); set(gca,'FontSize',14);
figure,h=plot(vals(1:end)/(3e6)/pi*180,mean(RMS,2));xlabel('Scanning uncertainty (1\sigma) / degree'),ylabel('RMS object phase difference / rad');
set(h,'LineWidth',2); set(gca,'FontSize',14);

%% save data
save([pie.cAppPath,sprintf('/../../data/error analysis/xyStageError-study_%s.mat', regexprep(datestr(now, 31), ':', '.'))], 'vals', 'RMS');
