classdef PIE_Analyze < mic.Base
    
    
    properties (Constant)
        dWidth  = 1600;
        dHeight =  1000;
        
        % Axes tab IDS
        U8DATA              = 1
        U8PROBEOBJECT       = 2
        U8GUESS             = 3
        U8RECONSTRUCTION    = 4
        U8ANALYSIS          = 5
        U8LOG    = 6
        
        % Meta flags
        U8DOMAIN_REAL       = 1
        U8DOMAIN_FFT        = 2
        U8AUXPLOTOFF        = 3
        U8PLOTUNSHIFTED     = 4
        
        % Flags for reconstruction process:
        U8STATE_INITIAL             = 1
        U8STATE_DATA_LOADED         = 2
        U8STATE_PHASE_PROCESSED     = 3
        U8STATE_GUESS_PROCESSED     = 4
        U8STATE_RECONSTRUCTED       = 5
        U8STATE_ANALYZED       = 6
        
        U8MAXSTATES         = 6
        
        ceValidForAnalysis = {'Image capture'}; % list of valid scan types
        
        dColorActive        = [.8, 1, 1]
        dColorForbid        = [1, .8, .8]
        dColorDefault       = [.94, .94, .94]
    end
    
    properties
        
        cAppPath = fileparts(mfilename('fullpath'))
        
        % Graphical elements
        hFigure     % Main figure (not overwritable)
        
        u8State % Keeps track of where we are in reconstruction
        ceUIBStateControlled % List of state-controlled uibuttons
        
        % Main diffraction pattern(s)
        ceInt_0 % reset int values
        ceIntR = {} % processed (cropped, binned rotated)
        ceInt = {} % processed (cropped, binned rotated)
        
        
        tSelectedLog % table of log values
        ceImageLog % cell of log values that are image captures
        ceScanPaths
        ceSeriesPaths
        ceIntMeta = {} % meta data associated with images
        ceAnalysisPara={}
        ceAnalysisTable={}
        cSeries
        
        dAnalysisRegion = []
        dAnalysisRegion2 = []
        dAnalysisMask = []
        
        dReconstructed
        dZernike
        dZernikeResidual
        
        dResult
        dProbe
        dObject
        dProbeGuess
        dObjectGuess
        dProbeRecon
        dObjectRecon
        dSelectedObject
        dUnit_mm
        dError
        dPos_mm
        dLambda_um
        dCTF
        do_um
        dc_um
        dMaxObjectLen = 2000;
        
        
        % pupil guides
        dBeamWidthEstPx = 1
        dObstructionWidthEstPx
        
        %% Axes tab group
        uitgAxesDisplay     % displays axes: diffraction pattern, probe and pbject, pnitial guess, reconstruction, analysis, log
        
        % Axes: Data
        hsaInterferogram
        uibLeft
        uibRight
        uitStack
        uitMetaInfo
        uicbAutoCenter
        u8NumInt
        u8ActiveIntIdx
        
        
        % Axes: Probe and object
        haProbeAmp
        haProbePha
        haObjectAmp
        haObjectPha
        uicbRemoveLinearPhase
        
        % Axes: Initial guess
        haGuessProbeAmp
        haGuessProbePha
        haGuessObjectAmp
        haGuessObjectPha
        
        % Axes: Reconstruction
        haReconProbeAmp
        haReconProbePha
        haReconObjectAmp
        haReconObjectPha
        
        % Axes: Analysis
        haAnalysis
        haZernikeDecomp
        uitRMS
        
        % Axes: Log
        htLog
        
        %% Controls
        hpControls
        hpAnalysisSetup
        hpLoadInterferogram
        hpPhase
        hpAnalysis
        
        % Controls:Experiment Setup
        uieLambda
        uieNA
        uieScanRange
        uiez2
        uieGratTilt
        uieDetTilt
        prControlsSetup
        uieGlobalRot
        uieDetSize
        uieCenterObstruction
        uieBinning
        
        % Controls:Data
        hpData
        uitgSelectDataSource
        uicbFourierPtychography
        uieCCDCenter
        uieObsOffset
        uibSelectCenter
        uibResetCenter
        
        uipSelectMask
        uibLoadMask
        
        % Controls:Data:Load diffraction pattern single/stack
        uibLoadSingleFromLog
        uibLoadSingleFromImg
        uibLoadStackFromLog
        uipDataType
        uibLoadROI
        uipROIs
        
        % Keep these synchronized
        uieLogFileNameSingle
        uieLogFileNameStack
        
        uibSetLogFileSingle
        uibSetLogFileStack
        
        uilSingleList
        uilStackList
        
        % Controls:Data:Simulation
        uieRes
        uieZrn
        uiePhaseStepsSim
        uibLoadScanningPos
        uiez1
        uieNp % number of photos
        uibLoadZrn
        uibSimulate
        uibSimulateS
        uieScanSteps
        uipbExposureProgress
        
        
        % Controls:Data:Probe and Object
        uipProbeType
        uipObjectType
        uipPropagator
        uibLoadProbe
        uibLoadObject
        uibGenProbeObject
        uieRprobe
        uicbGuess
        uitOverlap
        uitSampling
        uieProbeOffset
        uieProbeAmp
        uibCopyProbe
        uibSimulatePO
        
        % Controls:Data:FPM
        uieNAo
        uieLo
        uieMag
        uieScanAngles
        
        % Controls:Data:SimSegments
        uieSegmentPath
        uibLoadSegment
        ceSegments = [];
        
        % Controls:Data:SimStochastics
        uibCustomSim
        uibReset
        
        uiePhaseShiftingError
        uieZLinearDrift
        uieXLinearDrift
        uieYLinearDrift
        
        uieGratingTiltErr
        uieDetectorTiltErr
        
        uieShotToShot
        uie2ndOrderStrength
        uie11orderStrength
        uieDetCurve
        uieFlareLevel
        uieMSFR
        uieNonlinearity
        uieAirflow
        
        % Controls:Reconstruction
        uitgAnalysisDomain
        uipUnwrapEngine
        uibComputePhase
        uibStop
        uieFDZ1
        uitGramRot
        uitIteration
        uieAlpha
        uieBeta
        uieMaxIteration
        uieAccuracy
        uicbGPU
        
        
        % Controls: Reconstruction: rPIE
        uieGamma
        uicbCorrectPos
        uipCorrectMethod
        uieModeNumber
        uilSelectMode
        uibShowMode
        
        % Controls: Reconstruction: RAAR
        uieDelta
        
        % Controls: Reconstruction: WDD
        uieEpsilon
        
        % Controls: Reconstruction: ML
        uipLikelihoodType
        uieRegularization
        
        % Controls:Analysis
        uipSelectRegion
        uipSelectObject
        uibAnalyze
        uibLoadAnalysisRegion
        uieSigma
        uibDataCursor
        
    end
    
    properties (SetAccess = private)
        
    end
    
    methods
        function this = PIE_Analyze()
            this.init()
        end
        
        
        function init(this)
            
            this.uitgAxesDisplay = ...
                mic.ui.common.Tabgroup('ceTabNames', {'Data', 'Probe and object','Initial guess', 'Reconstruction','Analysis','Log'});
            
            %% Axes:
            this.hsaInterferogram = mic.ui.axes.ScalableAxes(...
                'fhOnDomainChange', @(cDomain)this.replot(this.U8DATA,  cDomain));
            this.uibLeft        = mic.ui.common.Button('cText', '<', 'fhDirectCallback', @(src, evt)this.cb(src));
            this.uibRight       = mic.ui.common.Button('cText', '>', 'fhDirectCallback', @(src, evt)this.cb(src));
            this.uitStack       = mic.ui.common.Text('cVal', '0/0');
            this.uitMetaInfo    = mic.ui.common.Text('cVal', '');
            
            this.uitRMS                 = mic.ui.common.Text('cVal','RMS:');
            
            
            %% Controls: Experiment setup panel
            this.uieLambda      = mic.ui.common.Edit('cLabel', 'Lambda (nm)', 'cType', 'd','fhDirectCallback', @(src, evt)this.cb(src), 'lNotifyOnProgrammaticSet', false);
            this.uieScanRange   = mic.ui.common.Edit('cLabel', 'Scan range (mm)', 'cType', 'd', 'fhDirectCallback', @(src, evt)this.cb(src), 'lNotifyOnProgrammaticSet', false);
            this.uieNA          = mic.ui.common.Edit('cLabel', 'NA', 'cType', 'd', 'fhDirectCallback', @(src, evt)this.cb(src), 'lNotifyOnProgrammaticSet', false);
            this.uiez2          = mic.ui.common.Edit('cLabel', 'Z_2 (mm)', 'cType', 'd','fhDirectCallback', @(src, evt)this.cb(src), 'lNotifyOnProgrammaticSet', false);
            this.uieGratTilt    = mic.ui.common.Edit('cLabel', 'Gr. Tilt (deg)', 'cType', 'd');
            this.uieDetTilt     = mic.ui.common.Edit('cLabel', 'Dt. Tilt (deg)', 'cType', 'd');
            this.uieGlobalRot   = mic.ui.common.Edit('cLabel', 'CCD Rot (deg)', 'cType', 'd');
            this.uieDetSize     = mic.ui.common.Edit('cLabel', 'Dt. Size (mm)', 'cType', 'd','fhDirectCallback', @(src, evt)this.cb(src), 'lNotifyOnProgrammaticSet', false);
            this.uieCenterObstruction     = mic.ui.common.Edit('cLabel', 'Central.Obs', 'cType', 'd');
            this.uieBinning     = mic.ui.common.Edit('cLabel', 'Binning', 'cType', 'd', 'fhDirectCallback', @(src, evt)this.cb(src));
            
            
            this.uieLambda.set(13.5);
            this.uieScanRange.set(1);
            this.uieNA.set(0);
            this.uiez2.set(200);
            this.uieGratTilt.set(0);
            this.uieDetTilt.set(0);
            this.uieGlobalRot.set(0);
            this.uieDetSize.set(26);
            this.uieCenterObstruction.set(0);
            this.uieBinning.set(512);
            
            this.prControlsSetup = mic.ui.common.PositionRecaller(...
                'cConfigPath', fullfile(this.cAppPath, '+config'), ...
                'cName', 'config', ...
                'hGetCallback', @()this.prGetters(this.prControlsSetup), ...
                'hSetCallback', @(dRecall)this.prSetters(this.prControlsSetup, dRecall));
            
            %% Controls: Data panel
            this.uitgSelectDataSource = ...
                mic.ui.common.Tabgroup('ceTabNames', {'Load patterns', 'Load P/S Series', 'Simulation','Probe and object','FPM','Sim segments', 'Sim stochastics'});
            this.uieCCDCenter       = mic.ui.common.Edit('cLabel', 'CCD Center pixel', 'cType', 'c', ...
                'fhDirectCallback', @(src, evt)this.cb(src), 'lNotifyOnProgrammaticSet', false);
            this.uicbAutoCenter        = mic.ui.common.Checkbox('cLabel', 'Auto Center',  'fhDirectCallback', @(src, evt)this.cb(src));
            this.uieObsOffset       = mic.ui.common.Edit('cLabel', 'Obstruction offset pixel', 'cType', 'c', ...
                'fhDirectCallback', @(src, evt)this.cb(src), 'lNotifyOnProgrammaticSet', false);
            this.uibSelectCenter    = mic.ui.common.Button('cText', 'Select center', 'fhDirectCallback', @(src, evt)this.cb(src));
            this.uibResetCenter     = mic.ui.common.Button('cText', 'Set default center', 'fhDirectCallback', @(src, evt)this.cb(src));
            this.uicbFourierPtychography = mic.ui.common.Checkbox('cLabel', 'Fourier ptychography',  'fhDirectCallback', @(src, evt)this.cb(src));
            this.uipSelectMask      = mic.ui.common.Popup('cLabel', 'Select mask', 'ceOptions', {'Compute from det. geom.','Compute from MET5','Compute from eliptical geom.'}, ...
                'fhDirectCallback',@(src, evt)this.cb(src), 'lShowLabel', true);
            this.uibLoadMask        = mic.ui.common.Button('cText', 'Load mask', 'fhDirectCallback', @(src, evt)this.cb(src));
            
            this.uieCCDCenter.set('[]');
            this.uieObsOffset.set('[]');
            this.uicbAutoCenter.set(false);
            this.uicbFourierPtychography.set(false);
            
            % Controls:Data:LI/LIStack
            this.uieLogFileNameSingle   = mic.ui.common.Edit('cLabel', 'Log file name', 'fhDirectCallback', @(src, evt)this.cb(src));
            this.uieLogFileNameStack    = mic.ui.common.Edit('cLabel', 'Log file name', 'fhDirectCallback', @(src, evt)this.cb(src));
            
            this.uibSetLogFileSingle    = mic.ui.common.Button('cText', 'Set log file', 'fhDirectCallback', @(src, evt)this.cb(src));
            this.uibSetLogFileStack     = mic.ui.common.Button('cText', 'Set log file', 'fhDirectCallback', @(src, evt)this.cb(src));
            
            this.uibLoadSingleFromLog   = mic.ui.common.Button('cText', 'Load from data log', 'fhDirectCallback', @(src, evt)this.cb(src));
            this.uibLoadSingleFromImg   = mic.ui.common.Button('cText', 'Load from file', 'fhDirectCallback', @(src, evt)this.cb(src));
            this.uibLoadStackFromLog    = mic.ui.common.Button('cText', 'Load from data log', 'fhDirectCallback', @(src, evt)this.cb(src));
            
            this.uipDataType            = mic.ui.common.Popup('cLabel', 'Data format', 'ceOptions', {'2 x 1D (X,Y)', '2 x 1D (Y,X)', '2D'}, ...
                'fhDirectCallback',@(src, evt)this.cb(src), 'lShowLabel', true);
            this.uilSingleList          = mic.ui.common.List('cLabel', 'Logged images', ...
                'lShowDelete', false, 'lShowMove', false, 'lShowRefresh', false);
            
            this.uilStackList           = mic.ui.common.List('cLabel', 'Logged series', ...
                'lShowDelete', false, 'lShowMove', false, 'lShowRefresh', false);
            
            this.uipDataType.setSelectedIndex(uint8(3));
            
            % Controls:Data:simulation
            this.uieRes         = mic.ui.common.Edit('cLabel', 'Res', 'cType', 'd', 'fhDirectCallback', @(src, evt)this.cb(src), 'lNotifyOnProgrammaticSet', false);
            this.uiez1          = mic.ui.common.Edit('cLabel', 'Z_1 (mm)', 'cType', 'd');
            this.uieNp          = mic.ui.common.Edit('cLabel', 'N Phtn', 'cType', 'd');
            this.uieZrn         = mic.ui.common.Edit('cLabel', 'Zernike couples vector [N X 2]', 'cType', 'c','fhDirectCallback', @(src, evt)this.cb(src), 'lNotifyOnProgrammaticSet', false);
            
            
            this.uieScanSteps    = mic.ui.common.Edit('cLabel', 'N steps', 'cType', 'd', 'fhDirectCallback', @(src, evt)this.cb(src), 'lNotifyOnProgrammaticSet', false);
            this.uiePhaseStepsSim  = mic.ui.common.Edit('cLabel', 'Scanning positions', 'cType', 'c', 'fhDirectCallback', @(src, evt)this.cb(src), 'lNotifyOnProgrammaticSet', false);
            this.uibLoadZrn     = mic.ui.common.Button('cText', 'Load Zrn File', 'fhDirectCallback', @(src, evt)this.cb(src));
            this.uibSimulate    = mic.ui.common.Button('cText', 'Simulate Single', 'fhDirectCallback', @(src, evt)this.cb(src));
            this.uibSimulateS   = mic.ui.common.Button('cText', 'Simulate Stack', 'fhDirectCallback', @(src, evt)this.cb(src));
            
            this.uibLoadScanningPos...
                = mic.ui.common.Button('cText', 'Load scanning pos', 'fhDirectCallback', @(src, evt)this.cb(src));
            this.uipbExposureProgress = mic.ui.common.ProgressBar(...
                'dColorFill', [.4, .4, .8], ...
                'dColorBg', [1, 1, 1], ...
                'dHeight', 15, ...
                'dWidth', 455);
            
            this.uiez1.set(1);
            this.uieNp.set(0);
            this.uieRes.set(256);
            this.uieScanSteps.set(4);
            dN = this.uieScanSteps.get();
            dL = this.uieScanRange.get();
            dPosString = sprintf('%0.6f:%0.6f/%d:%0.6f', -dL/2,dL, dN-1, dL/2);
            this.uiePhaseStepsSim.set(sprintf('[%s;%s]''', dPosString, dPosString));
            
            this.uieZrn.set('[]');
            
            
            % Controls:Data:Probe and object
            this.uipProbeType     = mic.ui.common.Popup('cLabel', 'Probe type', 'ceOptions', {'Defocus wave','Plane wave','From reconstruction'}, ...
                'fhDirectCallback',@(src, evt)this.cb(src), 'lShowLabel', true);
            this.uipObjectType     = mic.ui.common.Popup('cLabel', 'Object type', 'ceOptions', {'Vacuum','Cameraman','Low resolution',...
                'Load object','From reconstruction'}, ...
                'fhDirectCallback',@(src, evt)this.cb(src), 'lShowLabel', true);
            this.uipPropagator     = mic.ui.common.Popup('cLabel', 'Propagator', 'ceOptions', {'angular spectrum','fourier','fresnel'}, ...
                'fhDirectCallback',@(src, evt)this.cb(src), 'lShowLabel', true);
            this.uibLoadProbe    = mic.ui.common.Button('cText', 'Load probe', 'fhDirectCallback', @(src, evt)this.cb(src));
            this.uibLoadObject    = mic.ui.common.Button('cText', 'Load object', 'fhDirectCallback', @(src, evt)this.cb(src));
            this.uibGenProbeObject   = mic.ui.common.Button('cText', 'Generate', 'fhDirectCallback', @(src, evt)this.cb(src));
            this.uibCopyProbe   = mic.ui.common.Button('cText', 'Copy probe', 'fhDirectCallback', @(src, evt)this.cb(src));
            this.uieRprobe = mic.ui.common.Edit('cLabel', 'Probe radius on det (mm)', 'cType', 'd', 'fhDirectCallback', @(src, evt)this.cb(src), 'lNotifyOnProgrammaticSet', false);
            this.uieRprobe.set(this.uieDetSize.get()/4);
            this.uicbGuess = mic.ui.common.Checkbox('cLabel', 'Initial guess',  'fhDirectCallback', @(src, evt)this.cb(src));
            this.uicbGuess.set(true);
            this.uitOverlap                 = mic.ui.common.Text('cVal','Overlap:');
            this.uitSampling                 = mic.ui.common.Text('cVal','Sampling factor:');
            this.uieProbeOffset       = mic.ui.common.Edit('cLabel', 'Probe offset (mm)', 'cType', 'c', ...
                'fhDirectCallback', @(src, evt)this.cb(src), 'lNotifyOnProgrammaticSet', false);
            this.uieProbeAmp       = mic.ui.common.Edit('cLabel', 'Probe amplitude', 'cType', 'd', ...
                'fhDirectCallback', @(src, evt)this.cb(src), 'lNotifyOnProgrammaticSet', false);
            this.uibSimulatePO   = mic.ui.common.Button('cText', 'Simulate', 'fhDirectCallback', @(src, evt)this.cb(src));
            
            this.uieProbeOffset.set('[]');
            this.uieProbeAmp.set(1);
            this.uipPropagator.setSelectedIndex(uint8(2));
            % Controls:Data:FPM
            this.uieNAo = mic.ui.common.Edit('cLabel', 'Object NA', 'cType', 'd', 'fhDirectCallback', @(src, evt)this.cb(src),'lNotifyOnProgrammaticSet', false);
            this.uieLo = mic.ui.common.Edit('cLabel', 'Obj. Length (mm)', 'cType', 'd', 'fhDirectCallback', @(src, evt)this.cb(src),'lNotifyOnProgrammaticSet', false);
            this.uieMag = mic.ui.common.Edit('cLabel', 'Magnification', 'cType', 'd', 'fhDirectCallback', @(src, evt)this.cb(src),'lNotifyOnProgrammaticSet', false);
            this.uieNAo.set(0);
            this.uieLo.set(3);
            
            this.uieScanAngles  = mic.ui.common.Edit('cLabel', 'Scanning angles', 'cType', 'c', 'fhDirectCallback', @(src, evt)this.cb(src), 'lNotifyOnProgrammaticSet', false);
            %             dL = this.uieScanRange.get();
            Lo_mm = this.uieLo.get();
            this.uieMag.set(1);
            dA =atan(eval(this.uiePhaseStepsSim.get())/Lo_mm)/pi*180;
            this.uieScanAngles.set(mat2str(dA));
            
            % Controls:Data:Sim segments
            this.uieSegmentPath    = mic.ui.common.Edit('cLabel', 'Segment path', 'fhDirectCallback', @(src, evt)this.cb(src));
            this.uibLoadSegment     = mic.ui.common.Button('cText', 'Load segments', 'fhDirectCallback', @(src, evt)this.cb(src));
            
            % Controls:Data:Sim stochastics
            this.uibCustomSim = mic.ui.common.Button('cText', 'Custom Sim', 'fhDirectCallback', @(src, evt)this.cb(src));
            this.uibReset = mic.ui.common.Button('cText', 'Reset', 'fhDirectCallback', @(src, evt)this.cb(src));
            
            this.uiePhaseShiftingError = mic.ui.common.Edit('cLabel', 'Stage uncertainty (nm)', 'cType', 'd', 'fhDirectCallback', @(src, evt)this.cb(src), 'lNotifyOnProgrammaticSet', false);
            
            this.uieZLinearDrift    = mic.ui.common.Edit('cLabel', 'Z Lin Drift (nm)', 'cType', 'd', 'fhDirectCallback', @(src, evt)this.cb(src), 'lNotifyOnProgrammaticSet', false);
            this.uieXLinearDrift   = mic.ui.common.Edit('cLabel', 'X Lin Drift (nm)', 'cType', 'd', 'fhDirectCallback', @(src, evt)this.cb(src), 'lNotifyOnProgrammaticSet', false);
            this.uieYLinearDrift   = mic.ui.common.Edit('cLabel', 'Y Lin Drift (nm)', 'cType', 'd', 'fhDirectCallback', @(src, evt)this.cb(src), 'lNotifyOnProgrammaticSet', false);
            this.uieGratingTiltErr  = mic.ui.common.Edit('cLabel', '-Grat tilt error(deg)', 'cType', 'd', 'fhDirectCallback', @(src, evt)this.cb(src), 'lNotifyOnProgrammaticSet', false);
            this.uieDetectorTiltErr = mic.ui.common.Edit('cLabel', '-Det tilt error (deg)', 'cType', 'd', 'fhDirectCallback', @(src, evt)this.cb(src), 'lNotifyOnProgrammaticSet', false);
            this.uieShotToShot       = mic.ui.common.Edit('cLabel', 'Shot to shot (%)', 'cType', 'd', 'fhDirectCallback', @(src, evt)this.cb(src), 'lNotifyOnProgrammaticSet', false);
            this.uie2ndOrderStrength = mic.ui.common.Edit('cLabel', '-', 'cType', 'd', 'fhDirectCallback', @(src, evt)this.cb(src), 'lNotifyOnProgrammaticSet', false);
            this.uie11orderStrength = mic.ui.common.Edit('cLabel', '-', 'cType', 'd', 'fhDirectCallback', @(src, evt)this.cb(src), 'lNotifyOnProgrammaticSet', false);
            this.uieDetCurve       = mic.ui.common.Edit('cLabel', '-Dt. Curve(um)', 'cType', 'd', 'fhDirectCallback', @(src, evt)this.cb(src), 'lNotifyOnProgrammaticSet', false);
            this.uieFlareLevel       = mic.ui.common.Edit('cLabel', 'DC flare', 'cType', 'd', 'fhDirectCallback', @(src, evt)this.cb(src), 'lNotifyOnProgrammaticSet', false);
            this.uieMSFR       = mic.ui.common.Edit('cLabel', 'MSFR', 'cType', 'd', 'fhDirectCallback', @(src, evt)this.cb(src), 'lNotifyOnProgrammaticSet', false);
            this.uieNonlinearity       = mic.ui.common.Edit('cLabel', 'Nonlinearity', 'cType', 'd', 'fhDirectCallback', @(src, evt)this.cb(src), 'lNotifyOnProgrammaticSet', false);
            this.uieAirflow       = mic.ui.common.Edit('cLabel', 'Airflow', 'cType', 'd', 'fhDirectCallback', @(src, evt)this.cb(src), 'lNotifyOnProgrammaticSet', false);
            
            
            %% Controls: Reconstruction
            this.uitIteration               = mic.ui.common.Text('cVal','');
            this.uitgAnalysisDomain = mic.ui.common.Tabgroup('ceTabNames', {'rPIE', 'RAAR', 'WDD','ML'});
            this.uipUnwrapEngine        = mic.ui.common.Popup('cLabel', 'Unwrapping algorithm', 'ceOptions', {'Sorting reliability unwrap'}, ...
                'fhDirectCallback',@(src, evt)this.cb(src), 'lShowLabel', true);
            this.uibComputePhase        = mic.ui.common.Button('cText', 'Reconstruction',  'fhDirectCallback', @(src, evt)this.cb(src));
            this.uibStop        = mic.ui.common.Button('cText', 'Stop',  'fhDirectCallback', @(src, evt)this.cb(src));
            this.uieAlpha                = mic.ui.common.Edit('cLabel', 'Alpha (O)', 'cType', 'd', 'fhDirectCallback', @(src, evt)this.cb(src));
            this.uieBeta                = mic.ui.common.Edit('cLabel', 'Beta (P)', 'cType', 'd', 'fhDirectCallback', @(src, evt)this.cb(src));
            this.uieMaxIteration              = mic.ui.common.Edit('cLabel', 'Max iteration', 'cType', 'd', 'fhDirectCallback', @(src, evt)this.cb(src));
            this.uieAccuracy              = mic.ui.common.Edit('cLabel', 'Accuracy', 'cType', 'd', 'fhDirectCallback', @(src, evt)this.cb(src));
            this.uicbGPU = mic.ui.common.Checkbox('cLabel', 'GPU acceleration',  'fhDirectCallback', @(src, evt)this.cb(src));
            this.uicbGPU.set(false);
            this.uieAlpha.set(0.5);
            this.uieBeta.set(0.5);
            this.uieMaxIteration.set(1000);
            this.uieAccuracy.set(0.0001);
            
            % Controls: Reconstruction: rPIE
            this.uieGamma                = mic.ui.common.Edit('cLabel', 'Gamma', 'cType', 'd', 'fhDirectCallback', @(src, evt)this.cb(src));
            this.uieGamma.set(0.2);
            this.uicbCorrectPos = mic.ui.common.Checkbox('cLabel', 'Correct position',  'fhDirectCallback', @(src, evt)this.cb(src));
            this.uicbCorrectPos.set(false);
            this.uieModeNumber = mic.ui.common.Edit('cLabel', 'Mode Number', 'cType', 'd', 'fhDirectCallback', @(src, evt)this.cb(src));
            this.uieModeNumber.set(1);
            this.uipCorrectMethod     = mic.ui.common.Popup('cLabel', 'Correct method', 'ceOptions', {'pcPIE','e-pcPIE'}, ...
                'fhDirectCallback',@(src, evt)this.cb(src), 'lShowLabel', true);
            this.uilSelectMode =    mic.ui.common.List('cLabel', 'Select mode', ...
                'lShowDelete', false, 'lShowMove', false, 'lShowRefresh', false);
            this.uibShowMode        = mic.ui.common.Button('cText', 'Show',  'fhDirectCallback', @(src, evt)this.cb(src));
            this.uilSelectMode.setOptions({'Mode 1'});
            
            % Controls: Reconstruction: RAAR
            this.uieDelta                = mic.ui.common.Edit('cLabel', 'Delta', 'cType', 'd', 'fhDirectCallback', @(src, evt)this.cb(src));
            this.uieDelta.set(0.1);
            
            % Controls: Reconstruction: WDD
            this.uieEpsilon                = mic.ui.common.Edit('cLabel', '\epsilon', 'cType', 'd', 'fhDirectCallback', @(src, evt)this.cb(src));
            this.uieEpsilon.set(0.1);
            
            % Controls: Reconstruction: ML
            this.uipLikelihoodType     = mic.ui.common.Popup('cLabel', 'Likelihood Type', 'ceOptions', {'Poisson','amplitude'}, ...
                'fhDirectCallback',@(src, evt)this.cb(src), 'lShowLabel', true);
            this.uieRegularization                = mic.ui.common.Edit('cLabel', 'Regularization', 'cType', 'd', 'fhDirectCallback', @(src, evt)this.cb(src));
            this.uieRegularization.set(1);
            
            %% Controls:Analysis
            this.uipSelectObject      = mic.ui.common.Popup('cLabel', 'Select analysis object', 'ceOptions',...
                {'Scanning position','Residual errors','Object amplitude','Object phase','Probe amplitude','Probe phase',...
                'Object spectrum amplitude','Object spectrum phase','Probe spectrum amplitude','Probe spectrum phase',...
                'Object amplitude difference','Object phase difference','Probe amplitude difference','Probe phase difference'}, ...
                'fhDirectCallback',@(src, evt)this.cb(src), 'lShowLabel', true);
            this.uipSelectRegion      = mic.ui.common.Popup('cLabel', 'Select analysis region', 'ceOptions', ...
                {'Entire region','Compute from scanning positions','Compute from probe diameter','Load region'}, ...
                'fhDirectCallback',@(src, evt)this.cb(src), 'lShowLabel', true);
            this.uibAnalyze        = mic.ui.common.Button('cText', 'Analyze', 'fhDirectCallback', @(src, evt)this.cb(src));
            this.uibLoadAnalysisRegion        = mic.ui.common.Button('cText', 'Load analysis region', 'fhDirectCallback', @(src, evt)this.cb(src));
            this.uibDataCursor        = mic.ui.common.Toggle('cTextFalse', 'Data cursor off','cTextTrue', 'Data cursor on', 'fhDirectCallback', @(src, evt)this.cb(src));
            this.uieSigma = mic.ui.common.Edit('cLabel', 'Sigma filter', 'cType', 'd', 'fhDirectCallback', @(src, evt)this.cb(src));
            this.uieSigma.set(0);
            
            %% Add buttons to state flag lists
            this.ceUIBStateControlled = {};
            
            this.ceUIBStateControlled{this.U8STATE_INITIAL} = ...
                {this.uibSimulate, this.uibSimulateS, this.uibLoadSingleFromLog, ...
                this.uibLoadSingleFromImg, this.uibLoadStackFromLog};
            this.ceUIBStateControlled{this.U8STATE_DATA_LOADED} = ...
                {this.uibComputePhase};
            this.ceUIBStateControlled{this.U8STATE_PHASE_PROCESSED} = ...
                {};
            this.ceUIBStateControlled{this.U8STATE_RECONSTRUCTED} = ...
                {};
            this.ceUIBStateControlled{this.U8STATE_ANALYZED} = ...
                {};
            
            % Set state:
            this.setState(this.U8STATE_INITIAL);
        end
        
        % Callback handler
        function cb(this, src,evt)
            switch src
                %% Experiment setup
                case this.uieNA
                    NA = this.uieNA.get();
                    if NA>0
                        this.uipProbeType.setSelectedIndex(uint8(1));
                        z_mm = this.uiez2.get();
                        df_mm = this.uiez1.get();
                        Rc_mm = (z_mm+df_mm)*tan(asin(NA));
                        this.uieRprobe.set(Rc_mm);
                    else
                        this.uipProbeType.setSelectedIndex(uint8(2));
                    end
                    
                    if this.uicbFourierPtychography.get()
                        Lo_mm = this.uieLo.get();
                        Li_mm = this.uiez2.get();
                        NAi = this.uieNA.get();
                        NAo = sin(atan(Li_mm/Lo_mm*tan(asin(NAi))));
                        this.uieNAo.set(NAo);
                        lambda_mm =this.uieLambda.get()/1000000;
                        Rc_mm = 0.61*lambda_mm/NAi;
                        this.uieRprobe.set(Rc_mm);
                        this.uieMag.set(NAo/NAi);
                    end
                    
                    % redraw guide lines
                    try
                        this.replot(this.U8DATA, []);
                    catch
                    end
                case this.uieBinning
                    try
                        N = this.uieBinning.get();
                        this.uieRes.set(N);
                        propagator = this.uipPropagator.getOptions{this.uipPropagator.getSelectedIndex()};
                        lambda_um   = this.uieLambda.get()/1000;
                        z_um       = this.uiez2.get()*1000;
                        detSize_um  = this.uieDetSize.get()*1000;
                        this.dc_um = detSize_um/N;
                        u8ModeId = this.uilSelectMode.getSelectedIndexes();
                        if strcmp(propagator,'fourier')
                            this.do_um(u8ModeId) = lambda_um*z_um/N/this.dc_um;
                        else
                            this.do_um(u8ModeId) = this.dc_um; % object pixel pitch
                        end
                    catch
                    end
                    
                case {this.uiez2, this.uieLambda, this.uieDetSize,this.uipPropagator}
                    propagator = this.uipPropagator.getOptions{this.uipPropagator.getSelectedIndex()};
                    N           = this.uieRes.get();
                    lambda_um   = this.uieLambda.get()/1000;
                    z_um       = this.uiez2.get()*1000;
                    detSize_um  = this.uieDetSize.get()*1000;
                    this.dc_um = detSize_um/N;
                    u8ModeId = this.uilSelectMode.getSelectedIndexes();
                    if strcmp(propagator,'fourier')
                        this.do_um(u8ModeId) = lambda_um*z_um/N/this.dc_um;
                    else
                        this.do_um(u8ModeId) = this.dc_um; % object pixel pitch
                    end
                    % redraw guide lines
                    this.replot(this.U8DATA, []);
                case this.uieCenterObstruction
                    % redraw guide lines
                    this.replot(this.U8DATA, []);
                    
                case this.uieScanRange
                    dN = this.uieScanSteps.get();
                    dL = this.uieScanRange.get();
                    dPosString = sprintf('%0.6f:%0.6f/%d:%0.6f', -dL/2,dL, dN-1, dL/2);
                    dPos = eval(dPosString);
                    %                     k = 1;
                    %                     for i =1:dN
                    %                         for j =1:dN
                    %                             dPos2(k,:)=[dPos(i),dPos(j)];
                    %                             k=k+1;
                    %                         end
                    %                     end
                    %                     this.dPos_mm = dPos2;
                    [dm,dn] = meshgrid(dPos,dPos);
                    this.dPos_mm  = [dm(:),dn(:)];
                    this.uiePhaseStepsSim.set(sprintf('[%s;%s]''', dPosString, dPosString));
                    Lo_mm = this.uieLo.get();
                    dA =atan(eval(this.uiePhaseStepsSim.get())/Lo_mm)/pi*180;
                    this.uieScanAngles.set(mat2str(dA));
                    
                    if strcmp(this.uipSelectObject.getOptions{this.uipSelectObject.getSelectedIndex()},...
                            'Scanning position')&&strcmp(this.uitgAxesDisplay.getSelectedTabName(),'Analysis')
                        % Plot wavefronts on phase tab
                        this.replot(this.U8ANALYSIS, []);
                    end
                    
                    %% Data
                case this.uieObsOffset
                    this.validateCouplesEditBox(src, '[]');
                    % redraw guide lines
                    this.replot(this.U8DATA, []);
                    % Need to recenter data
                    this.handleLoadData();
                    
                case this.uibSelectCenter
                    [dSr, dSc] = size(this.ceInt{1});
                    
                    % hide guide lines
                    this.replot(this.U8DATA, this.U8AUXPLOTOFF);
                    drawnow
                    dPts = this.hsaInterferogram.ginput(this.dObstructionWidthEstPx);
                    this.uieCCDCenter.set(sprintf('[%d, %d]', round(dPts(1)), round(dPts(2))));
                    
                    
                    this.handleLoadData(this.ceInt_0,{'shifted'});
                    
                    % redraw guide lines
                    this.replot(this.U8DATA, []);
                    
                case this.uibResetCenter
                    this.uieCCDCenter.set('[]');
                    this.handleLoadData();
                    this.replot(this.U8DATA, []);% redraw guide lines
                    
                case this.uipSelectMask
                    SelectMask=this.uipSelectMask.getSelectedIndex();
                    switch SelectMask
                        case 1
                            this.computeDefaultMaskGeometry();
                        case 2
                            this.computeMET5Mask();
                        case 3
                            this.computeElipticalMask();
                    end
                    
                case this.uibLoadMask
                    
                case this.uieCCDCenter
                    this.validateCouplesEditBox(src, '[]');
                    % redraw guide lines
                    this.replot(this.U8DATA, []);
                    % Need to recenter data
                    this.handleLoadData();
                    
                    
                    % Load patterns
                case this.uibLoadSingleFromImg
                    cDataDir = fullfile(this.cAppPath, '..', '..', 'data','Aerial images');
                    [fn,pn]=uigetfile({'*.mat;*.SPE','PIE int bundle (*.mat) or WinView (*.SPE)'},'Loading',cDataDir);
                    fileformat=fn(end-2:end);
                    filename= strcat(pn,fn);
                    switch fileformat
                        case 'SPE'
                            ceInts=lsianalyze.utils.speread(filename);
                            this.handleLoadData(ceInts, {});
                            
                        case 'mat'
                            load(filename);
                            try
                                this.dPos_mm = dPos_mm;
                                L = length(aerialImages);
                                N = ceil(sqrt(L));
                                if N^2>L
                                    for i= L+1:N^2
                                        aerialImages{i} = zeros(size(aerialImages{1}));
                                    end
                                end
                                this.handleLoadData(reshape(aerialImages,N,N), {});
                            catch
                            end
                    end
                    
                case this.uibLoadSingleFromLog
                    u8loadIdx = this.uilSingleList.getSelectedIndexes();
                    load(this.ceScanPaths{u8loadIdx});
                    this.handleLoadData({dImg}, {stLog});
                    
                    % Load P/S Sries
                case {this.uibSetLogFileSingle, this.uibSetLogFileStack}
                    cDataDir = fullfile(this.cAppPath, '..', '..', '..', 'Data', '*.csv');
                    [d, p] = uigetfile(cDataDir);
                    this.uieLogFileNameSingle.set([p d]);
                    this.uieLogFileNameStack.set([p d]);
                    % Load Analysis log
                    %                     try
                    if strcmp( d(end-10:end),'scanlog.csv')
                        this.uitgAxesDisplay.selectTabByIndex(this.U8LOG);
                        %this.replot(this.U8LOG, []);
                        filename=[p d(1:end-11)];
                        filename=strcat(filename,'analysislog.csv');
                        [fid, isNewLogFile] = this.openOrCreateFile( filename);
                        fclose(fid);
                        if isNewLogFile
                            delete(filename);
                        else
                            ceLog= readtable(filename);
                            ceTable = table2cell(ceLog);
                            ceAnalysisTables=ceTable(:,32:end);
                            [M,~]=size(ceAnalysisTables);
                            for i=1:M
                                ceAnalysisTables{i,1}=char(ceAnalysisTables{i,1});
                            end
                            set(this.htLog,'Data',ceAnalysisTables);
                            this.ceAnalysisTable=ceAnalysisTables;
                        end
                    end
                    %                     catch
                    %                     end
                    % Load csv file and populate list boxes
                    this.parseLogFile([p d], p);
                    
                case {this.uieLogFileNameSingle, this.uieLogFileNameStack}
                    cPath = src.get();
                    
                    cDir = [fileparts(cPath) filesep];
                    this.parseLogFile(cPath, cDir);
                    
                    
                    
                case this.uibLoadStackFromLog
                    u8loadIdx = this.uilStackList.getSelectedIndexes();
                    % save selected Series
                    s=this.uilStackList.getOptions();
                    if length(u8loadIdx)==2
                        this.cSeries=strcat(s{u8loadIdx(1)}(findstr(s{u8loadIdx(1)},'Series'):findstr(s{u8loadIdx(1)},',')-1),...
                            '&',s{u8loadIdx(2)}(findstr(s{u8loadIdx(2)},'Series')+6:findstr(s{u8loadIdx(1)},',')-1));
                    else
                        this.cSeries=s{u8loadIdx(1)}(findstr(s{u8loadIdx(1)},'Series'):findstr(s{u8loadIdx(1)},',')-1);
                    end
                    ceInts = cell(0);
                    ceMeta = cell(0);
                    
                    
                    switch double(this.uipDataType.getSelectedIndex())
                        case {1, 2} % 2 X 1D
                            if numel(u8loadIdx) ~= 2
                                msgbox('Need 2 data sets for 2 X 1D option');
                                
                                return
                            end
                            for k = 1:length(this.ceSeriesPaths{u8loadIdx(1)})
                                fprintf('Loading X: %d of %d\n', k, length(this.ceSeriesPaths{u8loadIdx(1)}));
                                load(this.ceSeriesPaths{u8loadIdx(1)}{k});
                                
                                ceInts{k,1} = dImg;
                                ceMeta{k,1} = stLog;
                            end
                            for k = 1:length(this.ceSeriesPaths{u8loadIdx(2)})
                                fprintf('Loading Y: %d of %d\n', k, length(this.ceSeriesPaths{u8loadIdx(2)}));
                                load(this.ceSeriesPaths{u8loadIdx(2)}{k});
                                
                                ceInts{k,2} = dImg;
                                ceMeta{k,2} = stLog;
                            end
                            
                            dN = size(ceInts,1);
                            
                            %                             % crop off redundant value
                            %
                            %                             ceInts = ceInts(1:dN-21, :);
                            %                             ceMeta = ceMeta(1:dN-21, :);
                            
                            if double(this.uipDataType.getSelectedIndex()) == 2
                                ceInts = fliplr(ceInts);
                                ceMeta = fliplr(ceMeta);
                            end
                            
                            this.handleLoadData(ceInts, ceMeta);
                            
                        case 3 % 2D
                            
                            for k = 1:length(this.ceSeriesPaths{u8loadIdx})
                                fprintf('Loading %d of %d\n', k, length(this.ceSeriesPaths{u8loadIdx}));
                                load(this.ceSeriesPaths{u8loadIdx}{k});
                                
                                ceInts{k} = dImg;
                                ceMeta{k} = stLog;
                            end
                            
                            
                            % For now, let's stack images into a square:
                            dN = sqrt(length(ceInts));
                            ceInts = reshape(ceInts, dN, dN);
                            ceMeta = reshape(ceMeta, dN, dN);
                            
                            % These images are rastered, so we need to flip every other one:
                            
                            for k = 2:2:dN % New way after flipping ret fine dir
                                %for k = 1:2:dN % Old way (before 3/9/18)
                                ceInts(:, k) =  flipud(ceInts(:, k));
                                ceMeta(:, k) =  flipud(ceMeta(:, k));
                            end
                            % Now crop off redundant image
                            ceInts = ceInts(1:dN-1,1:dN-1);
                            ceMeta = ceMeta(1:dN-1,1:dN-1);
                            this.handleLoadData(ceInts, ceMeta);
                    end
                    if this.uicbAutoCenter.get()
                        dCenter = lsianalyze.utils.autoCenter(this.ceInt);
                        this.uieCCDCenter.set(['[',num2str(dCenter),']']);
                        % redraw guide lines
                        this.replot(this.U8DATA, []);
                        % Need to recenter data
                        this.handleLoadData();
                    end
                    
                    
                    % Simulation
                case this.uiePhaseStepsSim
                    [lValid, vals] = this.validateCouplesEditBox(src, '[]');
                    if lValid
                        [sr, sc] = size(vals);
                        this.uieScanSteps.set(sr);
                        dN = this.uieScanSteps.get();
                        Lo_mm = this.uieLo.get();
                        dA =atan(eval(this.uiePhaseStepsSim.get())/Lo_mm)/pi*180;
                        this.uieScanAngles.set(mat2str(dA));
                    end
                    
                case this.uieRes
                    N = this.uieRes.get();
                    this.uieBinning.set(N);
                    propagator = this.uipPropagator.getOptions{this.uipPropagator.getSelectedIndex()};
                    lambda_um   = this.uieLambda.get()/1000;
                    z_um       = this.uiez2.get()*1000;
                    detSize_um  = this.uieDetSize.get()*1000;
                    this.dc_um = detSize_um/N;
                    u8ModeId = this.uilSelectMode.getSelectedIndexes();
                    if strcmp(propagator,'fourier')
                        this.do_um(u8ModeId) = lambda_um*z_um/N/this.dc_um;
                    else
                        this.do_um(u8ModeId) = this.dc_um; % object pixel pitch
                    end
                    
                case this.uibLoadZrn
                    [p, d] = uigetfile();
                    zrn = load([d,p]);
                    this.setZrnString(zrn.zrn);
                    
                case this.uibSimulateS
                    this.simulateInteferograms(true);
                    
                case this.uibSimulate
                    this.simulateInteferograms(false);
                    
                case this.uibLoadScanningPos
                    cDataDir = fullfile(this.cAppPath,  '..', '..', 'Data','scanning');
                    [fn,pn]=uigetfile({'*.mat','Position (*.mat)'},'Loading',cDataDir);
                    fileformat=fn(end-2:end);
                    filename= strcat(pn,fn);
                    switch fileformat
                        case 'mat'
                            load(filename);
                    end
                    this.dPos_mm = dPos_mm;
                    this.uieScanSteps.set(ceil(sqrt(size(dPos_mm,1))));
                    this.uieScanRange.set(max(max(dPos_mm(:,1:2),[],1)-min(dPos_mm(:,1:2),[],1)));
                    
                    % Make phase tab active:
                    this.uitgAxesDisplay.selectTabByIndex(this.U8ANALYSIS);
                    this.uipSelectObject.setSelectedIndex(uint8(1));
                    % Plot wavefronts on phase tab
                    this.replot(this.U8ANALYSIS, []);
                    
                case this.uieScanSteps
                    dN = this.uieScanSteps.get();
                    dL = this.uieScanRange.get();
                    dPosString = sprintf('%0.6f:%0.6f/%d:%0.6f', -dL/2,dL, dN-1, dL/2);
                    dPos = eval(dPosString);
                    [dx,dy] = meshgrid(dPos,dPos);
                    this.dPos_mm  = [dx(:),dy(:)];
                    this.uiePhaseStepsSim.set(sprintf('[%s;%s]''', dPosString, dPosString));
                    Lo_mm = this.uieLo.get();
                    dA =atan(eval(this.uiePhaseStepsSim.get())/Lo_mm)/pi*180;
                    this.uieScanAngles.set(mat2str(dA));
                    if strcmp(this.uipSelectObject.getOptions{this.uipSelectObject.getSelectedIndex()},...
                            'Scanning position')&&strcmp(this.uitgAxesDisplay.getSelectedTabName(),'Analysis')
                        % Plot wavefronts on phase tab
                        this.replot(this.U8ANALYSIS, []);
                    end
                    
                    % Probe and Object
                case this.uibSimulatePO
                    this.simulateInteferograms(true);
                case this.uieRprobe
                    probeType = this.uipProbeType.getOptions{this.uipProbeType.getSelectedIndex()};
                    if strcmp(probeType,'Defocus wave')
                        Rc_mm = this.uieRprobe.get();
                        z_mm = this.uiez2.get();
                        df_mm = this.uiez1.get();
                        NA = sin(atan(Rc_mm / (z_mm+df_mm)));
                        this.uieNA.set(NA);
                    end
                case this.uieProbeOffset
                    this.validateCouplesEditBox(src, '[]');
                    
                case this.uipProbeType
                    probeType = this.uipProbeType.getOptions{this.uipProbeType.getSelectedIndex()};
                    switch probeType
                        case 'Defocus wave'
                            this.cb(this.uieRprobe);
                            this.uipPropagator.setSelectedIndex(uint8(2));
                        case 'Plane wave'
                            this.uieNA.set(0);
                            this.uipPropagator.setSelectedIndex(uint8(1));
                    end
                    
                case this.uibLoadProbe
                    cDataDir = fullfile(this.cAppPath,  '..', '..', 'Data','probe');
                    [fn,pn]=uigetfile({'*.mat','Probe (*.mat)'},'Loading',cDataDir);
                    u8ModeId = this.uilSelectMode.getSelectedIndexes();
                    modeNumber = this.uieModeNumber.get();
                    lambda_um   = this.uieLambda.get()/1000;
                    this.dLambda_um(u8ModeId) = lambda_um;
                    fileformat=fn(end-2:end);
                    filename= strcat(pn,fn);
                    switch fileformat
                        case 'mat'
                            load(filename);
                    end
                    N           = this.uieRes.get();
                    [m,n] = meshgrid(linspace(0,1,N));
                    [sr,sc]= size(probe);
                    [p,q] = meshgrid(linspace(0,1,sc),linspace(0,1,sr));
                    
                    initialGuess = this.uicbGuess.get();
                    if initialGuess
                        this.dProbeGuess(:,:,u8ModeId) = interp2(p,q,probe,m,n,'nearest');
                        for mn = 2:modeNumber
                            this.dProbeGuess(:,:,mn) = this.dProbeGuess(:,:,1);
                        end
                    else
                        this.dProbe(:,:,u8ModeId) = interp2(p,q,probe,m,n,'nearest');
                        for mn = 2:modeNumber
                            this.dProbe(:,:,mn) = this.dProbe(:,:,1);
                        end
                    end
                    if initialGuess
                        % Make phase tab active:
                        this.uitgAxesDisplay.selectTabByIndex(this.U8GUESS);
                        % Plot wavefronts on phase tab
                        this.replot(this.U8GUESS, []);
                        % Set state:
                        this.setState(this.U8STATE_GUESS_PROCESSED);
                    else
                        % Make phase tab active:
                        this.uitgAxesDisplay.selectTabByIndex(this.U8PROBEOBJECT);
                        % Plot wavefronts on phase tab
                        this.replot(this.U8PROBEOBJECT, []);
                        % Set state:
                        this.setState(this.U8STATE_PHASE_PROCESSED);
                    end
                    
                case this.uibLoadObject
                    cDataDir = fullfile(this.cAppPath,  '..', '..', 'Data','object');
                    [fn,pn]=uigetfile({'*.mat','Object (*.mat)'},'Loading',cDataDir);
                    fileformat=fn(end-2:end);
                    filename= strcat(pn,fn);
                    switch fileformat
                        case 'mat'
                            load(filename);
                    end
                    N           = this.uieRes.get();
                    lambda_um   = this.uieLambda.get()/1000;
                    u8ModeId = this.uilSelectMode.getSelectedIndexes();
                    FP = this.uicbFourierPtychography.get();
                    z_um       = this.uiez2.get()*1000;
                    scanRange_um  = this.uieScanRange.get()*1000;
                    probeOffset = eval(this.uieProbeOffset.get());
                    if isempty(probeOffset)
                        probeOffset=[0,0];
                    else
                        probeOffset=round(probeOffset*1000./this.do_um(u8ModeId));
                    end
                    propagator = this.uipPropagator.getOptions{this.uipPropagator.getSelectedIndex()};
                    u8ModeId = this.uilSelectMode.getSelectedIndexes();
                    modeNumber = this.uieModeNumber.get();
                    this.dLambda_um(u8ModeId) = lambda_um;
                    samplingFactor_det = lambda_um.*z_um/(this.dc_um*N*this.dc_um);
                    if samplingFactor_det>1&&~strcmp(propagator,'fourier')
                        fprintf('Please adjust configurations for propagation sampling\n');
                    end
                    
                    if FP
                        k0 = 2*pi/lambda_um;
                        %                         kxy = k0*sin(atan(this.dPos_mm*1000/z_um));
                        %                         kr = k0*sin(sqrt(this.dPos_mm(:,1).^2+this.dPos_mm(:,2).^2)*1000/z_um);
                        kr = k0*sin(atan(sqrt(this.dPos_mm(:,1).^2+this.dPos_mm(:,2).^2)/this.uieLo.get()));
                        phi = atan2(this.dPos_mm(:,1),this.dPos_mm(:,2));
                        kxy = [kr.*sin(phi),kr.*cos(phi)];
                        dkxy = 2*pi/this.dc_um/N*this.uieMag.get();
                        dPosShifts = kxy./dkxy;
                        dPosShifts = dPosShifts -min(dPosShifts,[],1);
                        dPosShifts = round(dPosShifts);
                        K = max(dPosShifts(:,1))+N;
                        L = max(dPosShifts(:,2))+N;
                        if abs(K-L)<=2
                            K=max(K,L);
                            L=max(K,L);
                        end
                    else
                        dPosShifts = round((this.dPos_mm(:,1:2)-min(this.dPos_mm(:,1:2),[],1))*1000/this.do_um(u8ModeId));
                        K = max(dPosShifts(:,1))+N;
                        L = max(dPosShifts(:,2))+N;
                        if abs(K-L)<=2
                            K=max(K,L);
                            L=max(K,L);
                        end
                    end
                    if u8ModeId==1
                        this.dMaxObjectLen = K;
                    elseif K<this.dMaxObjectLen
                        K = this.dMaxObjectLen;
                        L = this.dMaxObjectLen;
                    elseif K>this.dMaxObjectLen
                        fprintf('Mode 1 wavelength needs to be the smallest!\n');
                        return;
                    end
                    if K>2000
                        fprintf('object sampling: %d, please adjust scanning range\n',K);
                        return;
                    end
                    [m,n] = meshgrid(linspace(0,1,L),linspace(0,1,K));
                    [sr,sc]= size(object);
                    [p,q] = meshgrid(linspace(0,1,sc),linspace(0,1,sr));
                    object = interp2(p,q,object,m,n,'nearest');
                    if FP
                        object =  PIE.utils.Propagate (object,propagator,this.do_um(u8ModeId),lambda_um,-1);
                        object = circshift(object,[probeOffset(2), probeOffset(1)]);
                    end
                    initialGuess = this.uicbGuess.get();
                    if initialGuess
                        this.dObjectGuess(:,:,u8ModeId) = object;
                        for mn = 2:modeNumber
                            this.dObjectGuess(:,:,mn) = this.dObjectGuess(:,:,1);
                        end
                    else
                        this.dObject(:,:,u8ModeId) = object;
                        for mn = 2:modeNumber
                            this.dObject(:,:,mn) = this.dObject(:,:,1);
                        end
                    end
                    if initialGuess
                        % Make phase tab active:
                        this.uitgAxesDisplay.selectTabByIndex(this.U8GUESS);
                        % Plot wavefronts on phase tab
                        this.replot(this.U8GUESS, []);
                        % Set state:
                        this.setState(this.U8STATE_GUESS_PROCESSED);
                    else
                        % Make phase tab active:
                        this.uitgAxesDisplay.selectTabByIndex(this.U8PROBEOBJECT);
                        % Plot wavefronts on phase tab
                        this.replot(this.U8PROBEOBJECT, []);
                        % Set state:
                        this.setState(this.U8STATE_PHASE_PROCESSED);
                    end
                    
                case this.uibGenProbeObject
                    this.generateProbeObject();
                    
                case this.uibCopyProbe
                    if this.uicbGuess.get()
                        this.dProbe = this.dProbeGuess;
                    else
                        this.dProbeGuess = this.dProbe;
                    end
                    % FPM
                case this.uicbFourierPtychography
                    if this.uicbFourierPtychography.get()
                        this.uipPropagator.setSelectedIndex(uint8(2));
                    end
                    
                case this.uieNAo
                    if this.uicbFourierPtychography.get()
                        Lo_mm = this.uieLo.get();
                        Li_mm = this.uiez2.get();
                        NAo = this.uieNAo.get();
                        NAi = sin(atan(Lo_mm/Li_mm*tan(asin(NAo))));
                        this.uieNA.set(NAi);
                        lambda_mm =this.uieLambda.get()/1000000;
                        Rc_mm = 0.61*lambda_mm/NAi;
                        this.uieRprobe.set(Rc_mm);
                        this.uieMag.set(NAo/NAi);
                    end
                    
                case this.uieLo
                    if this.uicbFourierPtychography.get()
                        Lo_mm = this.uieLo.get();
                        Li_mm = this.uiez2.get();
                        NAo = this.uieNAo.get();
                        NAi = sin(atan(Lo_mm/Li_mm*tan(asin(NAo))));
                        this.uieNA.set(NAi);
                        lambda_mm =this.uieLambda.get()/1000000;
                        Rc_mm = 0.61*lambda_mm/NAi;
                        this.uieRprobe.set(Rc_mm);
                        this.uieMag.set(NAo/NAi);
                    end
                case this.uieMag
                    if this.uicbFourierPtychography.get()
                        Li_mm = this.uiez2.get();
                        NAo = this.uieNAo.get();
                        Mag = this.uieMag.get();
                        Lo_mm = Li_mm/Mag;
                        if Lo_mm - this.uieLo.get()>1e-4
                            this.uieLo.set(Lo_mm);
                            NAi = sin(atan(Lo_mm/Li_mm*tan(asin(NAo))));
                            this.uieNA.set(NAi);
                            lambda_mm =this.uieLambda.get()/1000000;
                            Rc_mm = 0.61*lambda_mm/NAi;
                            this.uieRprobe.set(Rc_mm);
                        end
                    end
                    
                case this.uieScanAngles
                    [lValid, vals] = this.validateCouplesEditBox(src, '[]');
                    if lValid
                        [sr, sc] = size(vals);
                        this.uieScanSteps.set(sr);
                        dL = this.uieScanRange.get();
                        Lo_mm = this.uieLo.get();
                        dP =tan(eval(this.uieScanAngles.get())/180*pi)*Lo_mm+dL/2;
                        this.uiePhaseStepsSim.set(mat2str(dP));
                    end
                    
                    % Sim Segments
                case this.uibLoadSegment
                    try
                        cDataDir = fullfile(this.cAppPath, '..', '+utils', '*.mat');
                        [d, p] = uigetfile(cDataDir);
                        this.uieSegmentPath.set([p d]);
                        load([p d]);
                        this.ceSegments = segs;
                        this.uieRes.set(length(segs{1}));
                    catch
                        this.ceSegments = [];
                    end
                    
                case this.uieSegmentPath
                    try
                        load(this.uieSegmentPath.get());
                        this.ceSegments = segs;
                    catch
                        this.ceSegments = [];
                    end
                    
                    % Sim stochastics
                case this.uibCustomSim
                    %                     this.customSimulation();
                    this.simStackAndReconstruct();
                    
                case this.uibReset
                    this.setSimParams('zDrift',0,'xDrift',0,'yDrift',0,'gratTilt',0,'detTilt',0,'shotToShot',0,...
                        '2ndOrderStrength',0,'11OrderStrength',0,'xyStageError',0,'detectorCurve',0,'flareLevel',0,...
                        'MSFR',0,'nonlinearity', 0,'airflow',0,'nPhotons',0);
                    
                case {this.uieZrn}
                    this.validateCouplesEditBox(src, '[]');
                    
                    
                    
                    %% Reconstruction
                case this.uibComputePhase
                    this.reconstruct(this.uitgAnalysisDomain.getSelectedTabName());
                    
                case this.uibStop
                    global stopSign
                    stopSign=1;
                    
                    % rPIE
                case this.uieModeNumber
                    modeNumber = this.uieModeNumber.get();
                    ceModeStr = cell(modeNumber,1);
                    for i= 1:modeNumber
                        ceModeStr{i} =['Mode ',num2str(i)];
                    end
                    try
                        this.uilSelectMode.setOptions(ceModeStr);
                    catch
                    end
                    
                case this.uibShowMode
                    u8ModeId = this.uilSelectMode.getSelectedIndexes();
                    this.uieLambda.set(this.dLambda_um(u8ModeId)*1e3);
                    switch this.uitgAxesDisplay.getSelectedTabName()
                        case 'Reconstruction'
                            % Plot wavefronts on phase tab
                            this.replot(this.U8RECONSTRUCTION, []);
                        case 'Probe and object'
                            this.replot(this.U8PROBEOBJECT, []);
                        case 'Initial guess'
                            this.replot(this.U8GUESS, []);
                        case 'Analysis'
                            this.replot(this.U8ANALYSIS, []);
                    end
                    
                    %% Analysis
                case this.uibAnalyze
                    this.analyze();
                    
                case this.uibLoadAnalysisRegion
                    cDataDir = fullfile(this.cAppPath,  '..', '..', 'Data','mask');
                    [fn,pn]=uigetfile({'*.mat','Mask (*.mat)'},'Loading',cDataDir);
                    fileformat=fn(end-2:end);
                    filename= strcat(pn,fn);
                    switch fileformat
                        case 'mat'
                            load(filename);
                    end
                    this.dAnalysisMask = mask;
                    imagesc(this.haAnalysis, mask);axis(this.haAnalysis,'xy','tight','equal');
                    this.haAnalysis.Title.String = 'Analysis region';
                    this.uipSelectRegion.setSelectedIndex(uint8(4));
                    
                case this.uibDataCursor
                    if this.uibDataCursor.get()
                        datacursormode(this.hFigure,'on');
                    else
                        datacursormode(this.hFigure,'off');
                    end
                    
                    %% Other
                case this.htLog
                    s=evt.Indices;
                    [~,N]=size(this.ceAnalysisTable);
                    str='[';
                    for i=1:N-4
                        str=strcat(str,num2str(i),',',num2str(this.ceAnalysisTable{s(1),i+4}),';');
                    end
                    str=[str(1:end-1),']'];
                    this.uieZrn.set(str);
                    
                case {this.uibLeft, this.uibRight}
                    if (src == this.uibLeft)
                        this.cycleInterferogram(-1);
                    else
                        this.cycleInterferogram(1);
                    end
                    
            end
        end
        
        function generateProbeObject(this)
            % load parameters
            probeType = this.uipProbeType.getOptions{this.uipProbeType.getSelectedIndex()};
            objectType = this.uipObjectType.getOptions{this.uipObjectType.getSelectedIndex()};
            probeOffset = eval(this.uieProbeOffset.get());
            probeAmp = this.uieProbeAmp.get();
            u8ModeId = this.uilSelectMode.getSelectedIndexes();
            modeNumber = this.uieModeNumber.get();
            initialGuess = this.uicbGuess.get();
            FP = this.uicbFourierPtychography.get();
            N           = this.uieRes.get();
            NA           = this.uieNA.get();
            Rc_um   = this.uieRprobe.get()*1000;
            lambda_um   = this.uieLambda.get()/1000;
            this.dLambda_um(u8ModeId) = lambda_um;
            df_um       = this.uiez1.get()*1000;% negative sign corresponds convergent
            z_um       = this.uiez2.get()*1000;
            scanRange_um  = this.uieScanRange.get()*1000;
            zernCouples = eval(this.uieZrn.get());
            scanSteps  = this.uieScanSteps.get();
            detSize_um  = this.uieDetSize.get()*1000;
            propagator = this.uipPropagator.getOptions{this.uipPropagator.getSelectedIndex()};
            CenterObstruction = this.uieCenterObstruction.get();
            this.dc_um= detSize_um/N;
            if strcmp(propagator,'fourier')
                this.do_um(u8ModeId) = lambda_um*z_um/N/this.dc_um;
            else
                this.do_um(u8ModeId) = this.dc_um; % object pixel pitch
            end
            samplingFactor_det = lambda_um.*z_um/(this.dc_um*N*this.dc_um);
            if samplingFactor_det>1&&~strcmp(propagator,'fourier')
                fprintf('Please adjust configurations for propagation sampling\n');
            end
            if isempty(probeOffset)
                probeOffset=[0,0];
            else
                probeOffset=round(probeOffset*1000./this.do_um(u8ModeId));
            end
            
            % generate probe phase based on zernike polynomials
            zfn         =   PIE.utils.generateZernikeFunction(zernCouples,N,1); % generate zernike function form based on zernike couples
            if ~FP
                [x_um,y_um] = meshgrid(linspace(-detSize_um/2,detSize_um/2,N));
                if NA>0
                    r = sin(atan(sqrt(x_um.^2+y_um.^2)/(z_um+df_um)))/NA;
                    th = atan2(y_um,x_um);
                else
                    r = sqrt(x_um.^2+y_um.^2)/Rc_um;
                    th = atan2(y_um,x_um);
                end
                probe_pha = 2*pi*zfn(r,th);
            end
            
            %% initial probe
            if ~FP
                switch probeType
                    case 'Defocus wave'
                        samplingFactor_obj = lambda_um.*z_um/(N*this.dc_um*this.do_um(u8ModeId));
                        Rc_um = (z_um+df_um)*tan(asin(NA));
                        Rprobe_um = abs(df_um)*tan(asin(NA));
                        [n1,n2]=meshgrid(1:N);
                        n1 = n1-N/2-1;
                        n2 = n2-N/2-1;
                        Hs= exp(-1i*pi*df_um*this.dc_um^2/lambda_um/z_um^2*(n1.^2+n2.^2));
                        CO = 1-pinhole(round(2*Rc_um/this.dc_um*CenterObstruction),N,N);
                        probe = PIE.utils.Propagate (pinhole(round(2*Rc_um/this.dc_um),N,N).*CO.*Hs.*exp(1i*probe_pha),propagator,...
                            this.do_um(u8ModeId),lambda_um,-df_um);
                    case 'Plane wave'
                        Rprobe_um = Rc_um;
                        samplingFactor_obj = lambda_um.*abs(df_um)/(N*this.do_um(u8ModeId)*this.do_um(u8ModeId));
                        CO = 1-pinhole(round(2*Rprobe_um/this.dc_um*CenterObstruction),N,N);
                        probe = PIE.utils.Propagate (pinhole(round(2*Rprobe_um/this.do_um(u8ModeId)),N,N).*CO.*exp(1i*probe_pha),propagator,...
                            this.do_um(u8ModeId),lambda_um,abs(df_um));
                    case 'From reconstruction'
                        probe = this.dProbeRecon;
                end
                probe = probeAmp*circshift(probe,[probeOffset(2), probeOffset(1)]);
            else % define probe in pupil space for fourier ptychography
                samplingFactor_obj = lambda_um.*z_um/(N*this.dc_um*this.do_um(u8ModeId));
                Rprobe_um = z_um*tan(asin(NA));
                k0=2*pi/lambda_um;
                cutoff = NA*k0;
                kmax=pi/this.dc_um;
                [kxm,kym] =meshgrid(linspace(-kmax,kmax,N));
                kr = sqrt(kxm.^2+kym.^2);
                phi = atan2(kym,kxm);
                this.dCTF = (kr<cutoff);
                this.dCTF = circshift(this.dCTF,[probeOffset(2), probeOffset(1)]);
                NAo = this.uieNAo.get();
                kzm = sqrt(k0^2-(kxm/NA*NAo).^2-(kym/NA*NAo).^2);
                defocus_pha = exp(1i.*df_um.*real(kzm)).*exp(abs(df_um).*abs(imag(kzm)));
                probe_pha = 2*pi*zfn(kr/cutoff,phi);
                probe = probeAmp*this.dCTF.*exp(1i*probe_pha).*defocus_pha;
                
                %                 [xp_um,yp_um] = meshgrid(linspace(-do_um*N/2,do_um*N/2,N));
                %                 Lo_um = this.uieLo.get()*1000;
                %                 NAo = this.uieNAo.get();
                %                 probe_amp = zeros(N);
                %                 probe_amp((xp_um.^2+yp_um.^2<= Rprobe_um.^2)&(xp_um.^2+yp_um.^2> (Rprobe_um*CenterObstruction).^2)) = 1;
                %                 kr = sin(atan(sqrt(xp_um.^2+yp_um.^2)/Lo_um));
                %                 phi = atan2(yp_um,xp_um);
                %                 kx = kr.*cos(phi);
                %                 ky = kr.*sin(phi);
                %                 kz = 2*pi/lambda_um*sqrt(1-kx.^2-ky.^2);
                %                 defocus_pha = df_um.*real(kz);
            end
            if initialGuess
                if u8ModeId==1
                    this.dProbeGuess = single(probe);
                    for m = 2:modeNumber
                        this.dProbeGuess (:,:,m) = single(probe);
                    end
                else
                    this.dProbeGuess(:,:,u8ModeId) = single(probe);
                end
            else
                if u8ModeId==1
                    this.dProbe = single(probe);
                    for m = 2:modeNumber
                        this.dProbe(:,:,m) = single(probe);
                    end
                else
                    this.dProbe(:,:,u8ModeId) = single(probe);
                end
            end
            
            if samplingFactor_obj>1
                fprintf('Please adjust configurations for propagation sampling\n');
            end
            airyR_um = 0.61*lambda_um/this.uieNA.get();
            
            if ~isempty(this.ceSegments)
                NSeg = length(this.ceSegments);
                segFactor = sqrt(NSeg)/N;
            else
                segFactor = 1;
            end
            if size(this.dPos_mm,2)==3
                Nz = length(unique(this.dPos_mm(:,3)));
            else
                Nz=1;
            end
            dR = scanRange_um/(sqrt(scanSteps^2/Nz)-1); % scanning interval
            if Rprobe_um> airyR_um
                overlap = PIE.utils.overlapRatio(Rprobe_um,dR); % overlap ratio of two circles
                %                 samplingFactor = 2*Rprobe_um/dR*segFactor*sqrt(Nz);
            else
                overlap = PIE.utils.overlapRatio(airyR_um,dR); % overlap ratio of two circles
                %                 samplingFactor = 2*airyR_um/dR*segFactor*sqrt(Nz);
            end
            samplingFactor = this.do_um(u8ModeId)/dR*segFactor*N/2*sqrt(Nz);
            this.uitOverlap.set(['Overlap: ',num2str(round(overlap*100)),'%']); % check overlap
            
            this.uitSampling.set(['Sampling factor: ',num2str(round(samplingFactor*10)/10)]); % check sampling
            
            %% initial object
            if FP
                k0 = 2*pi/lambda_um;
                %                 kxy = k0*sin(atan(this.dPos_mm*1000/z_um));
                %                 kr = k0*sin(sqrt(this.dPos_mm(:,1).^2+this.dPos_mm(:,2).^2)*1000/z_um);
                kr = k0*sin(atan(sqrt(this.dPos_mm(:,1).^2+this.dPos_mm(:,2).^2)/this.uieLo.get()));
                phi = atan2(this.dPos_mm(:,1),this.dPos_mm(:,2));
                kxy = [kr.*sin(phi),kr.*cos(phi)];
                dkxy = 2*pi/this.dc_um/N*this.uieMag.get();
                dPosShifts = kxy./dkxy;
                dPosShifts = dPosShifts -min(dPosShifts,[],1);
                dPosShifts = round(dPosShifts);
                K = max(dPosShifts(:,1))+N;
                L = max(dPosShifts(:,2))+N;
                if abs(K-L)<=2
                    K=max(K,L);
                    L=max(K,L);
                end
            else
                dPosShifts = round((this.dPos_mm(:,1:2)-min(this.dPos_mm(:,1:2),[],1))*1000/this.do_um(u8ModeId));
                K = max(dPosShifts(:,1))+N;
                L = max(dPosShifts(:,2))+N;
                if abs(K-L)<=2
                    K=max(K,L);
                    L=max(K,L);
                end
            end
            if u8ModeId==1
                this.dMaxObjectLen = K;
            elseif K<this.dMaxObjectLen
                K = this.dMaxObjectLen;
                L = this.dMaxObjectLen;
            elseif K>this.dMaxObjectLen
                fprintf('Mode 1 wavelength needs to be the smallest!\n');
                return;
            end
            %             K = round(scanRange_um/this.do_um)+N;
            %             L = round(scanRange_um/this.do_um)+N; % size of object [K,L]
            if round(scanRange_um/this.do_um(u8ModeId))<scanSteps/sqrt(Nz)-1
                this.uitSampling.set(['Over scanning ',num2str(round(scanRange_um/this.do_um(u8ModeId)))]); % check sampling
            end
            if K>2000
                fprintf('object sampling: %d, please adjust scanning range\n',K);
                return;
            end
            switch objectType
                case 'Vacuum'
                    object = single(ones(K,L));
                case 'Cameraman'
                    I = single(flipud(imread('cameraman.tif')));
                    [m,n] = meshgrid(linspace(0,1,L),linspace(0,1,K));
                    [sr,sc,~]= size(I);
                    [p,q] = meshgrid(linspace(0,1,sc),linspace(0,1,sr));
                    object_amp = interp2(p,q,I(:,:,1),m,n,'nearest');
                    object_amp = mat2gray(object_amp)*0.8+0.2; % object amplitude
                    %                     object_amp = ones(K,L);
                    I =single(imread('pears.png'));
                    [sr,sc,~]= size(I);
                    [p,q] = meshgrid(linspace(0,1,sc),linspace(0,1,sr));
                    object_pha = interp2(p,q,I(:,:,1),m,n,'nearest');
                    object_pha=mat2gray(object_pha);
                    object_pha = (object_pha-0.5)*1*pi; % object phase
                    object = object_amp.*exp(1i*object_pha);
                case 'Low resolution'
                    exitWave = this.dObject.*pad2(probe,K,L);
                    object = PIE.utils.postPropagate (exitWave,propagator,1,1);
                case 'Load object'
                    if initialGuess
                        object = this.dObjectGuess(:,:,u8ModeId);
                    else
                        object = this.dObject(:,:,u8ModeId);
                    end
                case 'From reconstruction'
                    object = this.dObjectRecon(:,:,u8ModeId);
                    if FP
                        object =  PIE.utils.Propagate (object,propagator,this.do_um(u8ModeId),lambda_um,1);
                    end
            end
            if FP
                object =  PIE.utils.Propagate (object,propagator,this.do_um(u8ModeId),lambda_um,-1);
                object = circshift(object,[probeOffset(2), probeOffset(1)]);
            end
            
            if initialGuess
                if u8ModeId==1
                    this.dObjectGuess =  object;
                    for m = 2:modeNumber
                        this.dObjectGuess(:,:,m) = object;
                    end
                else
                    this.dObjectGuess(:,:,u8ModeId) = object;
                end
            else
                if u8ModeId==1
                    this.dObject = object;
                    for m = 2:modeNumber
                        this.dObject(:,:,m) = object;
                    end
                else
                    this.dObject(:,:,u8ModeId) = object;
                end
            end
            
            if initialGuess
                % Make phase tab active:
                this.uitgAxesDisplay.selectTabByIndex(this.U8GUESS);
                % Plot wavefronts on phase tab
                this.replot(this.U8GUESS, []);
                % Set state:
                this.setState(this.U8STATE_GUESS_PROCESSED);
            else
                % Make phase tab active:
                this.uitgAxesDisplay.selectTabByIndex(this.U8PROBEOBJECT);
                % Plot wavefronts on phase tab
                this.replot(this.U8PROBEOBJECT, []);
                % Set state:
                this.setState(this.U8STATE_PHASE_PROCESSED);
            end
        end
        
        
        function setZrnString(this, zrn)
            % build string:
            zrn = zrn(:);
            zrns = (0:length(zrn) - 1);
            zrns(2,:) = zrn';
            this.uieZrn.set(mat2str(zrns'));
        end
        
        
        function [fid, isCreated] = openOrCreateFile(this, fullFilePath)
            [d p e] = fileparts(fullFilePath);
            
            % Check if dir exists:
            saFls = dir(d);
            if isempty(saFls)
                % make the dir:
                mkdir(d);
            end
            
            % now check if a file exists:
            fid = fopen(fullFilePath, 'r');
            
            if (fid == -1)
                isCreated = true;
            else
                fclose(fid);
                isCreated = false;
            end
            fid = fopen(fullFilePath, 'a');
            
        end
        
        % validates whether a char edit box evaluates to a Nx2 matrix,
        % colors accordingly.  Empty value is changed to []
        function [lOut, vals] = validateCouplesEditBox(~, src, cDefaultVal)
            lOut = true;
            vals = [];
            if isempty(src.get())
                src.styleDefault();
                src.set(cDefaultVal);
                return
            end
            try
                vals = eval(src.get());
                [~, sc] = size(vals);
                if (sc == 2 || sc == 0)
                    src.styleDefault();
                    lOut = true;
                else
                    src.styleBad();
                    lOut = false;
                end
            catch
                % can't read this edit box
                src.styleBad();
                lOut = false;
            end
        end
        
        
        function setState(this, dStateFlag)
            this.u8State = dStateFlag;
            
            % all uib elements below this state flag go to default color:
            for k = 1:dStateFlag - 1
                for m = 1:length(this.ceUIBStateControlled{k})
                    this.ceUIBStateControlled{k}{m}.setColor(this.dColorDefault);
                end
            end
            
            % all uib elements at this state flag to active color:
            for m = 1:length(this.ceUIBStateControlled{dStateFlag})
                this.ceUIBStateControlled{dStateFlag}{m}.setColor(this.dColorActive);
            end
            
            % all uib elements above this state flag go to forbid color:
            for k = dStateFlag + 1:this.U8MAXSTATES
                for m = 1:length(this.ceUIBStateControlled{k})
                    this.ceUIBStateControlled{k}{m}.setColor(this.dColorForbid);
                end
            end
            
        end
        
        
        % reads log CSV and populates the uil lists:
        function parseLogFile(this, cPath, cDir)
            this.ceScanPaths = {};
            this.ceSeriesPaths = {{}};
            
            try
                this.tSelectedLog   = readtable(cPath);
            catch
                msgbox(sprintf('Unable to read file %s', cPath));
                return
            end
            
            ceTable             = table2cell(this.tSelectedLog);
            
            seriesColIdx    = (strcmp(this.tSelectedLog.Properties.VariableNames, 'seriesIndex'));
            scanColIdx      = (strcmp(this.tSelectedLog.Properties.VariableNames, 'scanIndex'));
            dateColIdx      = (strcmp(this.tSelectedLog.Properties.VariableNames, 'timeStamp'));
            scanAxesIdx     = (strcmp(this.tSelectedLog.Properties.VariableNames, 'scanAxes'));
            scanOutputIdx   = (strcmp(this.tSelectedLog.Properties.VariableNames, 'scanOutput'));
            fileNameIdx     = (strcmp(this.tSelectedLog.Properties.VariableNames, 'fileName'));
            
            % Filter ceTable for scanOutput == 'Image capture'
            isValidScanType = @(type) any(cellfun(@(x) strcmp(x, type), this.ceValidForAnalysis));
            ceTable(cellfun(@(x) ~isValidScanType(x), ceTable(:,scanOutputIdx)), :) = [];
            
            % find unique series numbers
            dSeriesNumbers = unique(cell2mat(ceTable(:, seriesColIdx)));
            
            dImgsInSeries = [];
            for k = 1:length(dSeriesNumbers)
                dImgsInSeries(k) = size(ceTable(cellfun(@(x) x == dSeriesNumbers(k), ceTable(:,seriesColIdx)), :), 1); %#ok<AGROW>
            end
            
            
            fhBuildSeriesLabel = @(cDateStr, seriesIdx, dNum, cScanIndex) sprintf('%s - Series %d, %d images: %s', cDateStr, seriesIdx, dNum, cScanIndex);
            fhBuildScanLabel = @(cDateStr, seriesIdx, scanIdx, dNum) sprintf('%s - Series %d: %d/%d', cDateStr, seriesIdx, scanIdx, dNum);
            fhBUildGenericLabel = @(cDateStr, cFileName) sprintf('%s - %s', cDateStr, cFileName);
            % Generate scan metadata:
            
            ceSeriesDates = {};
            ceSeriesScanAxisType = {};
            
            ceListOptionsScan = {};
            ceListOptionsSeries = {};
            
            for k = 1:size(ceTable,1)
                % If not a series, just add it
                if ~any(seriesColIdx) || isempty(ceTable{k, seriesColIdx})
                    ceListOptionsScan{k} = ...
                        fhBUildGenericLabel(ceTable{k, dateColIdx}, ceTable{k, fileNameIdx});
                    [~, d, ~] = fileparts(ceTable{k, fileNameIdx});
                    this.ceScanPaths{k} =  fullfile(cDir, [d '.mat']);
                else
                    dSeriesIdx = find(dSeriesNumbers == ceTable{k, seriesColIdx});
                    if length(ceSeriesDates) < dSeriesIdx
                        ceSeriesDates{end + 1} = ceTable{k, dateColIdx}; %#ok<*AGROW>
                        ceSeriesScanAxisType{end + 1} = ceTable{k, scanAxesIdx};
                    end
                    this.ceScanPaths{k} = sprintf('%sseries_%0.3d/%s.mat', cDir, ceTable{k, seriesColIdx}, ceTable{k, fileNameIdx});
                    
                    if length(this.ceSeriesPaths) < dSeriesIdx
                        this.ceSeriesPaths{dSeriesIdx} = cell(0);
                    end
                    dSeriesNum = length(this.ceSeriesPaths{dSeriesIdx});
                    this.ceSeriesPaths{dSeriesIdx}{dSeriesNum + 1} = sprintf('%sseries_%0.3d/%s.mat', cDir, ceTable{k, seriesColIdx}, ceTable{k, fileNameIdx});
                    dImgsInThisSeries = dImgsInSeries(dSeriesNumbers == ceTable{k, seriesColIdx});
                    ceListOptionsScan{k} = fhBuildScanLabel(...
                        ceTable{k, dateColIdx}, ...
                        ceTable{k, seriesColIdx}, ...
                        ceTable{k, scanColIdx}, ...
                        dImgsInThisSeries);
                end
            end
            
            this.uilSingleList.setOptions(ceListOptionsScan);
            
            for k = 1:length(dSeriesNumbers)
                ceListOptionsSeries{k} = fhBuildSeriesLabel(...
                    ceSeriesDates{k}, dSeriesNumbers(k), ...
                    dImgsInSeries(k), ceSeriesScanAxisType{k});
            end
            
            this.uilStackList.setOptions(ceListOptionsSeries);
            
            
            fprintf('Parsed log file %s\n', cPath);
        end
        
        function reconstruct(this, cAnalysisDomain)
            if isempty(this.ceInt)
                msgbox('Please load/simulate diffraction patterns first!', 'Error');
                return;
            end
            global stopSign
            stopSign = 0; % used for stopping iteration
            
            % load parameter
            scanSteps               = this.uieScanSteps.get();
            lambda_um = this.uieLambda.get()/1000;
            FP = this.uicbFourierPtychography.get();
            z_um = this.uiez2.get()*1000;
            modeNumber = this.uieModeNumber.get();
            N = length(this.ceInt{1,1});
            propagator = this.uipPropagator.getOptions{this.uipPropagator.getSelectedIndex()};
            
            % precompute FFT
            for j = 1:modeNumber
                H{j} = PIE.utils.prePropagate (this.dProbeGuess(:,:,j),propagator,this.do_um(j),this.dLambda_um(j),z_um,1);
                Hm{j} = PIE.utils.prePropagate (this.dProbeGuess(:,:,j),propagator,this.do_um(j),this.dLambda_um(j),-z_um,1);
                
                % generate scanning coordinates
                if FP
                    k0 = 2*pi/this.dLambda_um(j);
                    %                 kxy = k0*sin(atan(this.dPos_mm*1000/z_um));
                    kr = k0*sin(atan(sqrt(this.dPos_mm(:,1).^2+this.dPos_mm(:,2).^2)/this.uieLo.get()));
                    phi = atan2(this.dPos_mm(:,1),this.dPos_mm(:,2));
                    kxy = [kr.*sin(phi),kr.*cos(phi)];
                    dkxy = 2*pi/this.dc_um/N*this.uieMag.get();
                    Rpix = kxy./dkxy;
                    Rpix = Rpix -min(Rpix,[],1);
                    Rpix = round(Rpix);
                    %                 [K,L] = size(this.dObjectGuess);
                else
                    Rpix = round((this.dPos_mm(:,1:2)-min(this.dPos_mm(:,1:2),[],1))*1000/this.do_um(j));
                end
                if size(this.dPos_mm,2) ==3
                    Rpix(:,3) = this.dPos_mm(:,3)*1000;
                    zMax = this.do_um(j)^2*N/lambda_um;
                    DoF = lambda_um/(this.uieNA.get())^2;
                    fprintf('DoF is %0.3fum, and maximum propagation distance for 3-D scanning is %0.3fum, while the input is %0.3fum \n',DoF,zMax,abs(max(Rpix(:,3))));
                end
                ceRpix{j} = Rpix;
            end
            
            nInt = size(this.dPos_mm,1);
            for j =1:nInt
                for id = 1:modeNumber
                    mRpix(id,:,j) = ceRpix{id}(j,:);
                end
            end
            %             dPosShifts = eval(this.uiePhaseStepsSim.get());
            %             Rm = round(dPosShifts(:,1)*1000/do_um);
            %             Rn = round(dPosShifts(:,2)*1000/do_um);
            
            % generate scanning position
            
            if isempty(this.dProbeGuess)
                this.dProbeGuess = pinhole(round(N/2),N,N);
            end
            if isempty(this.dObjectGuess)
                KL = round(this.uieScanRange.get()*1000./this.do_um(1))+N;
                this.dObjectGuess = ones(KL);
            else
                KL = length(this.dObjectGuess);
            end
            this.dProbeRecon =this.dProbeGuess;
            this.dObjectRecon =this.dObjectGuess;
            sqrtInt  = single(zeros(N,N,scanSteps^2));
            Is = 0;
            ceInts = this.ceInt(:);
            for j =1:nInt
                sqrtInt(:,:,j)=sqrt(ceInts{j});
                Is = Is + ceInts{j};
            end
            
            totalI =0;
            Iseg = zeros(length(this.ceSegments),scanSteps^2);
            if ~isempty(this.ceSegments)
                for k =1:length(this.ceSegments)
                    for j =1:nInt
                        Iseg(k,j) = mean ( ceInts{j}(this.ceSegments{k}==1));
                    end
                    totalI = totalI+mean(Is(this.ceSegments{k}==1));
                end
            else
                totalI =totalI+sum(Is(:));
            end
            
            % initial reconstruction parameters
            iteration = this.uieMaxIteration.get();
            errors = zeros(iteration,1);
            alpha = this.uieAlpha.get(); % weight factor for updating object
            beta = this.uieBeta.get(); % weight factor for updating probe
            gamma = this.uieGamma.get(); % weight factor for rPIE, 1 for ePIE
            delta = this.uieDelta.get(); % weight factor for RAAR, 1 for DM
            epsilon = this.uieEpsilon.get(); % weight factor for WDD
            correctPos = this.uicbCorrectPos.get();
            GPUacceleration = this.uicbGPU.get();
            correctMethod =  this.uipCorrectMethod.getOptions{this.uipCorrectMethod.getSelectedIndex()};
            if GPUacceleration == 0
                % iterations for reconstruction
                for i = 1:iteration
                    tempError = 0;
                    switch cAnalysisDomain % reconstruction method
                        case 'rPIE' % scanning solution
                            for j =1:nInt
                                if correctPos==1 % apply position correction
                                    if i==1 % inital parameters
                                        Cpix = zeros(scanSteps^2,2);
                                        rot_rad =0;
                                        scale = 0;
                                        nC = 10;
                                        centralPix =max(Rpix(:,1))-min(Rpix(:,1));
                                        rCpix0 = floor(centralPix/(scanSteps-1)); % random position searching radius
                                        maxRot_rad0 = 1/180*pi; % maximum rotation
                                        maxScale0 =0.01;% maximum scale
                                    end
                                    decreaseFactor = (iteration-i)/iteration;
                                    [this.dObjectRecon,this.dProbeRecon,tempError,Cpix(j,:),rot_rad,scale] = ...
                                        PIE.utils.pcPIE(this.dObjectRecon,this.dProbeRecon,sqrtInt(:,:,j),...
                                        tempError,alpha,beta,gamma,propagator,H{1},Hm{1},N,KL,Rpix(j,:),Cpix(j,:),...
                                        rot_rad,scale, nC,rCpix0,maxRot_rad0,maxScale0,decreaseFactor,centralPix,correctMethod);
                                else % normal rPIE without position calibration
                                    if ~isempty( this.ceSegments)
                                        [this.dObjectRecon,this.dProbeRecon,tempError] = PIE.utils.sPIE(this.dObjectRecon,this.dProbeRecon,...
                                            sqrtInt(:,:,j),Iseg(:,j),mRpix(:,:,j),N,propagator,H,Hm,alpha,beta, gamma, tempError,modeNumber,this.ceSegments,this.do_um,this.dLambda_um);
                                    else
                                        [this.dObjectRecon,this.dProbeRecon,tempError] = PIE.utils.rPIE(this.dObjectRecon,this.dProbeRecon,...
                                            sqrtInt(:,:,j),mRpix(:,:,j),N,propagator,H,Hm,alpha,beta, gamma, tempError, modeNumber,this.do_um,this.dLambda_um,FP,this.dCTF);
                                    end
                                end
                                drawnow;
                                if strcmp(this.uipSelectObject.getOptions{this.uipSelectObject.getSelectedIndex()},...
                                        'Scanning position')&&strcmp(this.uitgAxesDisplay.getSelectedTabName(),'Analysis')
                                    % Plot wavefronts on phase tab
                                    this.replot(this.U8ANALYSIS, j);
                                end
                            end
                            
                        case 'RAAR' % batch scanning solution
                            if i==1 % initial exitWaves
                                exitWaves = zeros(N,N,scanSteps^2);
                                for j =1:scanSteps^2
                                    reconBox = this.dObjectRecon(Rpix(j,1)+[1:N],Rpix(j,2)+[1:N]);
                                    exitWaves(:,:,j) = this.dProbeRecon.*reconBox;
                                end
                            end
                            [this.dObjectRecon,this.dProbeRecon,tempError] = PIE.utils.RAAR(this.dObjectRecon,this.dProbeRecon,...
                                sqrtInt,exitWaves,tempError,alpha,beta,delta,Rpix,N,propagator,H{1},Hm{1},scanSteps);
                            
                        case 'WDD' % Wigner distribution deconvolution
                            [this.dObjectRecon,this.dProbeRecon] = PIE.utils.WDD(this.dObjectGuess,this.dProbeGuess,sqrtInt,N,Rpix,scanSteps,epsilon);
                            stopSign=1;
                            
                        case 'ML' % maximum-likelihood method
                            likelihoodType = this.uipLikelihoodType.getOptions{this.uipLikelihoodType.getSelectedIndex()};
                            regularization = this.uieRegularization.get(); % regularization factor for ML
                            aj = 1/2/N^2*ones(scanSteps^2,1);
                            for j =1:scanSteps^2
                                if ~isempty( this.ceSegments)
                                    [this.dObjectRecon,this.dProbeRecon,tempError,aj(j),alpha(j),beta(j)] = PIE.utils.sML(this.dObjectRecon,this.dProbeRecon,...
                                        sqrtInt(:,:,j),Iseg(:,j),Rpix(j,:),N,propagator,H{1},Hm{1}, gamma, tempError, likelihoodType,regularization, aj(j),this.ceSegments);
                                else
                                    [this.dObjectRecon,this.dProbeRecon,tempError,aj(j),alpha(j),beta(j)] = PIE.utils.ML(this.dObjectRecon,this.dProbeRecon,...
                                        sqrtInt(:,:,j),Rpix(j,:),N,propagator,H{1},Hm{1}, gamma, tempError, likelihoodType,regularization, aj(j));
                                end
                                drawnow;
                                if strcmp(this.uipSelectObject.getOptions{this.uipSelectObject.getSelectedIndex()},...
                                        'Scanning position')&&strcmp(this.uitgAxesDisplay.getSelectedTabName(),'Analysis')
                                    % Plot wavefronts on phase tab
                                    this.replot(this.U8ANALYSIS, j);
                                end
                            end
                            this.uieAlpha.set(mean(abs(alpha)));
                            this.uieBeta.set(mean(abs(beta)));
                    end
                    % error evaluation
                    errors(i) = sum(tempError(:))/totalI;
                    
                    iterationStr = sprintf('%d iterations finished,residual error: %0.5f',i,errors(i));
                    this.uitIteration.set(iterationStr);drawnow;
                    
                    if stopSign==1||(i>1&&errors(i)<this.uieAccuracy.get())||i==iteration
                        this.dError = errors;
                        % Make phase tab active:
                        this.uitgAxesDisplay.selectTabByIndex(this.U8RECONSTRUCTION);
                        % Plot wavefronts on phase tab
                        this.replot(this.U8RECONSTRUCTION, []);
                        break;
                    elseif strcmp(this.uitgAxesDisplay.getSelectedTabName(),'Reconstruction')
                        % Plot wavefronts on phase tab
                        this.replot(this.U8RECONSTRUCTION, []);
                    end
                end
            else
                % iterations for reconstruction
                switch cAnalysisDomain % reconstruction method
                    case 'rPIE' % scanning solution
                        % set gpu array
                        gpuObject = gpuArray(this.dObjectRecon);
                        gpuProbe = gpuArray(this.dProbeRecon);
                        gpuSqrtInt = gpuArray(sqrtInt);
                        for i =1:iteration
                            tempError = 0;
                            for j =1: nInt
                                reconBox = gpuObject(Rpix(j,1)+[1:N],Rpix(j,2)+[1:N]);
                                exitWave = reconBox.*gpuProbe;
                                detectorWave = PIE.utils.postPropagate (exitWave,propagator,H{1},1);
                                correctedWave = gpuSqrtInt(:,:,j).*detectorWave./(abs(detectorWave)+eps);
                                exitWaveNew = PIE.utils.postPropagate (correctedWave,propagator,Hm{1},1);
                                tempProbe = gpuProbe;
                                denomO = gamma*max(abs(tempProbe(:)).^2) + (1-gamma)*abs(tempProbe).^2;
                                newReconBox = reconBox + alpha*conj(tempProbe).*(exitWaveNew-exitWave)./denomO;
                                denomP = gamma*max(abs(reconBox(:)).^2) + (1-gamma).*abs(reconBox).^2;
                                gpuProbe = gpuProbe + beta*conj(reconBox).*(exitWaveNew-exitWave)./denomP;
                                gpuObject(Rpix(j,1)+[1:N],Rpix(j,2)+[1:N]) = newReconBox;
                                tempError = tempError + abs(gpuSqrtInt(:,:,j)-abs(detectorWave)).^2;
                            end
                            % error evaluation
                            errors(i) =gather( sum(tempError(:)))/totalI;
                            
                            iterationStr = sprintf('%d iterations finished,residual error: %0.5f',i,errors(i));
                            this.uitIteration.set(iterationStr);drawnow;
                            if stopSign==1||(i>1&&((abs(errors(i)-errors(i-1))<0e-7)||errors(i)<this.uieAccuracy.get()))
                                this.dError = errors;
                                % gather data from gpu to memory
                                this.dObjectRecon = gather(gpuObject);
                                this.dProbe = gather(gpuProbe);
                                % Make phase tab active:
                                this.uitgAxesDisplay.selectTabByIndex(this.U8RECONSTRUCTION);
                                % Plot wavefronts on phase tab
                                this.replot(this.U8RECONSTRUCTION, []);
                                break;
                            elseif strcmp(this.uitgAxesDisplay.getSelectedTabName(),'Reconstruction')
                                % gather data from gpu to memory
                                this.dObjectRecon = gather(gpuObject);
                                this.dProbe = gather(gpuProbe);
                                % Plot wavefronts on phase tab
                                this.replot(this.U8RECONSTRUCTION, []);
                            end
                        end
                end
            end
            
            % Set state:
            this.setState(this.U8STATE_PHASE_PROCESSED);
        end
        
        
        function analyze(this)
            % load parameters
            u8ModeId = this.uilSelectMode.getSelectedIndexes();
            selectedObject = this.uipSelectObject.getOptions{this.uipSelectObject.getSelectedIndex()};
            selectedRegion = this.uipSelectRegion.getOptions{this.uipSelectRegion.getSelectedIndex()};
            FP = this.uicbFourierPtychography.get();
            probeType = this.uipProbeType.getOptions{this.uipProbeType.getSelectedIndex()};
            N           = this.uieRes.get();
            NA           = this.uieNA.get();
            Rc_um   = this.uieRprobe.get()*1000;
            df_um       = this.uiez1.get()*1000;% negative sign corresponds convergent
            z_um       = this.uiez2.get()*1000;
            Magnification = this.uieMag.get();
            lambda_um = this.dLambda_um(u8ModeId);
            % calculate probe radius
            if ~FP
                switch probeType
                    case 'Defocus wave'
                        Rprobe_um = abs(df_um)*tan(asin(NA));
                    case 'Plane wave'
                        Rprobe_um = Rc_um;
                end
                RAiry_um = 0.61*lambda_um/NA;
                if Rprobe_um<RAiry_um
                    Rprobe_um = RAiry_um;
                end
            else
                Rprobe_um = z_um*tan(asin(NA));
            end
            Rpixs = round(Rprobe_um/this.do_um(u8ModeId));
            if FP
                objectRecon = fftshift(fft2(fftshift(this.dObjectRecon(:,:,u8ModeId))));
                probeRecon = fftshift(fft2(fftshift(this.dProbeRecon(:,:,u8ModeId))));
                objectSpectrum = this.dObjectRecon(:,:,u8ModeId);
                probeSpectrum = this.dProbeRecon(:,:,u8ModeId);
                objectOri = fftshift(fft2(fftshift(this.dObject(:,:,u8ModeId))));
                probeOri = fftshift(fft2(fftshift(this.dProbe(:,:,u8ModeId))));
            else
                objectRecon = this.dObjectRecon(:,:,u8ModeId);
                probeRecon = this.dProbeRecon(:,:,u8ModeId);
                objectSpectrum = fftshift(fft2(fftshift(this.dObjectRecon(:,:,u8ModeId))));
                probeSpectrum = fftshift(fft2(fftshift(this.dProbeRecon(:,:,u8ModeId))));
                objectOri = this.dObject(:,:,u8ModeId);
                probeOri = this.dProbe(:,:,u8ModeId);
            end
            %              objectRecon = circshift(objectRecon,[1,1]);
            switch selectedObject
                case 'Object amplitude'
                    this.dSelectedObject = abs(objectRecon);
                    if FP
                        this.dUnit_mm = this.dc_um/1000/Magnification ;
                    else
                        this.dUnit_mm = this.do_um(u8ModeId)/1000;
                    end
                case 'Object phase'
                    this.dSelectedObject = atan2(imag(objectRecon),real(objectRecon));
                    if FP
                        this.dUnit_mm = this.dc_um/1000/Magnification ;
                    else
                        this.dUnit_mm = this.do_um(u8ModeId)/1000;
                    end
                case 'Probe amplitude'
                    this.dSelectedObject = abs(probeRecon);
                    if FP
                        this.dUnit_mm = this.dc_um/1000/Magnification ;
                    else
                        this.dUnit_mm = this.do_um(u8ModeId)/1000;
                    end
                case 'Probe phase'
                    this.dSelectedObject = atan2(imag(probeRecon),real(probeRecon));
                    if FP
                        this.dUnit_mm = this.dc_um/1000/Magnification ;
                    else
                        this.dUnit_mm = this.do_um(u8ModeId)/1000;
                    end
                case 'Object spectrum amplitude'
                    this.dSelectedObject = abs(objectSpectrum);
                    if FP
                        this.dUnit_mm = this.do_um(u8ModeId)/1000;
                    else
                        this.dUnit_mm = this.dc_um/1000;
                    end
                case 'Object spectrum phase'
                    this.dSelectedObject = atan2(imag(objectSpectrum),real(objectSpectrum));
                    if FP
                        this.dUnit_mm = this.do_um(u8ModeId)/1000;
                    else
                        this.dUnit_mm = this.dc_um/1000;
                    end
                case 'Probe spectrum amplitude'
                    this.dSelectedObject = abs(probeSpectrum);
                    if FP
                        this.dUnit_mm = this.do_um(u8ModeId)/1000;
                    else
                        this.dUnit_mm = this.dc_um/1000;
                    end
                case 'Probe spectrum phase'
                    this.dSelectedObject = atan2(imag(probeSpectrum),real(probeSpectrum));
                    if FP
                        this.dUnit_mm = this.do_um(u8ModeId)/1000;
                    else
                        this.dUnit_mm = this.dc_um/1000;
                    end
                case 'Object amplitude difference'
                    this.dSelectedObject = abs(objectRecon)./mean(abs(objectRecon(this.dAnalysisMask==1))) - abs(objectOri)./mean(abs(objectOri(this.dAnalysisMask==1)));
                    if FP
                        this.dUnit_mm = this.dc_um/1000/Magnification*N/length(objectRecon) ;
                    else
                        this.dUnit_mm = this.do_um(u8ModeId)/1000;
                    end
                case 'Object phase difference'
                    %                     this.dSelectedObject = atan2(imag(objectRecon),real(objectRecon))-...
                    %                         atan2(imag(objectOri),real(objectOri));
                    s=objectRecon./objectOri;
                    this.dSelectedObject = atan2(imag(s),real(s));
                    dph = 1.2;
                    %                                         this.dSelectedObject(this.dSelectedObject<min(min(this.dSelectedObject))+dph)=this.dSelectedObject(this.dSelectedObject<min(min(this.dSelectedObject))+dph)+2*pi;
                    %                                         this.dSelectedObject(this.dSelectedObject>max(max(this.dSelectedObject))-dph)=this.dSelectedObject(this.dSelectedObject>max(max(this.dSelectedObject))-dph)-2*pi;
                    this.dSelectedObject =this.dSelectedObject- mean(mean(this.dSelectedObject));
                    %                     this.dSelectedObject =PIE.utils.DelTilt(this.dSelectedObject);
                    if FP
                        this.dUnit_mm = this.dc_um/1000/Magnification*N/length(objectRecon) ;
                    else
                        this.dUnit_mm = this.do_um(u8ModeId)/1000;
                    end
                case 'Probe amplitude difference'
                    this.dSelectedObject = abs(probeRecon)./mean(abs(probeRecon(pinhole(2*Rpixs,N,N)==1))) - abs(probeOri)./mean(abs(probeOri(pinhole(2*Rpix,N,N)==1)));
                    if FP
                        this.dUnit_mm = this.dc_um/1000/Magnification ;
                    else
                        this.dUnit_mm = this.do_um(u8ModeId)/1000;
                    end
                case 'Probe phase difference'
                    %                     this.dSelectedObject = atan2(imag(probeRecon),real(probeRecon))-...
                    %                         atan2(imag(probeOri),real(probeOri));
                    s=probeRecon./probeOri;
                    this.dSelectedObject = atan2(imag(s),real(s));
                    this.dSelectedObject =this.dSelectedObject- mean(mean(this.dSelectedObject));
                    if FP
                        this.dUnit_mm = this.dc_um/1000/Magnification ;
                    else
                        this.dUnit_mm = this.do_um(u8ModeId)/1000;
                    end
            end
            
            [K,L] = size(this.dSelectedObject);
            if FP
                k0 = 2*pi/lambda_um;
                %                 kxy = k0*sin(atan(this.dPos_mm*1000/z_um));
                kr = k0*sin(atan(sqrt(this.dPos_mm(:,1).^2+this.dPos_mm(:,2).^2)/this.uieLo.get()));
                phi = atan2(this.dPos_mm(:,1),this.dPos_mm(:,2));
                kxy = [kr.*sin(phi),kr.*cos(phi)];
                dkxy = 2*pi/this.dc_um/N*this.uieMag.get();
                Rpix = kxy./dkxy;
                Rpix = Rpix -min(Rpix,[],1);
                Rpix = round(Rpix);
                %                 [K,L] = size(this.dObjectGuess);
            else
                Rpix = round((this.dPos_mm(:,1:2)-min(this.dPos_mm(:,1:2),[],1))*1000/this.do_um(u8ModeId));
            end
            
            switch selectedRegion
                case 'Entire region'
                    this.dAnalysisMask = ones(K,L);
                case 'Compute from scanning positions'
                    this.dAnalysisMask = zeros(K,L);
                    Rminm = min(Rpix(:,1));
                    Rminn = min(Rpix(:,2));
                    Rmaxm = max(Rpix(:,1));
                    Rmaxn = max(Rpix(:,2));
                    probeRegion = pinhole(2*Rpixs,N,N);
                    if K>N
                        this.dAnalysisMask = zeros(K,L);
                        this.dAnalysisMask(Rminm+N/2+1:Rmaxm+N/2,Rminn+N/2+1:Rmaxn+N/2) = 1;
                    else
                        this.dAnalysisMask = probeRegion;
                    end
                case 'Compute from probe diameter'
                    this.dAnalysisMask = zeros(K,L);
                    probeRegion = pinhole(2*Rpixs,N,N);
                    if K>N
                        for i = 1:length(this.dPos_mm)
                            this.dAnalysisMask(Rpix(i,1)+[1:N],Rpix(i,2)+[1:N]) = ...
                                probeRegion + this.dAnalysisMask(Rpix(i,1)+[1:N],Rpix(i,2)+[1:N]);
                        end
                        this.dAnalysisMask(this.dAnalysisMask~=0)=1;
                    else
                        this.dAnalysisMask = probeRegion;
                    end
                case 'Load region'
            end
            sigma = this.uieSigma.get();
            if sigma >0 %% using sigma filter remove noise
                selectedData = this.dSelectedObject(this.dAnalysisMask==1);
                rms = std(selectedData);
                selectedData(abs(selectedData-mean(selectedData))>sigma*rms)=NaN;
                this.dSelectedObject(this.dAnalysisMask==1) = selectedData;
            end
            % Make phase tab active:
            this.uitgAxesDisplay.selectTabByIndex(this.U8ANALYSIS);
            
            % Plot wavefronts on phase tab
            this.replot(this.U8ANALYSIS, []);
            
            % Set state:
            this.setState(this.U8STATE_ANALYZED);
        end
        
        function SaveAnalysisLog(this)
            filename=this.uieLogFileNameStack.get();
            if strcmp( filename(end-10:end),'scanlog.csv')
                filename=strcat(filename(1:end-11),'analysislog.csv');
                [fid, isNewLogFile] = this.openOrCreateFile( filename);
                cWriteStr = '';
                nl = java.lang.System.getProperty('line.separator').char;
                ceFieldNames=['Lambda/um','T/um','NA','z2/mm','Gr. Tilt/Deg','Dt. Tilt/Deg','CCD Rot/Deg','Dt. Size'...
                    'Center Obs.','Binning','Data format','CCD Center','Obs. Offset','Mask Selection','AnalysisDomain',...
                    'UnwrapEngine','z1/um','PhaseSteps','FT Type','Filter Width','Filter Type','ReconstructionType',...
                    'Number of Zernikes','Rimmer Type','Auto Load Rimmer','Auto Load DB','Auto Load FT','Zernikes Basis Number',...
                    'Remove X Tilt','Remove Y Tilt','Remove Defocus','Scaled Null Wave',this.htLog.ColumnName'];
                if isNewLogFile
                    for k = 1:length(ceFieldNames)
                        cWriteStr = sprintf('%s%s,',cWriteStr, ceFieldNames{k});
                    end
                    cWriteStr(end) = [];
                    cWriteStr = [cWriteStr nl];
                end
                AnalysisLog=[this.ceAnalysisPara,this.ceAnalysisTable{end,:}];
                % Write structure fields
                for k = 1:length(ceFieldNames)
                    cWriteStr = sprintf('%s%s,', cWriteStr, num2str(AnalysisLog{k}));
                end
                cWriteStr(end) = [];
                cWriteStr = [cWriteStr nl];
                fwrite(fid, cWriteStr);
                fclose(fid);
            end
        end
        
        
        function customSimulation(this)
            this.uieCCDCenter.set('[]');
            this.uieObsOffset.set('[]');
            
            N           = this.uieRes.get();
            lambda_um   = this.uieLambda.get()/1000;
            T_um        = this.uieScanRange.get();
            z1_um       = this.uiez1.get()*1000;
            z2_mm       = this.uiez2.get();
            alpha       = this.uieGratTilt.get() * pi/180;
            gamma       = this.uieDetTilt.get() * pi/180;
            zernCouples = eval(this.uieZrn.get());
            NA          = this.uieNA.get();
            Np          = this.uieNp.get();
            detSize     = this.uieDetSize.get();
            CenterObstruction     = this.uieCenterObstruction.get();
            
            SelectMask  = this.uipSelectMask.getSelectedIndex();
            
            % Make phase shifts:
            N = 128;
            npx = 1;
            npy = 100;
            nT = 2;
            
            phIdxX = linspace(0, 2*pi * nT,  npx + 1);
            phIdxY = linspace(0,  2*pi * nT, npy + 1);
            
            phIdxX = phIdxX(1:end-1);
            phIdxY = phIdxY(1:end-1);
            
            [dPhX, dPhY] = meshgrid(phIdxX, phIdxY);
            
            % Rotation
            dCoup = [dPhX(:)';dPhY(:)'];
            dRot = 0; %mrad
            dCoupR = [cos(dRot/1000), -sin(dRot/1000); sin(dRot/1000), cos(dRot/1000)]*dCoup;
            dPhX = reshape(dCoupR(1,:), npy, npx);
            dPhY = reshape(dCoupR(2,:), npy, npx);
            
            
            
            % lets try a z1 shift:
            zxdelta = 0.01;
            zydelta = 0.01;
            
            isX = false;
            
            zidxX = linspace(-zxdelta/2, zxdelta/2, npx);
            zidxY = linspace(-zydelta/2, zydelta/2, npy);
            [Zx, Zy] = meshgrid(zidxX, zidxY);
            
            
            
            
            for k = 1:npx
                for m = 1:npy
                    fprintf('[m, k] = [%d, %d]\n', m, k);
                    Wint = lsianalyze.utils.simulateTwoRayMethod1D(...
                        N, lambda_um, NA, T_um, z1_um + (Zx(m,k) + Zy(m,k)), z2_mm, ...
                        alpha, gamma, zernCouples, Np, detSize,CenterObstruction,SelectMask,[dPhX(m,k), dPhY(m, k)], isX);
                    simInts{k, m} = Wint; %#ok<AGROW>
                    
                    this.uipbExposureProgress.set(m/npy+(k-1/npx));
                    drawnow
                end
                
            end
            
            this.handleLoadData(simInts, {'sim'});
            this.uipbExposureProgress(1);
            
            
        end
        
        
        function simulateInteferograms(this, lComputePSStack)
            % Set center pixel to middle:
            this.uieCCDCenter.set('[]');
            this.uieObsOffset.set('[]');
            
            
            % Generate diffraction patterns:
            %             scanRange_um = this.uieScanRange.get()*1000;
            N           = this.uieRes.get();
            dMaxPhoton  = this.uieNp.get();
            FP = this.uicbFourierPtychography.get();
            s2s         = this.uieShotToShot.get();
            modeNumber = this.uieModeNumber.get();
            lambda_um   = this.uieLambda.get()/1000;
            z_um        =this.uiez2.get()*1000;
            propagator = this.uipPropagator.getOptions{this.uipPropagator.getSelectedIndex()};
            for j = 1:modeNumber
                H{j} = PIE.utils.prePropagate (this.dProbe(:,:,j),propagator,this.do_um(j),this.dLambda_um(j),z_um,1);
            end
            if isempty(this.dProbe)
                this.dProbe = pinhole(round(N/2),N,N);
            end
            if isempty(this.dObject)
                KL = round(this.uieScanRange.get()*1000./this.do_um(1))+N;
                this.dObject = ones(KL);
            end
            
            % load error sources
            dZDriftTotal = this.uieZLinearDrift.get() / 1000; % (um)
            dXDriftTotal = this.uieXLinearDrift.get() / 1000000; % (mm)
            dYDriftTotal = this.uieYLinearDrift.get() / 1000000; % (mm)
            dNonlinearity = this.uieNonlinearity.get()*randn(N);
            dStageUncertainty = this.uiePhaseShiftingError.get() / 1000000; % (mm)
            airflow     = this.uieAirflow.get();
            MSFR        = this.uieMSFR.get();
            dcFlare     =this.uieFlareLevel.get();
            dPosSim_mm = this.dPos_mm;
            xL = max(this.dPos_mm(:,2))-min(this.dPos_mm(:,2));
            dPosSim_mm(:,2) = dPosSim_mm(:,2)+dXDriftTotal*(dPosSim_mm(:,2)-min(dPosSim_mm(:,2)))/xL...
                + dStageUncertainty*randn(size(dPosSim_mm,1),1);
            yL = max(this.dPos_mm(:,1))-min(this.dPos_mm(:,1));
            dPosSim_mm(:,1) = dPosSim_mm(:,1)+dYDriftTotal*(dPosSim_mm(:,1)-min(dPosSim_mm(:,1)))/yL...
                + dStageUncertainty*randn(size(dPosSim_mm,1),1);
            for j = 1:modeNumber
                if FP
                    k0 = 2*pi/this.dLambda_um(j);
                    %                 kxy = k0*sin(atan(this.dPos_mm*1000/z_um));
                    %                 kr = k0*sin(atan(sqrt(this.dPos_mm(:,1).^2+this.dPos_mm(:,2).^2)*1000/z_um));
                    kr = k0*sin(atan(sqrt(dPosSim_mm(:,1).^2+dPosSim_mm(:,2).^2)/this.uieLo.get()));
                    phi = atan2(dPosSim_mm(:,1),dPosSim_mm(:,2));
                    kxy = [kr.*sin(phi),kr.*cos(phi)];
                    %                 dkxy = 2*pi/this.dc_um/N;
                    dkxy = 2*pi/this.dc_um/N*this.uieMag.get();
                    dPosShifts = kxy./dkxy;
                    %                 kxy = round(dPosShifts).*dkxy;
                    dPosShifts = dPosShifts -min(dPosShifts,[],1);
                    dPosShifts = round(dPosShifts);
                else
                    dPosShifts = round((dPosSim_mm(:,1:2)-min(dPosSim_mm(:,1:2),[],1))*1000/this.do_um(j));
                end
                [K,L,~] = size(this.dObject);
                dPosShifts(dPosShifts(:,1)+N>K,1) = K-N;
                dPosShifts(dPosShifts(:,1)<0,1) = 0;
                dPosShifts(dPosShifts(:,2)+N>L,2) = L-N;
                dPosShifts(dPosShifts(:,2)<0,2) = 0;
                if size(this.dPos_mm,2) ==3
                    dPosShifts(:,3) = this.dPos_mm(:,3)*1000;
                    zL = max(dPosShifts(:,3))-min(dPosShifts(:,3));
                    dPosShifts(:,3) = dPosShifts(:,3)+dZDriftTotal*(dPosShifts(:,3)-min(dPosShifts(:,3)))/zL;
                end
                cePosShifts{j} = dPosShifts;
            end
            % add MSFR
            dMSFR_Phase =   PIE.utils.generateMSFN(MSFR,N,10,100); % generate MSFR
            if lComputePSStack % True if simulating "stack"
                this.uipbExposureProgress.set(0);
                nSteps = length(dPosShifts);
                nStep2 = ceil(sqrt(nSteps))^2;
                simInts = cell(nStep2, 1); 
                minValue=0;
                maxValue=0;
                for m = 1:nSteps
                    % add airflow
                    dAirflow_Phase =   PIE.utils.generateMSFN(airflow,N,1,10); % generate airflow
                    dAdditionalPhase = dAirflow_Phase+dMSFR_Phase;
                    %                     if FP
                    for j = 1:modeNumber
                        mPosShifts(j,:) = cePosShifts{j}(m,:);
                    end
                    sqrtInt = PIE.utils.simulateDiffractionPattern(this.dProbe.*exp(1i*2*pi*dAdditionalPhase),...
                        this.dObject,this.ceSegments,modeNumber,N,propagator,mPosShifts,H,1,this.do_um,this.dLambda_um);
                    %                     else
                    %                         [K,L] = size(this.dObject);
                    %                         xo_um = linspace(-L/2,L/2-1,L)*this.dc_um*N/L/this.uieMag.get(); % object coordinates
                    %                         yo_um = linspace(-K/2,K/2-1,K)*this.dc_um*N/K/this.uieMag.get(); % object coordinates
                    %                         [x,y] = meshgrid(xo_um,yo_um);
                    %                         tilt = exp(1i*kxy(m,2).*x+1i*kxy(m,1).*y);
                    %                         object =  PIE.utils.Propagate (this.dObject,propagator,this.do_um,lambda_um,1);
                    %                         object = object.*tilt;
                    %                         object =  PIE.utils.Propagate (object,propagator,this.do_um,lambda_um,-1);
                    %                         sqrtInt = PIE.utils.simulateFPMPattern(this.dProbe,object,this.ceSegments,modeNumber,N,...
                    %                             propagator,dPosShifts(m,:),H,1,this.do_um,lambda_um);
                    %                     end
                    % add systematic error
                    Int = sqrtInt.^2;
                    simInts{m} = PIE.utils.addSystematicError(Int,dMaxPhoton,s2s,dcFlare,this.ceSegments);
                    % calculate the minimum and maximum instensity to normalize
                    % intensities
                    if -min(simInts{m}(:))>minValue
                        minValue=min(simInts{m}(:));
                    end
                    if max(simInts{m}(:))>maxValue
                        maxValue=max(simInts{m}(:));
                    end
                    this.uipbExposureProgress.set(m/nSteps);
                end
                
                % add ccd nonlinearity
                if minValue<0
                    for m = 1:nSteps
                        simInts{m}=((simInts{m}-minValue)/(maxValue-minValue)).^(1+dNonlinearity);
                    end
                else
                    for m = 1:nSteps
                        simInts{m}=(simInts{m}/maxValue).^(1+dNonlinearity);
                    end
                end
                
                
                % make sure data set is N by N
                for m = nSteps+1:nStep2
                    simInts{m} = zeros(N);
                end
                %
                simInts =reshape(simInts,sqrt(nStep2),sqrt(nStep2));
                
                this.handleLoadData(simInts, {'sim'});
                this.dAnalysisRegion(simInts{1,1}==0)=0;
                this.uipbExposureProgress(1);
                drawnow
            else % Simulate a singe image:
                tic
                for j = 1:modeNumber
                    mPosShifts(j,1:2) = cePosShifts{j}(1,1:2);
                end
                % add airflow
                dAirflow_Phase =   PIE.utils.generateMSFN(airflow,N,1,10); % generate airflow
                dAdditionalPhase = dAirflow_Phase+dMSFR_Phase;
                %
                sqrtInt = PIE.utils.simulateDiffractionPattern(this.dProbe.*exp(1i*2*pi*dAdditionalPhase),this.dObject,...
                    this.ceSegments,modeNumber,N,propagator,mPosShifts,H,1,this.do_um,this.dLambda_um);
                % add systematic error
                Int = sqrtInt.^2;
                Int = PIE.utils.addSystematicError(Int,dMaxPhoton,s2s,dcFlare,this.ceSegments);
                
                % add ccd nonlinearity
                if min(Int(:))<0
                    Int=((Int-min(Int(:)))/(max(Int(:))-min(Int(:)))).^(1+dNonlinearity);
                else
                    Int=(Int/max(Int(:))).^(1+dNonlinearity);
                end
                this.handleLoadData({Int}, {'sim'});
                this.dAnalysisRegion(sqrtInt==0)=0;
                fprintf('Single image took %s\n', s2f(toc));
            end
            
            % load segment mask
            if ~isempty(this.ceSegments)
                segMask = zeros(length(this.ceSegments{1}));
                for k=1:length(this.ceSegments)
                    segMask= segMask+this.ceSegments{k};
                end
                %                 this.dAnalysisRegion=ones(length(segs{1}));
                this.dAnalysisRegion(segMask==0)=0;
            end
            
            % Set state:
            this.setState(this.U8STATE_DATA_LOADED);
            if this.uicbAutoCenter.get()
                dCenter = lsianalyze.utils.autoCenter(this.ceInt);
                this.uieCCDCenter.set(['[',num2str(dCenter),']']);
                % redraw guide lines
                this.replot(this.U8DATA, []);
                % Need to recenter data
                this.handleLoadData();
            end
        end
        
        function computeDefaultMaskGeometry(this)
            lambda_um   = this.uieLambda.get()/1000;
            T_um        = this.uieScanRange.get();
            z2_mm       = this.uiez2.get();
            z1_mm       = this.uiez1.get();
            NA          = this.uieNA.get();
            detSize     = this.uieDetSize.get();
            dObsOffset = eval(this.uieObsOffset.get());
            CenterObstruction     = this.uieCenterObstruction.get();
            [sr, sc]    = size(this.ceInt{1});
            
            if NA>0&&~this.uicbFourierPtychography.get()
                this.dBeamWidthEstPx    = (tan(asin(NA))*(z2_mm) * sr / (detSize));
            else
                this.dBeamWidthEstPx    = this.uieRprobe.get()/detSize*sr;
            end
            
            this.dObstructionWidthEstPx = (tan(asin(CenterObstruction*NA))*(z2_mm) * sr / (detSize));
            %remove lines
            this.dAnalysisRegion = ones(sr,sc);
            this.dAnalysisRegion2=ones(sr,sc);
            
        end
        
        function computeMET5Mask(this)
            lambda_um   = this.uieLambda.get()/1000;
            T_um        = this.uieScanRange.get();
            z2_mm       = this.uiez2.get();
            z1_mm       = this.uiez1.get();
            NA          = this.uieNA.get();
            detSize     = this.uieDetSize.get();
            dObs     = this.uieCenterObstruction.get();
            [sr, sc]    = size(this.ceInt{1});
            
            idx = linspace(-detSize/2, detSize/2, sr);
            [dX, dY] = meshgrid(idx);
            
            dOutBeam    =  tan(asin(NA))*z1_mm+tan(asin(NA - lambda_um/T_um))*(z2_mm-z1_mm);
            dInBeam     =  tan(asin(NA * dObs) + lambda_um/T_um)*(z2_mm);
            
            
            dROI = zeros(sr, sc);
            
            dROI(dX.^2 + dY.^2 > dInBeam^2 & dX.^2 + dY.^2 < dOutBeam^2) = 1;
            
            dShear          = tan(asin(NA)) * z2_mm - dOutBeam;
            
            
            this.dAnalysisRegion = dROI;
            
            dROI( dY > 0 & (abs(dX - dY) < dShear | abs(dX + dY) < dShear)) = 0;
            
            this.dAnalysisRegion2 = dROI;
            
            
        end
        function computeElipticalMask(this)
            lambda_um   = this.uieLambda.get()/1000;
            T_um        = this.uieScanRange.get();
            z2_mm       = this.uiez2.get();
            z1_mm       = this.uiez1.get();
            NA          = this.uieNA.get();
            detSize     = this.uieDetSize.get();
            dObsOffset = eval(this.uieObsOffset.get());
            CenterObstruction     = this.uieCenterObstruction.get();
            [sr, sc]    = size(this.ceInt{1});
            this.dBeamWidthEstPx    = (tan(asin(NA))*(z2_mm) * sr / (detSize));
            this.dShearPix          = (this.dBeamWidthEstPx - ...
                (tan(asin(NA))*z1_mm+tan(asin(NA - lambda_um/T_um))*(z2_mm-z1_mm)) * sr / (detSize));
            this.dObstructionWidthEstPx = (tan(asin(CenterObstruction*NA))*(z2_mm) * sr / (detSize));
            this.dObstructionShearPix =(tan(asin(lambda_um/T_um))*z2_mm* sr / detSize);
            % Set default mask:
            if mod(sr,2)==0
                if this.dObstructionWidthEstPx~=0
                    if ~isempty(dObsOffset)
                        this.dAnalysisRegion = lsianalyze.utils.elipticalHole(2*round(this.dBeamWidthEstPx - this.dShearPix),2*round(0.5714*(this.dBeamWidthEstPx - this.dShearPix)), sr, sc)-...
                            circshift(lsianalyze.utils.elipticalHole(2*round(this.dObstructionWidthEstPx + this.dObstructionShearPix),2*round(0.5714*(this.dObstructionWidthEstPx + this.dObstructionShearPix)), sr, sc), ...
                            [ dObsOffset(2),  dObsOffset(1)]);
                    else
                        this.dAnalysisRegion = lsianalyze.utils.elipticalHole(2*round(this.dBeamWidthEstPx - this.dShearPix),2*round(0.5714*(this.dBeamWidthEstPx - this.dShearPix)), sr, sc)-...
                            lsianalyze.utils.elipticalHole(2*round(this.dObstructionWidthEstPx + this.dObstructionShearPix),2*round(0.5714*(this.dObstructionWidthEstPx + this.dObstructionShearPix)), sr, sc);
                    end
                else
                    this.dAnalysisRegion = lsianalyze.utils.elipticalHole(2*round(this.dBeamWidthEstPx - this.dShearPix), 2*round(0.5714*(this.dBeamWidthEstPx - this.dShearPix)),sr, sc);
                end
            else
                if this.dObstructionWidthEstPx~=0
                    if ~isempty(dObsOffset)
                        this.dAnalysisRegion = lsianalyze.utils.elipticalHole(2*ceil(this.dBeamWidthEstPx - this.dShearPix)-1,2*ceil(0.5714*(this.dBeamWidthEstPx - this.dShearPix))-1, sr, sc)-...
                            circshift(lsianalyze.utils.elipticalHole(2*floor(this.dObstructionWidthEstPx + this.dObstructionShearPix)+1,2*floor(0.5714*(this.dObstructionWidthEstPx + this.dObstructionShearPix))+1, sr, sc), ...
                            [ dObsOffset(2),  dObsOffset(1)]);
                    else
                        this.dAnalysisRegion = lsianalyze.utils.elipticalHole(2*ceil(this.dBeamWidthEstPx - this.dShearPix)-1,2*ceil(0.5714*(this.dBeamWidthEstPx - this.dShearPix))-1, sr, sc)-...
                            lsianalyze.utils.elipticalHole(2*floor(this.dObstructionWidthEstPx + this.dObstructionShearPix)+1,2*floor(0.5714*(this.dObstructionWidthEstPx + this.dObstructionShearPix))+1, sr, sc);
                    end
                else
                    this.dAnalysisRegion = lsianalyze.utils.elipticalHole(2*ceil(this.dBeamWidthEstPx - this.dShearPix)-1, 2*ceil(0.5714*(this.dBeamWidthEstPx - this.dShearPix))-1,sr, sc);
                end
            end
            
            %remove lines
            this.dAnalysisRegion2=ones(sr,sc);
            
        end
        
        % Position recall getters
        function store = prGetters(this, src)
            store = [];
            switch src
                case this.prControlsSetup
                    ceFields = {'uieLambda', 'uieScanRange', 'uieNA', 'uiez2', ...
                        'uieGratTilt', 'uieDetTilt', 'uieGlobalRot', 'uieDetSize',...
                        'uieCenterObstruction','uieBinning','uiez1','uieScanSteps'};
                    for k = 1:length(ceFields)
                        store(end + 1) = this.(ceFields{k}).get(); %#ok<AGROW>
                    end
                    
                    % Now get center:
                    dCenter = eval(this.uieCCDCenter.get());
                    if numel(dCenter) == 0
                        dX = 0;
                        dY = 0;
                    else
                        dX = dCenter(1);
                        dY = dCenter(2);
                    end
                    store(end+1) =  dX;
                    store(end+1) =  dY;
                    
            end
            
        end
        
        % Position recall setters
        function prSetters(this, src, dRecall)
            switch src
                case this.prControlsSetup
                    ceFields = {'uieLambda', 'uieScanRange', 'uieNA', 'uiez2', ...
                        'uieGratTilt', 'uieDetTilt', 'uieGlobalRot', 'uieDetSize',...
                        'uieCenterObstruction','uieBinning','uiez1','uieScanSteps'};
                    for k = 1:length(ceFields)
                        this.(ceFields{k}).set(dRecall(k));
                    end
                    
                    % Now set center:
                    dCenter = dRecall(end-1:end);
                    if ~all(dCenter==0)
                        this.uieCCDCenter.set(sprintf('[%d,%d]', dCenter(1), dCenter(2)));
                    else
                        this.uieCCDCenter.set('[]');
                    end
                    % update scanning position
                    this.cb(this.uieScanRange);
                    this.cb(this.uieNA);
            end
            
            % Need to reload data after a parameter change:
            this.handleLoadData();
        end
        
        % Call this when data is loaded.  This will overwrite any current
        % data.  If cell data is Nx2, assume first col is X and 2nd is Y
        function handleLoadData(this, ceData, ceMeta)
            Binning=this.uieBinning.get();
            % if load data is called empty, then reset back to original
            % images
            if (nargin == 1)
                ceData = this.ceInt_0;
                ceMeta = this.ceIntMeta;
            else
                % Clear ints:
                this.ceInt_0 = ceData;
                if isempty(ceMeta) || strcmp(ceMeta{1},'shifted')~=1
                    this.ceInt  = cell(size(this.ceInt_0));
                end
            end
            [sr,sc]=size(ceData);
            for k = 1:sr
                for t=1:sc
                    if ~isempty(ceMeta) && strcmp(ceMeta{1},'sim')==1
                        % ceMeta={};
                    elseif ~isempty(ceMeta) && strcmp(ceMeta{1},'shifted')==1
                        ceData{k,t}=bin2(this.ceInt_0{k,t},Binning,Binning);
                    else
                        ceData{k,t}=bin2(ceData{k,t},Binning,Binning);
                    end
                    if length(ceMeta(:)) >= k*t
                        this.ceIntMeta{k,t} = ceMeta{k,t};
                    else
                        this.ceIntMeta{k,t} = struct();
                    end
                end
            end
            
            %[sr, sc] = size(ceData);
            this.u8NumInt       = sr * sc;
            this.u8ActiveIntIdx = 1;
            
            dGlobalRot = this.uieGlobalRot.get();
            
            % Process ints:
            dNumInts = sr * sc;
            for k = 1:sr
                for t = 1:sc
                    fprintf('Processing ints %d of %d\n', (k-1)*sc + t, dNumInts);
                    % this.ceInt_0{k}=ceData{k};
                    
                    if (dGlobalRot ~= 0)
                        % Rotate patterns here
                        this.ceIntR{k,t} = imrotate(ceData{k,t},dGlobalRot,'crop');
                    else
                        this.ceIntR{k,t} = ceData{k,t};
                    end
                    if ~isempty(ceMeta) && strcmp(ceMeta{1},'shifted')==1
                        this.ceIntR{k,t} = this.ceInt{k,t};
                        this.ceIntMeta={};
                    end
                    % normalize intensity
                    %                     if strcmp(ceMeta{1},'sim')~=1
                    FP = this.uicbFourierPtychography.get();
                    if k==1 && t==1
                        if ~isempty(this.ceSegments) % normalize intensity for segments detector
                            segInt=zeros(length(this.ceSegments),1);
                            detectorUsage = 0;
                            for ns=1:length(this.ceSegments)
                                seg = logical(this.ceSegments{ns});
                                segInt(ns) =mean(this.ceIntR{1,1}(seg));
                                detectorUsage = detectorUsage+sum(seg(:));
                            end
                            detectorUsage = detectorUsage/length(this.ceIntR{1,1})^2;
                            if ~FP
                                normalizeFactor = sum(abs(this.dProbeGuess(:)).^2)/sum(segInt)*detectorUsage;
                            else
                                normalizeFactor = sum(abs(this.dObjectGuess(:)).^2)/sum(segInt)*detectorUsage;
                            end
                        else
                            if ~FP
                                normalizeFactor = sum(abs(this.dProbeGuess(:)).^2)/sum(this.ceIntR{1,1}(:));
                            else
                                normalizeFactor = sum(abs(this.dObjectGuess(:)).^2)/sum(this.ceIntR{1,1}(:));
                            end
                        end
                    end
                    this.ceIntR{k,t} = this.ceIntR{k,t}*normalizeFactor;
                    %                     end
                    
                    this.ceInt{k,t} = this.ceIntR{k,t};
                    
                    % circshift if necessary:
                    dCenter = eval(this.uieCCDCenter.get());
                    [srs, scs] = size(this.ceIntR{k,t});
                    if (~isempty(dCenter))
                        this.ceInt{k,t} = circshift(this.ceIntR{k,t}, ...
                            [(round(srs/2) - dCenter(2)), round(scs/2) - dCenter(1)]);
                        dMetaFlags = [];
                    else
                        dMetaFlags=this.U8PLOTUNSHIFTED;
                    end
                    
                end
            end
            %             try
            % Compute analysis region and shear region
            SelectMask=this.uipSelectMask.getSelectedIndex();
            switch SelectMask
                case 1
                    this.computeDefaultMaskGeometry();
                case 2
                    this.computeMET5Mask();
                case 3
                    this.computeElipticalMask();
            end
            
            this.replot(this.U8DATA, dMetaFlags);
            %             catch
            %             end
            
            % Make data tab active:
            this.uitgAxesDisplay.selectTabByIndex(this.U8DATA);
            
            % Set state:
            this.setState(this.U8STATE_DATA_LOADED);
            
        end
        
        function cycleInterferogram(this, dDir)
            [sr, sc] = size(this.ceInt);
            N = sr*sc;
            this.u8ActiveIntIdx = this.u8ActiveIntIdx + dDir;
            if this.u8ActiveIntIdx == 0
                this.u8ActiveIntIdx = N;
            elseif this.u8ActiveIntIdx == N + 1
                this.u8ActiveIntIdx = 1;
            end
            this.replot(this.U8DATA, []);
        end
        
        
        % Main redraw function. Pass tab indices to refresh axes
        function replot(this, dTabIdx, dMetaFlags)
            
            for k = 1:length(dTabIdx)
                switch dTabIdx
                    
                    case this.U8DATA
                        [sr, sc] = size(this.ceInt);
                        
                        % Plot main pattern
                        this.hsaInterferogram.setHoldState('off');
                        dObsOffset = eval(this.uieObsOffset.get());
                        
                        if any(dMetaFlags == this.U8PLOTUNSHIFTED)
                            this.hsaInterferogram.imagesc(this.ceIntR{this.u8ActiveIntIdx});
                        else
                            this.hsaInterferogram.imagesc(this.ceInt{this.u8ActiveIntIdx});
                        end
                        
                        % Update interferogram number display
                        if (sc >= 2)
                            if this.u8ActiveIntIdx <= sr*sc/2
                                cSeries = 'X:';
                            else
                                cSeries = 'Y:';
                            end
                            this.uitStack.set(sprintf('%s%d/%d', cSeries, mod(this.u8ActiveIntIdx - 1, sr) + 1, sr));
                        else
                            this.uitStack.set(sprintf('%d/%d', this.u8ActiveIntIdx, this.u8NumInt));
                        end
                        
                        % Build Auxillary graphics
                        if (all(dMetaFlags ~= this.U8AUXPLOTOFF) && all(dMetaFlags ~= this.U8DOMAIN_FFT))
                            [sr, sc] = size(this.ceInt{1});
                            this.hsaInterferogram.setHoldState('on');
                            % Always draw guide lines:
                            
                            dMidX = (sr + 1)/2;
                            dMidY = (sc + 1)/2;
                            this.hsaInterferogram.plot([1, sc], dMidY*[1, 1], 'm');
                            this.hsaInterferogram.plot(dMidX*[1,1], [1, sr], 'm');
                            if this.uipSelectMask.getSelectedIndex()~=3
                                % Build shear guides:
                                dTh = linspace(0, 2*pi, 41);
                                dXCirc1 = this.dBeamWidthEstPx * cos(dTh) + dMidX;
                                dYCirc1 = this.dBeamWidthEstPx * sin(dTh)*1 + dMidY;
                                
                                this.hsaInterferogram.plot(dXCirc1, dYCirc1, 'g', 'lineWidth', 2);
                                
                                if isempty(dObsOffset)
                                    dObsOffset=[0,0];
                                end
                                if this.uieCenterObstruction.get() > 0
                                    dXCirc3 = this.dObstructionWidthEstPx * cos(dTh) + dMidX+dObsOffset(1);
                                    dYCirc3 = this.dObstructionWidthEstPx * sin(dTh)*1 + dMidY+dObsOffset(2);
                                    
                                    this.hsaInterferogram.plot(dXCirc3, dYCirc3, 'r', 'lineWidth', 2);
                                    
                                end
                            else
                                dTh = linspace(0, 2*pi, 41);
                                dXCirc1 = this.dBeamWidthEstPx * cos(dTh) + dMidX;
                                dYCirc1 = this.dBeamWidthEstPx * sin(dTh)*0.5714 + dMidY;
                                
                                this.hsaInterferogram.plot(dXCirc1, dYCirc1, 'g', 'lineWidth', 2);
                                
                                if isempty(dObsOffset)
                                    dObsOffset=[0,0];
                                end
                                if this.uieCenterObstruction.get() > 0
                                    dXCirc3 = this.dObstructionWidthEstPx * cos(dTh) + dMidX+dObsOffset(1);
                                    dYCirc3 = this.dObstructionWidthEstPx * sin(dTh)*1 + dMidY+dObsOffset(2);
                                    
                                    this.hsaInterferogram.plot(dXCirc3, dYCirc3, 'r', 'lineWidth', 2);
                                    
                                end
                                
                            end
                        end
                        
                        
                        % show meta information if it exists:
                        if length(this.ceIntMeta(:)) >= this.u8ActiveIntIdx&&~strcmp(this.ceIntMeta{1},'sim')
                            % Built meta string:
                            ceFieldNames = fieldnames(this.ceIntMeta{this.u8ActiveIntIdx});
                            cStr = '';
                            for k = 1:length(ceFieldNames)
                                cStr = [cStr ceFieldNames{k} ': ' this.ceIntMeta{this.u8ActiveIntIdx}.(ceFieldNames{k}) '    '];
                            end
                            
                        else
                            cStr = 'No meta information about this image';
                        end
                        
                        this.uitMetaInfo.set(cStr);
                        this.hsaInterferogram.setHoldState('off');
                        
                    case  this.U8PROBEOBJECT
                        propagator = this.uipPropagator.getOptions{this.uipPropagator.getSelectedIndex()};
                        lambda_um =this.uieLambda.get()/1000;
                        u8ModeId = this.uilSelectMode.getSelectedIndexes();
                        N           = this.uieRes.get();
                        [K,L] = size(this.dObject(:,:,u8ModeId));
                        xp_mm = linspace(-N/2,N/2,N)*this.do_um(u8ModeId)/1000; % object coordinates
                        Magnification = this.uieMag.get();
                        if this.uicbFourierPtychography.get()
                            object = PIE.utils.Propagate (this.dObject(:,:,u8ModeId),propagator,this.do_um(u8ModeId),lambda_um,1);
                            xo_mm = linspace(-L/2,L/2,L)*this.dc_um*N/L/1000/Magnification; % object coordinates
                            yo_mm = linspace(-K/2,K/2,K)*this.dc_um*N/K/1000/Magnification; % object coordinates
                        else
                            object = this.dObject(:,:,u8ModeId);
                            xo_mm = linspace(-L/2,L/2,L)*this.do_um(u8ModeId)/1000; % object coordinates
                            yo_mm = linspace(-K/2,K/2,K)*this.do_um(u8ModeId)/1000; % object coordinates
                        end
                        
                        imagesc(this.haProbeAmp, xp_mm,xp_mm,abs(this.dProbe(:,:,u8ModeId)));colorbar(this.haProbeAmp);axis(this.haProbeAmp,'xy');
                        this.haProbeAmp.Title.String = 'Probe amplitude';this.haProbeAmp.XLabel.String = 'mm';this.haProbeAmp.YLabel.String = 'mm';
                        imagesc(this.haProbePha, xp_mm,xp_mm,atan2(imag(this.dProbe(:,:,u8ModeId)),real(this.dProbe(:,:,u8ModeId))));colorbar(this.haProbePha);axis(this.haProbePha,'xy');
                        this.haProbePha.Title.String = 'Probe phase';this.haProbePha.XLabel.String = 'mm';this.haProbePha.YLabel.String = 'mm';
                        imagesc(this.haObjectAmp, xo_mm,yo_mm,abs(object));colorbar(this.haObjectAmp);axis(this.haObjectAmp,'xy');
                        this.haObjectAmp.Title.String = 'Object amplitude';this.haObjectAmp.XLabel.String = 'mm';this.haObjectAmp.YLabel.String = 'mm';
                        imagesc(this.haObjectPha, xo_mm,yo_mm,atan2(imag(object),real(object)));colorbar(this.haObjectPha);axis(this.haObjectPha,'xy');
                        this.haObjectPha.Title.String = 'Object phase';this.haObjectPha.XLabel.String = 'mm';this.haObjectPha.YLabel.String = 'mm';
                    case  this.U8GUESS
                        propagator = this.uipPropagator.getOptions{this.uipPropagator.getSelectedIndex()};
                        %                         detSize_um  = this.uieDetSize.get()*1000;
                        u8ModeId = this.uilSelectMode.getSelectedIndexes();
                        lambda_um =this.uieLambda.get()/1000;
                        %                         z_um =this.uiez2.get()*1000;
                        N           = this.uieRes.get();
                        [K,L] = size(this.dObjectGuess(:,:,u8ModeId));
                        xp_mm = linspace(-N/2,N/2,N)*this.do_um(u8ModeId)/1000; % object coordinates
                        Magnification = this.uieMag.get();
                        if this.uicbFourierPtychography.get()
                            object = PIE.utils.Propagate (this.dObjectGuess(:,:,u8ModeId),propagator,this.do_um(u8ModeId),lambda_um,1);
                            xo_mm = linspace(-L/2,L/2,L)*this.dc_um*N/L/1000/Magnification; % object coordinates
                            yo_mm = linspace(-K/2,K/2,K)*this.dc_um*N/K/1000/Magnification; % object coordinates
                        else
                            object = this.dObjectGuess(:,:,u8ModeId);
                            xo_mm = linspace(-L/2,L/2,L)*this.do_um(u8ModeId)/1000; % object coordinates
                            yo_mm = linspace(-K/2,K/2,K)*this.do_um(u8ModeId)/1000; % object coordinates
                        end
                        imagesc(this.haGuessProbeAmp, xp_mm,xp_mm,abs(this.dProbeGuess(:,:,u8ModeId)));colorbar(this.haGuessProbeAmp);axis(this.haGuessProbeAmp,'xy');
                        this.haGuessProbeAmp.Title.String = 'Guessed probe amplitude';this.haGuessProbeAmp.XLabel.String = 'mm';this.haGuessProbeAmp.YLabel.String = 'mm';
                        imagesc(this.haGuessProbePha, xp_mm,xp_mm,atan2(imag(this.dProbeGuess(:,:,u8ModeId)),real(this.dProbeGuess(:,:,u8ModeId))));colorbar(this.haGuessProbePha);axis(this.haGuessProbePha,'xy');
                        this.haGuessProbePha.Title.String = 'Guessed probe phase';this.haGuessProbePha.XLabel.String = 'mm';this.haGuessProbePha.YLabel.String = 'mm';
                        imagesc(this.haGuessObjectAmp, xo_mm,yo_mm,abs(object));colorbar(this.haGuessObjectAmp);axis(this.haGuessObjectAmp,'xy');
                        this.haGuessObjectAmp.Title.String = 'Guessed object amplitude';this.haGuessObjectAmp.XLabel.String = 'mm';this.haGuessObjectAmp.YLabel.String = 'mm';
                        imagesc(this.haGuessObjectPha, xo_mm,yo_mm,atan2(imag(object),real(object)));colorbar(this.haGuessObjectPha);axis(this.haGuessObjectPha,'xy');
                        this.haGuessObjectPha.Title.String = 'Guessed object phase';this.haGuessObjectPha.XLabel.String = 'mm';this.haGuessObjectPha.YLabel.String = 'mm';
                    case this.U8RECONSTRUCTION
                        % intermedium object and probe
                        propagator = this.uipPropagator.getOptions{this.uipPropagator.getSelectedIndex()};
                        lambda_um =this.uieLambda.get()/1000;
                        u8ModeId = this.uilSelectMode.getSelectedIndexes();
                        %                         z_um =this.uiez2.get()*1000;
                        N           = this.uieRes.get();
                        [K,L] = size(this.dObjectGuess(:,:,u8ModeId));
                        xp_mm = linspace(-N/2,N/2,N)*this.do_um(u8ModeId)/1000; % object coordinates
                        Magnification = this.uieMag.get();
                        if this.uicbFourierPtychography.get()
                            object = PIE.utils.Propagate (this.dObjectRecon(:,:,u8ModeId),propagator,this.do_um(u8ModeId),lambda_um,1);
                            xo_mm = linspace(-L/2,L/2,L)*this.dc_um*N/L/1000/Magnification; % object coordinates
                            yo_mm = linspace(-K/2,K/2,K)*this.dc_um*N/K/1000/Magnification; % object coordinates
                        else
                            object = this.dObjectRecon(:,:,u8ModeId);
                            xo_mm = linspace(-L/2,L/2,L)*this.do_um(u8ModeId)/1000; % object coordinates
                            yo_mm = linspace(-K/2,K/2,K)*this.do_um(u8ModeId)/1000; % object coordinates
                        end
                        object_pha = atan2(imag(object),real(object));
                        %                         object_pha = unwrap(unwrap(object_pha,[],1),[],2);
                        imagesc(this.haReconProbeAmp, xp_mm,xp_mm,abs(this.dProbeRecon(:,:,u8ModeId)));colorbar(this.haReconProbeAmp);axis(this.haReconProbeAmp,'xy');
                        this.haReconProbeAmp.Title.String = 'Reconstructed probe amplitude';this.haReconProbeAmp.XLabel.String = 'mm';this.haReconProbeAmp.YLabel.String = 'mm';
                        imagesc(this.haReconProbePha, xp_mm,xp_mm,atan2(imag(this.dProbeRecon(:,:,u8ModeId)),real(this.dProbeRecon(:,:,u8ModeId))));colorbar(this.haReconProbePha);axis(this.haReconProbePha,'xy');
                        this.haReconProbePha.Title.String = 'Reconstructed probe phase';this.haReconProbePha.XLabel.String = 'mm';this.haReconProbePha.YLabel.String = 'mm';
                        imagesc(this.haReconObjectAmp, xo_mm,yo_mm,abs(object));colorbar(this.haReconObjectAmp);axis(this.haReconObjectAmp,'xy');
                        this.haReconObjectAmp.Title.String = 'Reconstructed object amplitude';this.haReconObjectAmp.XLabel.String = 'mm';this.haReconObjectAmp.YLabel.String = 'mm';
                        imagesc(this.haReconObjectPha, xo_mm,yo_mm,object_pha);colorbar(this.haReconObjectPha);axis(this.haReconObjectPha,'xy');
                        this.haReconObjectPha.Title.String = 'Reconstructed object phase';this.haReconObjectPha.XLabel.String = 'mm';this.haReconObjectPha.YLabel.String = 'mm';
                        drawnow;
                        
                    case this.U8ANALYSIS
                        selectedObject = this.uipSelectObject.getOptions{this.uipSelectObject.getSelectedIndex()};
                        if strcmp(selectedObject,'Scanning position') ==1
                            if size(this.dPos_mm,2)==2
                                hold( this.haAnalysis,'off');
                                h = plot(this.haAnalysis,this.dPos_mm(:,2),this.dPos_mm(:,1),'.');
                                axis(this.haAnalysis,'xy','equal','tight');
                                set(h,'MarkerSize',20);
                                this.haAnalysis.Title.String = selectedObject;
                                this.haAnalysis.FontSize = 12;
                                this.haAnalysis.XLabel.String = 'x/mm';this.haAnalysis.YLabel.String = 'y/mm';
                                if ~isempty(dMetaFlags)
                                    hold( this.haAnalysis,'on');
                                    hs = plot(this.haAnalysis,this.dPos_mm(dMetaFlags,2),this.dPos_mm(dMetaFlags,1),'.');
                                    set(hs,'MarkerSize',40);
                                    hold( this.haAnalysis,'on');
                                end
                            else
                                hold( this.haAnalysis,'off');
                                h = plot3(this.haAnalysis,this.dPos_mm(:,3),this.dPos_mm(:,2),this.dPos_mm(:,1),'.');
                                axis(this.haAnalysis,'xy','square','tight');
                                set(h,'MarkerSize',20);
                                this.haAnalysis.Title.String = selectedObject;
                                this.haAnalysis.FontSize = 12;
                                this.haAnalysis.XLabel.String = 'z/mm';this.haAnalysis.YLabel.String = 'x/mm';this.haAnalysis.ZLabel.String = 'y/mm';
                                if ~isempty(dMetaFlags)
                                    hold( this.haAnalysis,'on');
                                    hs = plot3(this.haAnalysis,this.dPos_mm(dMetaFlags,3),this.dPos_mm(dMetaFlags,2),this.dPos_mm(dMetaFlags,1),'.');
                                    set(hs,'MarkerSize',40);
                                    hold( this.haAnalysis,'on');
                                end
                            end
                            
                        elseif strcmp(selectedObject,'Residual errors') ==1
                            iteration = length(this.dError);
                            h = plot(this.haZernikeDecomp,1:iteration,log10(this.dError));
                            this.haZernikeDecomp.Title.String = selectedObject;
                            set(h,'lineWidth',2);
                            this.haZernikeDecomp.FontSize = 12;
                            this.haZernikeDecomp.XLabel.String = 'Iteration times';this.haZernikeDecomp.YLabel.String = 'Residual errors (log10)';
                        else
%                             u8ModeId = this.uilSelectMode.getSelectedIndexes();
                            object = this.dSelectedObject;
%                                                         object = unwrap(unwrap(object,[],1),[],2);
%                                                          object=PIE.utils.UnwrapPhaseBySortingReliabilityWithMask(object,255*ones(size(object)));
                            [K,L] = size(object);
                            x_mm = this.dUnit_mm*linspace(-L/2,L/2,L);
                            y_mm = this.dUnit_mm*linspace(-K/2,K/2,K);
                            if max(this.dAnalysisMask(:))>1
                                mask = this.dAnalysisMask;
                                mask(mask~=0)=1;
                                object(mask==0) = NaN;
                                step = mean(object(this.dAnalysisMask==2&~isnan(object)))-mean(object(this.dAnalysisMask==1&~isnan(object)));
                                RMSStr = ['Step(rad): ',num2str(step)];
                                imagesc(this.haAnalysis, x_mm,y_mm,object);colorbar(this.haAnalysis);axis(this.haAnalysis,'xy');
                                this.haAnalysis.Title.String = selectedObject; this.haAnalysis.XLabel.String = 'x/mm';this.haAnalysis.YLabel.String = 'y/mm';
                            else
                              
                                object(this.dAnalysisMask==0)=NaN;
                                object =PIE.utils.DelTilt(object);
                                RMSStr = ['RMS(rad): ',num2str(std(object(this.dAnalysisMask==1&~isnan(object))))];
                                %
                                % crop region
                                Nc = sum(this.dAnalysisMask(round(K/2),:))-2;
                                % find centroid
                                NM=regionprops(this.dAnalysisMask,'Centroid');
                                nm=NM.Centroid;
                                nn=round(nm(1))-round(K/2);
                                mm=round(nm(2))-round(L/2);
%                                 object = circshift(object,-[nn,mm]);
                                object = crop2(object,Nc,Nc);
                                x_um = this.dUnit_mm*linspace(-Nc/2,Nc/2,Nc)*1000;
                                y_um = this.dUnit_mm*linspace(-Nc/2,Nc/2,Nc)*1000;
                                imagesc(this.haAnalysis, x_um,y_um,object);colorbar(this.haAnalysis);axis(this.haAnalysis,'xy');
                                this.haAnalysis.Title.String = selectedObject; this.haAnalysis.XLabel.String = 'x/mm';this.haAnalysis.YLabel.String = 'y/mm';
                                this.haAnalysis.XLabel.String = 'x/um';this.haAnalysis.YLabel.String = 'y/um';
                                this.haAnalysis.FontSize = 14;
                            end
                            this.uitRMS.set(RMSStr);
                        end
                        drawnow;
                        
                    case this.U8LOG
                        this.ceAnalysisTable(end+1,1)={datestr(now, 31)};
                        this.ceAnalysisTable(end,2)={this.cSeries};
                        this.ceAnalysisTable(end,3:4)=num2cell(this.dResult);
                        this.ceAnalysisTable(end,5:4+length(this.dZernike))=num2cell(this.dZernike);
                        ColumnName={'Time','Series','RMS','RMSFit'};
                        CZ=num2cell(1:length(this.dZernike));
                        for i=1:length(this.dZernike)
                            CZ{i}=strcat('Z',num2str( CZ{i}));
                        end
                        ColumnName(5:4+length(this.dZernike))=CZ;
                        set(this.htLog,'Data',this.ceAnalysisTable,'ColumnName',ColumnName);
                        if ~iscell(this.htLog.ColumnWidth)
                            this.htLog.ColumnWidth={'Auto'};
                        end
                        this.htLog.ColumnWidth{1}=110;
                        this.htLog.ColumnWidth{2}=70;
                        this.htLog.ColumnWidth(3:4+length(this.dZernike))={50};
                end
            end
            
        end
        
        
        % Given parameters, simulates stack, processes phase, and
        % reconstructs
        function simStackAndReconstruct(this)
            % Simulate stack:
            this.simulateInteferograms(true);
            
            % Process phase:
            this.reconstruct(this.uitgAnalysisDomain.getSelectedTabName());
            
            
        end
        
        function setSimParams(this, varargin)
            for k = 1:2:length(varargin)
                paramName = varargin{k};
                paramValue = varargin{k+1};
                switch paramName
                    case 'nPhSteps'
                        this.uieScanSteps.set(paramValue);
                        dN = paramValue;
                        dPhseString = sprintf('0:2*pi/%d:%d*pi/%d', dN, dN*2 - 2, dN);
                        this.uiePhaseStepsSim.set(sprintf('[%s;%s]''', dPhseString, dPhseString));
                    case 'nPhotons'
                        this.uieNp.set(paramValue);
                    case 'zDrift'
                        this.uieZLinearDrift.set(paramValue);
                    case 'xDrift'
                        this.uieXLinearDrift.set(paramValue);
                    case 'yDrift'
                        this.uieYLinearDrift.set(paramValue);
                    case 'nonlinearity'
                        this.uieNonlinearity.set(paramValue);
                    case 'gratTilt'
                        this.uieGratingTiltErr.set(paramValue);
                    case 'detTilt'
                        this.uieDetectorTiltErr.set(paramValue);
                    case 'shotToShot'
                        this.uieShotToShot.set(paramValue);
                    case '2ndOrderStrength'
                        this.uie2ndOrderStrength.set(paramValue);
                    case '11OrderStrength'
                        this.uie11orderStrength.set(paramValue);
                    case 'xyStageError'
                        this.uiePhaseShiftingError.set(paramValue);
                    case 'simDetectorRes'
                        this.uieRes.set(paramValue);
                    case 'detectorCurve'
                        this.uieDetCurve.set(paramValue);
                    case 'flareLevel'
                        this.uieFlareLevel.set(paramValue);
                    case 'MSFR'
                        this.uieMSFR.set(paramValue);
                    case 'airflow'
                        this.uieAirflow.set(paramValue);
                    case 'gratingPitch'
                        this.uieScanRange.set(paramValue);
                    case 'Z_1'
                        this.uiez1.set(paramValue);
                end
            end
        end
        
        function build(this, hFigure, dOffsetX, dOffsetY)
            if nargin <3
                dOffsetX = 0;
                dOffsetY = 0;
            elseif nargin == 3
                dOffsetY = dOffsetX;
                dOffsetX =  hFigure;
            end
            
            %% build the main window
            if nargin == 2||nargin == 4
                this.hFigure = hFigure;
            else
                this.hFigure = figure(...
                    'name', 'PIE analysis GUI v1.20190129',...
                    'Units', 'pixels',...
                    'Position', [5 - dOffsetX, 5 - dOffsetY,  this.dWidth, this.dHeight],...
                    'handlevisibility','off',... %out of reach gcf
                    'numberTitle','off',...
                    'Toolbar','none',...
                    'Menubar','none');
            end
            
            
            % Build all containers first:
            
            
            drawnow
            
            %% Axes
            dTgPx = 20;
            dTgPy = 20;
            this.uitgAxesDisplay.build(this.hFigure, dTgPx, dTgPy, 950, 950);
            
            % Axes:Data
            uitData = this.uitgAxesDisplay.getTabByName('Data');
            drawnow
            
            
            this.hsaInterferogram.build     (uitData, this.hFigure, 50, 80, 850, 770);
            
            this.uitStack.build             (uitData, 430, 17, 120, 40);
            this.uibLeft.build              (uitData, 370, 22, 35, 35);
            this.uibRight.build             (uitData, 570, 22, 35, 35);
            this.uitMetaInfo.build          (uitData, 50, 850, 850, 50);
            
            
            this.uitStack.setFontSize(28);
            this.uitStack.setColor([0, 0, .4]);
            this.uitStack.setAlign('center');
            
            this.uitMetaInfo.setFontSize(11);
            this.uitMetaInfo.setColor([0, .4, .1]);
            this.uitMetaInfo.setAlign('left');
            this.uitMetaInfo.setBackgroundColor([1, 1, 1]);
            
            
            % Axes:Probe and object
            uitPhase = this.uitgAxesDisplay.getTabByName('Probe and object');
            this.haProbeAmp = axes('Parent', uitPhase, ...
                'Units', 'pixels', ...
                'Position', [60, 500, 400, 340], ...
                'XTick', [], 'YTick', []);
            this.haProbePha = axes('Parent', uitPhase, ...
                'Units', 'pixels', ...
                'Position', [520, 500, 400, 340], ...
                'XTick', [], 'YTick', []);
            this.haObjectAmp = axes('Parent', uitPhase, ...
                'Units', 'pixels', ...
                'Position', [60, 50, 400, 340], ...
                'XTick', [], 'YTick', []);
            this.haObjectPha = axes('Parent', uitPhase, ...
                'Units', 'pixels', ...
                'Position', [520, 50, 400, 340], ...
                'XTick', [], 'YTick', []);
            
            % Axes:Initial guess
            uitGuess = this.uitgAxesDisplay.getTabByName('Initial guess');
            this.haGuessProbeAmp = axes('Parent', uitGuess, ...
                'Units', 'pixels', ...
                'Position', [60, 500, 400, 340], ...
                'XTick', [], 'YTick', []);
            this.haGuessProbePha = axes('Parent', uitGuess, ...
                'Units', 'pixels', ...
                'Position', [520, 500, 400, 340], ...
                'XTick', [], 'YTick', []);
            this.haGuessObjectAmp = axes('Parent', uitGuess, ...
                'Units', 'pixels', ...
                'Position', [60, 50, 400, 340], ...
                'XTick', [], 'YTick', []);
            this.haGuessObjectPha = axes('Parent', uitGuess, ...
                'Units', 'pixels', ...
                'Position', [520, 50, 400, 340], ...
                'XTick', [], 'YTick', []);
            
            % Axes:Reconstruction
            uitRecon = this.uitgAxesDisplay.getTabByName('Reconstruction');
            this.haReconProbeAmp = axes('Parent', uitRecon, ...
                'Units', 'pixels', ...
                'Position', [60, 500, 400, 340], ...
                'XTick', [], 'YTick', []);
            this.haReconProbePha = axes('Parent', uitRecon, ...
                'Units', 'pixels', ...
                'Position', [520, 500, 400, 340], ...
                'XTick', [], 'YTick', []);
            this.haReconObjectAmp = axes('Parent', uitRecon, ...
                'Units', 'pixels', ...
                'Position', [60, 50, 400, 340], ...
                'XTick', [], 'YTick', []);
            this.haReconObjectPha = axes('Parent', uitRecon, ...
                'Units', 'pixels', ...
                'Position', [520, 50, 400, 340], ...
                'XTick', [], 'YTick', []);
            
            % Axes:Analysis
            uitAnalysis = this.uitgAxesDisplay.getTabByName('Analysis');
            this.haAnalysis = axes('Parent', uitAnalysis, ...
                'Units', 'pixels', ...
                'Position', [200, 360, 560, 480], ...
                'XTick', [], 'YTick', []);
            this.haZernikeDecomp = axes('Parent', uitAnalysis, ...
                'Units', 'pixels', ...
                'Position', [80, 100, 780, 180], ...
                'XTick', [], 'YTick', []);
            
            this.uitRMS.build           (uitAnalysis, 80, 880, 200, 20);
            
            % Axes:Log
            uitLog = this.uitgAxesDisplay.getTabByName('Log');
            this.htLog = uitable(uitLog,'Position',[10 10 930 900],'CellSelectionCallback',@(src, evt)this.cb(src,evt));
            
            %% Controls main:
            this.hpControls = uipanel(...
                'Parent', this.hFigure,...
                'Units', 'pixels',...
                'Title', 'Controls',...
                'FontWeight', 'Bold',...
                'Clipping', 'on',...
                'BorderWidth',1, ...
                'Position', [970 20 620 950] ...
                );
            
            
            %% Controls: Experiment setup:
            this.hpAnalysisSetup = uipanel(...
                'Parent', this.hpControls,...
                'Units', 'pixels',...
                'Title', 'Experiment setup',...
                'FontWeight', 'Bold',...
                'Clipping', 'on',...
                'BorderWidth',1, ...
                'Position', [10 700 600 235] ...
                );
            drawnow
            Offset0=0;
            this.uieLambda.build    (this.hpAnalysisSetup, 20, 20+Offset0, 90, 20);
            this.uieScanRange.build         (this.hpAnalysisSetup, 125, 20+Offset0, 90, 20);
            this.uieNA.build        (this.hpAnalysisSetup, 20, 60+Offset0, 90, 20);
            this.uiez2.build        (this.hpAnalysisSetup, 125, 60+Offset0, 90, 20);
            this.uieGratTilt.build  (this.hpAnalysisSetup, 20, 100+Offset0, 90, 20);
            this.uieDetTilt.build   (this.hpAnalysisSetup, 125, 100+Offset0, 90, 20);
            this.uieGlobalRot.build (this.hpAnalysisSetup, 20, 140+Offset0, 90, 20);
            this.uieDetSize.build   (this.hpAnalysisSetup, 125, 140+Offset0, 90, 20);
            this.uieCenterObstruction.build   (this.hpAnalysisSetup, 20, 180+Offset0, 90, 20);
            this.uieBinning.build   (this.hpAnalysisSetup, 125, 180+Offset0, 90, 20);
            
            this.prControlsSetup.build(this.hpAnalysisSetup, 235, 15, 350, 200);
            
            %% Controls: Data panel:
            this.hpData = uipanel(...
                'Parent', this.hpControls,...
                'Units', 'pixels',...
                'Title', 'Data',...
                'FontWeight', 'Bold',...
                'Clipping', 'on',...
                'BorderWidth',1, ...
                'Position', [10 390 600 295] ...
                );
            drawnow
            
            Offset1=5;
            this.uitgSelectDataSource.build ...
                (this.hpData, 10, 10 + Offset1, 575, 200);
            
            this.uieCCDCenter.build     (this.hpData, 10, 215+Offset1, 120, 20);
            this.uieObsOffset.build     (this.hpData, 140, 215+Offset1, 120, 20);
            this.uibSelectCenter.build  (this.hpData, 10, 260+Offset1, 80, 20);
            this.uibResetCenter.build (this.hpData, 100, 260+Offset1, 110, 20);
            this.uicbAutoCenter.build  (this.hpData, 220, 260+Offset1, 100, 20);
            this.uicbFourierPtychography.build (this.hpData, 460, 260+Offset1, 140, 20);
            this.uipSelectMask.build    (this.hpData, 340, 215+Offset1, 240, 20)
            this.uibLoadMask.build      (this.hpData, 340, 260+Offset1, 100, 20);
            this.uitgSelectDataSource.selectTabByName('Probe and object');
            drawnow;
            
            % Controls:Data:LI
            Offset2=-30;
            uitLISingle = this.uitgSelectDataSource.getTabByName('Load patterns');
            this.uieLogFileNameSingle.build (uitLISingle, 20, 35+Offset2, 380, 20);
            this.uibSetLogFileSingle.build  (uitLISingle, 420, 35 + 13+Offset2, 100, 20);
            
            this.uibLoadSingleFromLog.build (uitLISingle, 420, 130+Offset2, 100, 20);
            this.uibLoadSingleFromImg.build (uitLISingle, 420, 100+Offset2, 100, 20);
            
            this.uilSingleList.build(uitLISingle, 20, 70+Offset2, 380, 100);
            
            % Controls:Data:LIstack
            Offset3=-30;
            uitLIStack = this.uitgSelectDataSource.getTabByName('Load P/S Series');
            this.uieLogFileNameStack.build  (uitLIStack, 20, 35+Offset3, 380, 20);
            this.uibSetLogFileStack.build   (uitLIStack, 420, 35 + 13+Offset3, 100, 20);
            
            this.uipDataType.build          (uitLIStack, 420, 70+Offset3, 100, 20);
            
            this.uibLoadStackFromLog.build (uitLIStack, 420, 130+Offset3, 120, 20);
            
            this.uilStackList.build(uitLIStack, 20, 70+Offset3, 380, 100);
            
            
            % Controls: Data panel:Simulation
            Offset4=-40;
            uitSim = this.uitgSelectDataSource.getTabByName('Simulation');drawnow
            this.uieRes.build       (uitSim, 10, 50+Offset4, 50, 20);
            this.uiez1.build        (uitSim, 70, 50+Offset4, 50, 20);
            this.uieNp.build        (uitSim, 130, 50+Offset4, 50, 20);
            this.uieScanSteps.build  (uitSim, 190, 50+Offset4, 50, 20);
            this.uiePhaseStepsSim.build  ...
                (uitSim, 250, 50+Offset4, 170, 20);
            this.uibLoadScanningPos.build ...
                (uitSim, 435, 50 + 13+Offset4, 100, 20);
            this.uieZrn.build       (uitSim, 10, 90+Offset4, 410, 20);
            
            this.uibLoadZrn.build   (uitSim, 435, 103+Offset4, 100, 20);
            this.uibSimulate.build  (uitSim, 150, 140+Offset4, 100, 20);
            this.uibSimulateS.build (uitSim, 275, 140+Offset4, 100, 20);
            
            this.uipbExposureProgress.build ...
                (uitSim, 10, 170+Offset4, 200, 20);
            % Probe and object
            Offsetp=0;
            uitProbe = this.uitgSelectDataSource.getTabByName('Probe and object');drawnow
            this.uipProbeType.build    (uitProbe, 10, 10+Offsetp, 100, 20);
            this.uipObjectType.build    (uitProbe, 250, 10+Offsetp, 100, 20);
            this.uieRprobe.build    (uitProbe, 120, 10+Offsetp, 100, 20);
            this.uipPropagator.build    (uitProbe, 10, 120+Offsetp, 140, 20);
            this.uibLoadProbe.build    (uitProbe, 10, 90+Offsetp, 100, 20);
            this.uibLoadObject.build    (uitProbe, 250, 90+Offsetp, 100, 20);
            this.uibGenProbeObject.build    (uitProbe, 250, 130+Offsetp, 100, 20);
            this.uibCopyProbe.build    (uitProbe, 120, 90+Offsetp, 100, 20);
            this.uicbGuess.build (uitProbe, 430, 10+Offsetp, 100, 20);
            this.uitOverlap.build (uitProbe, 430, 90+Offsetp, 120, 20);
            this.uitSampling.build (uitProbe, 430, 50+Offsetp, 120, 20);
            this.uieProbeOffset.build     (uitProbe, 10, 50+Offsetp, 100, 20);
            this.uieProbeAmp.build     (uitProbe, 120, 50+Offsetp, 100, 20);
            this.uibSimulatePO.build   (uitProbe, 430, 130+Offsetp, 100, 20);
            
            % FPM
            uitFPM = this.uitgSelectDataSource.getTabByName('FPM');drawnow
            this.uieNAo.build          (uitFPM, 10, 10, 90, 20);
            this.uieLo.build          (uitFPM, 10, 50, 90, 20);
            this.uieMag.build          (uitFPM, 10, 90, 90, 20);
            this.uieScanAngles.build  (uitFPM, 10, 130, 170, 20);
            
            % Sim segment
            uitSimSegment = this.uitgSelectDataSource.getTabByName('Sim segments');drawnow
            this.uieSegmentPath.build  (uitSimSegment, 20, 25, 380, 20);
            this.uibLoadSegment.build   (uitSimSegment, 420, 38, 100, 20);
            
            % Custom Sim
            uitCSim = this.uitgSelectDataSource.getTabByName('Sim stochastics');drawnow
            Offset8=27;
            this.uibCustomSim.build             (uitCSim, 410+Offset8, 180+Offset4, 90, 20);
            this.uibReset.build                 (uitCSim, 280+Offset8, 180+Offset4, 90, 20);
            this.uieXLinearDrift.build          (uitCSim, 20+Offset8, 50, 90, 25);
            this.uieYLinearDrift.build          (uitCSim, 20+Offset8, 90, 90, 25);
            this.uieZLinearDrift.build          (uitCSim, 20+Offset8, 10, 90, 25);
            this.uiePhaseShiftingError.build    (uitCSim, 150+Offset8, 10, 90, 25);
            this.uieGratingTiltErr.build        (uitCSim, 150+Offset8, 50, 90, 25);
            this.uieDetectorTiltErr.build       (uitCSim, 150+Offset8, 90, 90, 25);
            this.uieShotToShot.build            (uitCSim, 280+Offset8, 10, 90, 25);
            this.uie2ndOrderStrength.build      (uitCSim, 280+Offset8, 50, 90, 25);
            this.uie11orderStrength.build       (uitCSim, 280+Offset8, 90, 90, 25);
            this.uieDetCurve.build              (uitCSim, 410+Offset8, 10, 90, 25);
            this.uieFlareLevel.build            (uitCSim, 410+Offset8, 50, 90, 25);
            this.uieMSFR.build                  (uitCSim, 410+Offset8, 90, 90, 25);
            this.uieNonlinearity.build          (uitCSim, 20+Offset8, 130, 90, 25);
            this.uieAirflow.build               (uitCSim, 150+Offset8, 130, 90, 25);
            
            %% Controls: Reconstruction
            this.hpPhase = uipanel(...
                'Parent', this.hpControls,...
                'Units', 'pixels',...
                'Title', 'Reconstruction',...
                'FontWeight', 'Bold',...
                'Clipping', 'on',...
                'BorderWidth',1, ...
                'Position', [10 190 600 190] ...
                );
            drawnow
            %             Offset5=20;
            this.uitgAnalysisDomain.build(this.hpPhase, 10, 20, 400, 160);drawnow
            %
            uitRPIE = this.uitgAnalysisDomain.getTabByName('rPIE');
            uitRAAR = this.uitgAnalysisDomain.getTabByName('RAAR');
            uitWDD = this.uitgAnalysisDomain.getTabByName('WDD');
            uitML = this.uitgAnalysisDomain.getTabByName('ML');
            this.uibComputePhase.build  (this.hpPhase, 415, 160, 160, 20);
            this.uibStop.build  (this.hpPhase, 415, 130, 160, 20);
            this.uitIteration.build           (this.hpPhase, 180, 20, 250, 20);
            this.uieAlpha.build           (this.hpPhase, 420, 90, 70, 20);
            this.uieBeta.build           (this.hpPhase, 500, 90, 70, 20);
            this.uieMaxIteration.build           (this.hpPhase, 420, 50, 70, 20);
            this.uieAccuracy.build           (this.hpPhase, 500, 50, 70, 20);
            this.uicbGPU.build           (this.hpPhase, 420, 16, 120, 20);
            
            % Control: Reconstruction: rPIE
            this.uieGamma.build           (uitRPIE, 20, 10, 70, 20);
            this.uicbCorrectPos.build     (uitRPIE, 20, 50, 100, 20);
            this.uipCorrectMethod.build (uitRPIE, 20, 80, 150, 20);
            this.uieModeNumber.build  (uitRPIE,100,10,70,20);
            this.uilSelectMode.build (uitRPIE,200,10,80,95);
            this.uibShowMode.build (uitRPIE,300,100,80,20);
            
            
            % Control: Reconstruction: RAAR
            this.uieDelta.build           (uitRAAR, 20, 20, 100, 20);
            
            % Control: Reconstruction: WDD
            this.uieEpsilon.build           (uitWDD, 20, 20, 100, 20);
            
            % Control: Reconstruction: ML
            this.uipLikelihoodType.build    (uitML, 200, 10, 100, 20);
            this.uieRegularization.build           (uitML, 20, 10, 70, 20);
            
            %% Controls: Analysis
            this.hpAnalysis = uipanel(...
                'Parent', this.hpControls,...
                'Units', 'pixels',...
                'Title', 'Analysis',...
                'FontWeight', 'Bold',...
                'Clipping', 'on',...
                'BorderWidth', 1, ...
                'Position', [10 10 600 170] ...
                );
            drawnow
            
            this.uipSelectObject.build     (this.hpAnalysis, 20, 20, 180, 20);
            this.uipSelectRegion.build       (this.hpAnalysis, 20, 70, 180, 20);
            this.uibAnalyze.build   (this.hpAnalysis, 415, 140, 160, 20);
            this.uibLoadAnalysisRegion.build   (this.hpAnalysis, 20, 120, 120, 20);
            this.uibDataCursor.build   (this.hpAnalysis, 415, 20, 160, 20);
            this.uieSigma.build   (this.hpAnalysis, 220, 20, 80, 20);
            
            % Set hsa offset:
            this.hsaInterferogram.setAxesOffset([dTgPx + uitData.Position(1), dTgPy + uitData.Position(2)]);
            
            % refresh button state:
            this.setState(this.u8State);
        end
        
    end
    
end

