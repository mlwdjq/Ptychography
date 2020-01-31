classdef PIE_Analyze < mic.Base
    
    
    properties (Constant)
        dWidth  = 1600;
        dHeight =  1000;
        
        % Axes tab IDS
        U8DATA              = 1
        U8PROBEOBJECT       = 2
        U8GUESS             = 3
        U8RECONSTRUCTION    = 4
        U8LOG    = 5
        
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
        
        U8MAXSTATES         = 5
        
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
        
        % Main interferogram(s)
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
        dReconstructed
        dZernike
        dZernikeResidual
        dWx
        dWy
        dWxNoTilt
        dWyNoTilt
        dWxUnwrapped
        dWyUnwrapped
        dResult
        dProbe
        dObject
        dProbeGuess
        dObjectGuess
        
        % Shear guides
        dBeamWidthEstPx = 1
        dShearPix = 1
        dObstructionWidthEstPx
        dObstructionShearPix
        % Axes tab group
        uitgAxesDisplay     % displays axes: Interferogram, Phase, Reconstruction
        
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
        haReconstructed
        haZernikeDecomp
        uicbRemoveSphere
        uitxWFE
        uicbRemoveTiltX
        uicbRemoveTiltY
        uicbRemoveDefocus
        uicbScaledNW
        uitRMS
        uitRMSFit
        % Axes: Log
        htLog
        
        % Controls
        hpControls
        hpAnalysisSetup
        hpLoadInterferogram
        hpPhase
        hpReconstruction
        
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
        
        uieCCDCenter
        uieObsOffset
        uibSelectCenter
        uibResetCenter
        
        uipSelectMask
        uibLoadMask
        
        % Controls:Data:Load interferogram single/stack
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
        uibLoadPhaseStepsSim
        uiez1
        uieNp % number of photos
        uibLoadZrn
        uibSimulate
        uibSimulateS
        uieScanSteps
        uipbExposureProgress
        uicbMultiPhaseShifting
        
        % Controls:Data:Probe and Object
        uipProbeType
        uipObjectType
        uipPropagator
        uibLoadProbe
        uibLoadObject
        uibGenProbeObject
        uieRprobe
        uicbGuess
        
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
        uieMSFN
        uieNonlinearity
        uieAirflow
        
        % Controls:Phase
        uitgAnalysisDomain
        uipUnwrapEngine
        uibComputePhase
        uieFDZ1
        uitGramRot
        % Controls:Phase:FD
        
        uieFilterWidth
        uipWindowFunction
        
        % Controls:Phase:TD
        uiePhaseStepsTD
        uipFTType
        uieLowPass
        uieTDZ1
        uibTDLoadPhaseStepsTD
        
        % Controls:Reconstruction
        uitgReconstructionType
        uieNZernikes
        uibReconstruct
        uipReconUnit
        uicbResidualError
        
        
        % ...:Reconstruction:Rimmer
        uipRimmerType
        uicbAutoLoadRim
        % ...:Reconstruction:Derivative Basis
        
        uicbAutoLoadDB
        uieNZernikesBasis
        
        % ...:Reconstruction:Fourier
        uipFittingType
        uicbAutoLoadFT
        uicbOrthogonalization
        
        
        
        
        uitgLoadType
        uiPSSimulationRecaller
        
        
        hAx % cell array of axes
        
        
    end
    
    properties (SetAccess = private)
        
    end
    
    methods
        function this = PIE_Analyze()
            this.init()
        end
        
        
        
        function init(this)
            
            this.uitgAxesDisplay = ...
                mic.ui.common.Tabgroup('ceTabNames', {'Data', 'Probe and object','Initial guess', 'Reconstruction','Log'});
            
            % Axes:
            this.hsaInterferogram = mic.ui.axes.ScalableAxes(...
                'fhOnDomainChange', @(cDomain)this.replot(this.U8DATA,  cDomain));
            this.uibLeft        = mic.ui.common.Button('cText', '<', 'fhDirectCallback', @(src, evt)this.cb(src));
            this.uibRight       = mic.ui.common.Button('cText', '>', 'fhDirectCallback', @(src, evt)this.cb(src));
            this.uitStack       = mic.ui.common.Text('cVal', '0/0');
            this.uitMetaInfo    = mic.ui.common.Text('cVal', '');
            
            this.uitRMS                 = mic.ui.common.Text('cVal','RMS:');
            this.uitRMSFit                 = mic.ui.common.Text('cVal','RMSFit:');
            
            %this.uitRMS.cVal([]);
            
            % Controls: Experiment setup panel
            this.uieLambda      = mic.ui.common.Edit('cLabel', 'Lambda (nm)', 'cType', 'd');
            this.uieScanRange   = mic.ui.common.Edit('cLabel', 'Scan range (mm)', 'cType', 'd');
            this.uieNA          = mic.ui.common.Edit('cLabel', 'NA', 'cType', 'd');
            this.uiez2          = mic.ui.common.Edit('cLabel', 'Z_2 (mm)', 'cType', 'd');
            this.uieGratTilt    = mic.ui.common.Edit('cLabel', 'Gr. Tilt (deg)', 'cType', 'd');
            this.uieDetTilt     = mic.ui.common.Edit('cLabel', 'Dt. Tilt (deg)', 'cType', 'd');
            this.uieGlobalRot   = mic.ui.common.Edit('cLabel', 'CCD Rot (deg)', 'cType', 'd');
            this.uieDetSize     = mic.ui.common.Edit('cLabel', 'Dt. Size (mm)', 'cType', 'd');
            this.uieCenterObstruction     = mic.ui.common.Edit('cLabel', 'Center Obstruction', 'cType', 'd');
            this.uieBinning     = mic.ui.common.Edit('cLabel', 'Binning', 'cType', 'd');
            
            
            this.uieLambda.set(13.5);
            this.uieScanRange.set(0.1);
            this.uieNA.set(0.5);
            this.uiez2.set(20);
            this.uieGratTilt.set(1.12);
            this.uieDetTilt.set(1.12);
            this.uieGlobalRot.set(0);
            this.uieDetSize.set(26);
            this.uieCenterObstruction.set(0);
            this.uieBinning.set(512);
            
            this.prControlsSetup = mic.ui.common.PositionRecaller(...
                'cConfigPath', fullfile(this.cAppPath, '+config'), ...
                'cName', 'config', ...
                'hGetCallback', @()this.prGetters(this.prControlsSetup), ...
                'hSetCallback', @(dRecall)this.prSetters(this.prControlsSetup, dRecall));
            
            % Controls: Data panel
            this.uitgSelectDataSource = ...
                mic.ui.common.Tabgroup('ceTabNames', {'Load Interferogram', 'Load P/S Series', 'Simulation','Probe and object', 'Sim stochastics'});
            
            this.uieCCDCenter       = mic.ui.common.Edit('cLabel', 'CCD Center pixel', 'cType', 'c', ...
                'fhDirectCallback', @(src, evt)this.cb(src), 'lNotifyOnProgrammaticSet', false);
            this.uicbAutoCenter        = mic.ui.common.Checkbox('cLabel', 'Auto Center',  'fhDirectCallback', @(src, evt)this.cb(src));
            this.uieObsOffset       = mic.ui.common.Edit('cLabel', 'Obstruction offset pixel', 'cType', 'c', ...
                'fhDirectCallback', @(src, evt)this.cb(src), 'lNotifyOnProgrammaticSet', false);
            this.uibSelectCenter    = mic.ui.common.Button('cText', 'Select center', 'fhDirectCallback', @(src, evt)this.cb(src));
            this.uibResetCenter     = mic.ui.common.Button('cText', 'Set default center', 'fhDirectCallback', @(src, evt)this.cb(src));
            
            this.uipSelectMask      = mic.ui.common.Popup('cLabel', 'Select mask', 'ceOptions', {'Compute from det. geom.','Compute from MET5','Compute from eliptical geom.'}, ...
                'fhDirectCallback',@(src, evt)this.cb(src), 'lShowLabel', true);
            
            this.uibLoadMask        = mic.ui.common.Button('cText', 'Load mask', 'fhDirectCallback', @(src, evt)this.cb(src));
            
            this.uieCCDCenter.set('[]');
            this.uieObsOffset.set('[]');
            this.uicbAutoCenter.set(false);
            
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
            this.uieRes         = mic.ui.common.Edit('cLabel', 'Res', 'cType', 'd');
            this.uiez1          = mic.ui.common.Edit('cLabel', 'Z_1 (mm)', 'cType', 'd');
            this.uieNp          = mic.ui.common.Edit('cLabel', 'N Phtn', 'cType', 'd');
            this.uieZrn         = mic.ui.common.Edit('cLabel', 'Zernike couples vector [N X 2]', 'cType', 'c','fhDirectCallback', @(src, evt)this.cb(src), 'lNotifyOnProgrammaticSet', false);
            
            this.uicbMultiPhaseShifting = mic.ui.common.Checkbox('cLabel', 'Multi-shifting',  'fhDirectCallback', @(src, evt)this.cb(src));
            this.uieScanSteps    = mic.ui.common.Edit('cLabel', 'Scanning steps', 'cType', 'd', 'fhDirectCallback', @(src, evt)this.cb(src), 'lNotifyOnProgrammaticSet', false);
            this.uiePhaseStepsSim  = mic.ui.common.Edit('cLabel', 'Phase steps', 'cType', 'c', 'fhDirectCallback', @(src, evt)this.cb(src), 'lNotifyOnProgrammaticSet', false);
            this.uibLoadZrn     = mic.ui.common.Button('cText', 'Load Zrn File', 'fhDirectCallback', @(src, evt)this.cb(src));
            this.uibSimulate    = mic.ui.common.Button('cText', 'Simulate Single', 'fhDirectCallback', @(src, evt)this.cb(src));
            this.uibSimulateS   = mic.ui.common.Button('cText', 'Simulate Stack', 'fhDirectCallback', @(src, evt)this.cb(src));
            
            this.uibLoadPhaseStepsSim...
                = mic.ui.common.Button('cText', 'Load Phase steps', 'fhDirectCallback', @(src, evt)this.cb(src));
            this.uipbExposureProgress = mic.ui.common.ProgressBar(...
                'dColorFill', [.4, .4, .8], ...
                'dColorBg', [1, 1, 1], ...
                'dHeight', 15, ...
                'dWidth', 455);
            
            this.uiez1.set(-50);
            this.uieNp.set(30000);
            this.uieRes.set(256);
            this.uieScanSteps.set(4);
            this.uiePhaseStepsSim.set('[0:pi/2:3*pi/2;0:pi/2:3*pi/2]''')
            this.uieZrn.set('[]');
            this.uicbMultiPhaseShifting.set(true);
            
            % Controls:Data:Probe and object
            this.uipProbeType     = mic.ui.common.Popup('cLabel', 'Probe type', 'ceOptions', {'Defocus wave','Plane wave'}, ...
                'fhDirectCallback',@(src, evt)this.cb(src), 'lShowLabel', true);
            this.uipObjectType     = mic.ui.common.Popup('cLabel', 'Object type', 'ceOptions', {'Vacuum','Cameraman'}, ...
                'fhDirectCallback',@(src, evt)this.cb(src), 'lShowLabel', true);
            this.uipPropagator     = mic.ui.common.Popup('cLabel', 'Propagator', 'ceOptions', {'angular spectrum','fourier','fresnel'}, ...
                'fhDirectCallback',@(src, evt)this.cb(src), 'lShowLabel', true);
            this.uibLoadProbe    = mic.ui.common.Button('cText', 'Load probe', 'fhDirectCallback', @(src, evt)this.cb(src));
            this.uibLoadObject    = mic.ui.common.Button('cText', 'Load object', 'fhDirectCallback', @(src, evt)this.cb(src));
            this.uibGenProbeObject   = mic.ui.common.Button('cText', 'Generate', 'fhDirectCallback', @(src, evt)this.cb(src));
            this.uieRprobe = mic.ui.common.Edit('cLabel', 'Probe radius on det (mm)', 'cType', 'd', 'fhDirectCallback', @(src, evt)this.cb(src), 'lNotifyOnProgrammaticSet', false);
            this.uieRprobe.set(this.uieDetSize.get()/4); 
            this.uicbGuess = mic.ui.common.Checkbox('cLabel', 'Initial guess',  'fhDirectCallback', @(src, evt)this.cb(src));
            this.uicbGuess.set(true);
            
            % Controls:Data:Sim stochastics
            this.uibCustomSim = mic.ui.common.Button('cText', 'Custom Sim', 'fhDirectCallback', @(src, evt)this.cb(src));
            this.uibReset = mic.ui.common.Button('cText', 'Reset', 'fhDirectCallback', @(src, evt)this.cb(src));
            
            this.uiePhaseShiftingError = mic.ui.common.Edit('cLabel', 'Grat Phase Shift er (nm)', 'cType', 'd', 'fhDirectCallback', @(src, evt)this.cb(src), 'lNotifyOnProgrammaticSet', false);
            
            this.uieZLinearDrift    = mic.ui.common.Edit('cLabel', 'Z Lin Drift (nm)', 'cType', 'd', 'fhDirectCallback', @(src, evt)this.cb(src), 'lNotifyOnProgrammaticSet', false);
            this.uieXLinearDrift   = mic.ui.common.Edit('cLabel', 'X Lin Drift (nm)', 'cType', 'd', 'fhDirectCallback', @(src, evt)this.cb(src), 'lNotifyOnProgrammaticSet', false);
            this.uieYLinearDrift   = mic.ui.common.Edit('cLabel', 'Y Lin Drift (nm)', 'cType', 'd', 'fhDirectCallback', @(src, evt)this.cb(src), 'lNotifyOnProgrammaticSet', false);
            this.uieGratingTiltErr  = mic.ui.common.Edit('cLabel', 'Grat tilt error(deg)', 'cType', 'd', 'fhDirectCallback', @(src, evt)this.cb(src), 'lNotifyOnProgrammaticSet', false);
            this.uieDetectorTiltErr = mic.ui.common.Edit('cLabel', 'Det tilt error (deg)', 'cType', 'd', 'fhDirectCallback', @(src, evt)this.cb(src), 'lNotifyOnProgrammaticSet', false);
            this.uieShotToShot       = mic.ui.common.Edit('cLabel', 'Shot to shot (%)', 'cType', 'd', 'fhDirectCallback', @(src, evt)this.cb(src), 'lNotifyOnProgrammaticSet', false);
            this.uie2ndOrderStrength = mic.ui.common.Edit('cLabel', '2nd order strength)', 'cType', 'd', 'fhDirectCallback', @(src, evt)this.cb(src), 'lNotifyOnProgrammaticSet', false);
            this.uie11orderStrength = mic.ui.common.Edit('cLabel', '11 order strength)', 'cType', 'd', 'fhDirectCallback', @(src, evt)this.cb(src), 'lNotifyOnProgrammaticSet', false);
            this.uieDetCurve       = mic.ui.common.Edit('cLabel', 'Dt. Curve(um)', 'cType', 'd', 'fhDirectCallback', @(src, evt)this.cb(src), 'lNotifyOnProgrammaticSet', false);
            this.uieFlareLevel       = mic.ui.common.Edit('cLabel', 'DC flare', 'cType', 'd', 'fhDirectCallback', @(src, evt)this.cb(src), 'lNotifyOnProgrammaticSet', false);
            this.uieMSFN       = mic.ui.common.Edit('cLabel', 'MSFN', 'cType', 'd', 'fhDirectCallback', @(src, evt)this.cb(src), 'lNotifyOnProgrammaticSet', false);
            this.uieNonlinearity       = mic.ui.common.Edit('cLabel', 'Nonlinearity (a,b)', 'cType', 'c', 'fhDirectCallback', @(src, evt)this.cb(src), 'lNotifyOnProgrammaticSet', false);
            this.uieAirflow       = mic.ui.common.Edit('cLabel', 'Airflow', 'cType', 'd', 'fhDirectCallback', @(src, evt)this.cb(src), 'lNotifyOnProgrammaticSet', false);
            this.uieNonlinearity.set('[0, 0]');
            % Controls: Phase
            this.uitgAnalysisDomain = mic.ui.common.Tabgroup('ceTabNames', {'Time Domain', 'Fourier Domain'});
            
            this.uipUnwrapEngine        = mic.ui.common.Popup('cLabel', 'Unwrapping algorithm', 'ceOptions', {'Sorting reliability unwrap'}, ...
                'fhDirectCallback',@(src, evt)this.cb(src), 'lShowLabel', true);
            this.uibComputePhase        = mic.ui.common.Button('cText', 'Process Phase',  'fhDirectCallback', @(src, evt)this.cb(src));
            
            % Controls:Phase:TD
            this.uiePhaseStepsTD        = mic.ui.common.Edit('cLabel', 'Phase steps [N X 2]', 'cType', 'c', 'fhDirectCallback', @(src, evt)this.cb(src), 'lNotifyOnProgrammaticSet', false);
            this.uieLowPass             = mic.ui.common.Edit('cLabel', 'Low pass (N*Ws,0=none)', 'cType', 'd', 'fhDirectCallback', @(src, evt)this.cb(src));
            
            this.uipFTType              = mic.ui.common.Popup('cLabel', 'Fourier Transform Type', 'ceOptions', {'FFT (uniform steps)', 'DFT (custom steps)'}, ...
                'fhDirectCallback',@(src, evt)this.cb(src), 'lShowLabel', true);
            this.uibTDLoadPhaseStepsTD  = mic.ui.common.Button('cText', 'Load Phase steps', 'fhDirectCallback', @(src, evt)this.cb(src));
            this.uiePhaseStepsTD.set('[]');
            
            this.uieFDZ1                = mic.ui.common.Edit('cLabel', 'Z1', 'cType', 'd', 'fhDirectCallback', @(src, evt)this.cb(src));
            this.uitGramRot          = mic.ui.common.Text('cVal','Gram Rot:0degree');
            this.uieFDZ1.set(1);
            % Controls:Phase:FD
            this.uipWindowFunction      = mic.ui.common.Popup('cLabel', 'Filter type', 'ceOptions', {'Gauss', 'Hanning'}, ...
                'fhDirectCallback',@(src, evt)this.cb(src), 'lShowLabel', true);
            this.uieFilterWidth          = mic.ui.common.Edit('cLabel', 'Filter width', 'cType', 'd', 'fhDirectCallback', @(src, evt)this.cb(src));
            this.uieFilterWidth.set(0.1);
            
            % Controls:Reconstruction
            this.uitgReconstructionType = mic.ui.common.Tabgroup('ceTabNames', {'Rimmer', 'Derivative Basis','Fourier'});
            
            this.uipRimmerType          = mic.ui.common.Popup('cLabel', 'Rimmer type', 'ceOptions', {'2X downsample', 'Double shear'}, ...
                'fhDirectCallback',@(src, evt)this.cb(src), 'lShowLabel', true);
            
            this.uipFittingType         = mic.ui.common.Popup('cLabel', 'Fitting Type', 'ceOptions', {'Reconstructed wavefront', 'Shearing wavefronts'}, ...
                'fhDirectCallback',@(src, evt)this.cb(src), 'lShowLabel', true);
            
            this.uicbAutoLoadRim        = mic.ui.common.Checkbox('cLabel', 'Auto save/load matrix',  'fhDirectCallback', @(src, evt)this.cb(src));
            this.uicbAutoLoadDB         = mic.ui.common.Checkbox('cLabel', 'Auto save/load basis',  'fhDirectCallback', @(src, evt)this.cb(src));
            this.uicbAutoLoadFT         = mic.ui.common.Checkbox('cLabel', 'Auto save/load basis',  'fhDirectCallback', @(src, evt)this.cb(src));
            this.uicbOrthogonalization         = mic.ui.common.Checkbox('cLabel', 'Orthogonalization',  'fhDirectCallback', @(src, evt)this.cb(src));
            this.uieNZernikes           = mic.ui.common.Edit('cLabel', 'N Zernikes', 'cType', 'd', 'fhDirectCallback', @(src, evt)this.cb(src));
            this.uipReconUnit        = mic.ui.common.Popup('cLabel', 'Unit', 'ceOptions', {'Wave','mWave','nm'}, ...
                'fhDirectCallback',@(src, evt)this.cb(src), 'lShowLabel', true);
            this.uibReconstruct         = mic.ui.common.Button('cText', 'Reconstruct', 'fhDirectCallback', @(src, evt)this.cb(src));
            this.uieNZernikesBasis      = mic.ui.common.Edit('cLabel', 'Zernike basis order', 'cType', 'd', 'fhDirectCallback', @(src, evt)this.cb(src));
            this.uicbRemoveTiltX        = mic.ui.common.Checkbox('cLabel', 'Remove X tilt',  'fhDirectCallback', @(src, evt)this.cb(src));
            this.uicbRemoveTiltY        = mic.ui.common.Checkbox('cLabel', 'Remove Y tilt',  'fhDirectCallback', @(src, evt)this.cb(src));
            this.uicbRemoveDefocus      = mic.ui.common.Checkbox('cLabel', 'Remove defocus',  'fhDirectCallback', @(src, evt)this.cb(src));
            this.uicbResidualError      = mic.ui.common.Checkbox('cLabel', 'Residual error',  'fhDirectCallback', @(src, evt)this.cb(src));
            this.uicbScaledNW     = mic.ui.common.Checkbox('cLabel', 'Scaled null wavefront',  'fhDirectCallback', @(src, evt)this.cb(src));
            this.uicbRemoveTiltX.set(true);
            this.uicbRemoveTiltY.set(true);
            this.uicbRemoveDefocus.set(true);
            this.uicbResidualError.set(false);
            this.uicbScaledNW.set(true);
            this.uieNZernikesBasis.set(24);
            this.uieNZernikes.set(24);
            this.uicbAutoLoadRim.set(true);
            this.uicbAutoLoadDB.set(true);
            this.uicbAutoLoadFT.set(true);
            this.uicbOrthogonalization.set(false);
            
            % Look for most recent csv:
            
            
            % Add buttons to state flag lists
            this.ceUIBStateControlled = {};
            
            this.ceUIBStateControlled{this.U8STATE_INITIAL} = ...
                {this.uibSimulate, this.uibSimulateS, this.uibLoadSingleFromLog, ...
                this.uibLoadSingleFromImg, this.uibLoadStackFromLog};
            this.ceUIBStateControlled{this.U8STATE_DATA_LOADED} = ...
                {this.uibComputePhase};
            this.ceUIBStateControlled{this.U8STATE_PHASE_PROCESSED} = ...
                {this.uibReconstruct};
            this.ceUIBStateControlled{this.U8STATE_RECONSTRUCTED} = ...
                {};
            
            % Set state:
            this.setState(this.U8STATE_INITIAL);
        end
        
        % Callback handler
        function cb(this, src,evt)
            switch src
                
                case this.uibCustomSim
                    %                     this.customSimulation();
                    this.simStackAndReconstruct();
                    
                case this.uibReset
                    this.setSimParams('zDrift',0,'xDrift',0,'yDrift',0,'gratTilt',0,'detTilt',0,'shotToShot',0,...
                        '2ndOrderStrength',0,'11OrderStrength',0,'xyStageError',0,'detectorCurve',0,'flareLevel',0,...
                        'MSFN',0,'nonlinearity', '[0, 0]','airflow',0);
                    
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
                    
                    % redraw guide lines
                    this.replot(this.U8DATA, []);
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
                    
                    
                case {this.uieZrn, this.uiePhaseStepsTD, this.uieNonlinearity}
                    this.validateCouplesEditBox(src, '[]');
                    
                case this.uieCCDCenter
                    this.validateCouplesEditBox(src, '[]');
                    % redraw guide lines
                    this.replot(this.U8DATA, []);
                    % Need to recenter data
                    this.handleLoadData();
                case this.uibLoadProbe
                    
                case this.uibLoadObject
                    
                case this.uibGenProbeObject
                    probeType = this.uipProbeType.getOptions{this.uipProbeType.getSelectedIndex()};
                    objectType = this.uipObjectType.getOptions{this.uipObjectType.getSelectedIndex()};
                    propagator = this.uipPropagator.getOptions{this.uipPropagator.getSelectedIndex()};
                    initialGuess = this.uicbGuess.get();
                    N           = this.uieRes.get();
                    NA           = this.uieNA.get();
                    Rc_um   = this.uieRprobe.get()*1000;
                    lambda_um   = this.uieLambda.get()/1000;
                    df_um       = this.uiez1.get()*1000;% negative sign corresponds convergent
                    z_um       = this.uiez2.get()*1000;
                    scanRange_um  = this.uieScanRange.get()*1000;
                    scanSteps  = this.uieScanSteps.get();
                    detSize_um  = this.uieDetSize.get()*1000; 
                    dc_um       = detSize_um/N; % detector pixel pitch
                    samplingFactor_det = lambda_um.*z_um/(dc_um*N*dc_um);
                    if samplingFactor_det>1
                        fprintf('Please adjust configurations for propagation sampling\n');
                    end
                    if strcmp(propagator,'fourier')
                        do_um = lambda_um*z_um/N/dc_um;
                    else
                        do_um = dc_um; % object pixel pitch
                    end
                 %% initial probe
                    if NA >0
                            samplingFactor_obj = lambda_um.*z_um/(N*dc_um*do_um);
                            Rc_um = (z_um+df_um)*tan(asin(NA));
                            Rprobe_um = abs(df_um)*tan(asin(NA));
                            [n1,n2]=meshgrid(1:N);
                            n1 = n1-N/2-1;
                            n2 = n2-N/2-1;
                            probe = ifftshift(ifft2(ifftshift(pinhole(round(Rc_um/dc_um),N,N).*...
                                exp(-1i*pi*df_um*dc_um^2/lambda_um/z_um^2*(n1.^2+n2.^2)))));
                    else
                            Rprobe_um = Rc_um;
                            samplingFactor_obj = lambda_um.*abs(df_um)/(N*do_um*do_um);
                            probe = PIE.utils.Propagate (pinhole(round(2*Rprobe_um/do_um),N,N),'angular spectrum',...
                                do_um,lambda_um,abs(df_um));
                    end
                    if initialGuess
                        this.dProbeGuess = single(probe);
                    else
                        this.dProbe = single(probe);
                    end
                    if samplingFactor_obj>1
                        fprintf('Please adjust configurations for propagation sampling\n');
                    end
                    overlap = PIE.utils.overlapRatio(Rprobe_um,scanRange_um/(scanSteps-1)); % overlap ratio of two circles
                    if overlap<0.6 % check overlap
                        fprintf('Overlap ratio(%0.1f%%) is less than 60%%. Please adjust configurations\n',overlap*100);
                    end
                    %% initial object
                    K = round(scanRange_um/do_um)+N;
                    L = round(scanRange_um/do_um)+N; % size of object [K,L]
                    switch objectType
                        case 'Vacuum'
                            object = single(ones(K,L));
                        case 'Cameraman'
                            I = single(flipud(imread('cameraman.tif')));
                            [m,n] = meshgrid(linspace(0,1,L),linspace(0,1,K));
                            [sr,sc,~]= size(I);
                            [p,q] = meshgrid(linspace(0,1,sc),linspace(0,1,sr));
                            object_amp = interp2(p,q,I(:,:,1),m,n);
                            object_amp = mat2gray(object_amp)*0.8+0.2; % object amplitude
                             I =single(imread('pears.png'));
                             [sr,sc,~]= size(I);
                            [p,q] = meshgrid(linspace(0,1,sc),linspace(0,1,sr));
                            object_pha = interp2(p,q,I(:,:,1),m,n);
                            object_pha=mat2gray(object_pha);
                            object_pha = (object_pha-0.5)*2*pi; % object phase
                            object = object_amp.*exp(1i*object_pha);
                            
                    end
                    if initialGuess
                        this.dObjectGuess = single(object);
                    else
                        this.dObject = single(object);
                    end
                    % Make phase tab active:
                    this.uitgAxesDisplay.selectTabByIndex(this.U8PROBEOBJECT);
                    
                    
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
                    
                    
                    
                case this.uieObsOffset
                    this.validateCouplesEditBox(src, '[]');
                    % redraw guide lines
                    this.replot(this.U8DATA, []);
                    % Need to recenter data
                    this.handleLoadData();
                    
                case {this.uiez2, this.uieNA, this.uieLambda, this.uieScanRange, this.uieCenterObstruction}
                    % redraw guide lines
                    this.replot(this.U8DATA, []);
                    
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
                    
                case this.uibLoadSingleFromLog
                    u8loadIdx = this.uilSingleList.getSelectedIndexes();
                    load(this.ceScanPaths{u8loadIdx});
                    this.handleLoadData({dImg}, {stLog});
                    
                    
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
                    
                case this.uibLoadSingleFromImg
                    cDataDir = fullfile(this.cAppPath, '..', '..', '..', 'Data');
                    [fn,pn]=uigetfile({'*.mat;*.SPE','LSI int bundle (*.mat) or WinView (*.SPE)'},'Loading',cDataDir);
                    fileformat=fn(end-2:end);
                    filename= strcat(pn,fn);
                    switch fileformat
                        case 'SPE'
                            ceInts=lsianalyze.utils.speread(filename);
                            this.handleLoadData(ceInts, {});
                            
                        case 'mat'
                            load(filename);
                            try
                                this.handleLoadData({dImg}, {stLog});
                            catch
                                this.handleLoadData(ceInts, {});
                            end
                    end
                    
                    
                case {this.uibLeft, this.uibRight}
                    if (src == this.uibLeft)
                        this.cycleInterferogram(-1);
                    else
                        this.cycleInterferogram(1);
                    end
                case this.htLog
                    s=evt.Indices;
                    [~,N]=size(this.ceAnalysisTable);
                    str='[';
                    for i=1:N-4
                        str=strcat(str,num2str(i),',',num2str(this.ceAnalysisTable{s(1),i+4}),';');
                    end
                    str=[str(1:end-1),']'];
                    this.uieZrn.set(str);
                case this.uibLoadZrn
                    [p, d] = uigetfile();
                    zrn = load([d,p]);
                    this.setZrnString(zrn.zrn);
                    
                case this.uibSimulateS
                    this.simulateInteferograms(true);
                case this.uibSimulate
                    this.simulateInteferograms(false);
                    
                    % Load phase steps for simulation
                case this.uibLoadPhaseStepsSim
                    
                    % Load phase steps for TD analysis
                case this.uibTDLoadPhaseStepsTD
                    
                    
                    
                case this.uieScanSteps
                    dN = this.uieScanSteps.get();
                    dPhseString = sprintf('0:2*pi/%d:%d*pi/%d', dN, dN*2 - 2, dN);
                    this.uiePhaseStepsSim.set(sprintf('[%s;%s]''', dPhseString, dPhseString));
                    
                case this.uiePhaseStepsSim
                    [lValid, vals] = this.validateCouplesEditBox(src, '[]');
                    if lValid
                        [sr, sc] = size(vals);
                        this.uieScanSteps.set(sr);
                    end
                    
                case this.uibComputePhase
                    this.processPhase(this.uitgAnalysisDomain.getSelectedTabName());
                    
                case this.uibReconstruct
                    % this.reconstruct(this.uitgReconstructionType.getSelectedTabName());
                    % call reconstruct functions here
                    this.reconstruct(this.uitgReconstructionType.getSelectedTabName());
                    
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
        
        function processPhase(this, cAnalysisDomain)
            %             for i=1:4
            %                 for j=1:4
            %                     [I,map]=gray2ind(mat2gray(this.ceInt{i,j}),256);
            %                     if(i==1&&j==1)
            %                         imwrite(I,map,'movefig.gif','DelayTime',0.5,'LoopCount',Inf)
            %                     else
            %                         imwrite(I,map,'movefig.gif','WriteMode','append','DelayTime',0.5)
            %                     end
            %                 end
            %             end
            
            if isempty(this.ceInt)
                msgbox('Please load/simulate interferograms first!', 'Error');
                return;
            end
            u8UnwrapAlgorithm               = this.uipUnwrapEngine.getSelectedIndex();
            FilterWidth                     = this.uieFilterWidth.get();
            FilterType                      = this.uipWindowFunction.getSelectedIndex();
            detSize                         = this.uieDetSize.get();
            z2_mm       = this.uiez2.get();
            NA          = this.uieNA.get();
            switch cAnalysisDomain
                case 'Time Domain'
                    if length(this.ceInt)==1
                        msgbox('Please select Fourier Domain!', 'Error');
                        return;
                    end
                    dPhaseSteps             = eval(this.uiePhaseStepsTD.get());
                    u8FourierTransformType  = this.uipFTType.getSelectedIndex();
                    lambda_um   = this.uieLambda.get()/1000;
                    T_um        = this.uieScanRange.get();
                    shearingPercentage=lambda_um/T_um/tan(asin(NA));
                    cutFrequency=1/shearingPercentage;
                    Ws=this.uieLowPass.get();
                    dLowPass= round(cutFrequency*Ws);
                    
                    %                     [this.dWx, this.dWy, this.dWxNoTilt, this.dWyNoTilt, ...
                    %                         this.dWxUnwrapped, this.dWyUnwrapped,z1,CCDrot]  =  ...
                    %                         ...
                    %                             lsianalyze.utils.extractPhaseTDs(...
                    %                                 this.ceInt, dPhaseSteps, u8FourierTransformType, ...
                    %                                 u8UnwrapAlgorithm, this.dAnalysisRegion,this.dAnalysisRegion2, ....
                    %                                 this.uieNA.get(),this.uieScanRange.get(),z2_mm,detSize);
                    switch double(this.uipDataType.getSelectedIndex())
                        case {1, 2} % 1D
                            
                            [this.dWx, this.dWy, this.dWxNoTilt, this.dWyNoTilt, ...
                                this.dWxUnwrapped, this.dWyUnwrapped,z1,CCDrot]  =  ...
                                ...
                                lsianalyze.utils.extract2x1Dphase(...
                                this.ceInt, dPhaseSteps, u8FourierTransformType, ...
                                u8UnwrapAlgorithm, this.dAnalysisRegion,this.dAnalysisRegion2,dLowPass,FilterType, ....
                                NA,this.uieScanRange.get(),z2_mm,detSize);
                            
                            
                        case 3 % 2D
                            
                            [this.dWx, this.dWy, this.dWxNoTilt, this.dWyNoTilt, ...
                                this.dWxUnwrapped, this.dWyUnwrapped,z1,CCDrot]  =  ...
                                ...
                                lsianalyze.utils.extract2DTDphase(...
                                this.ceInt, this.ceIntMeta, dPhaseSteps, u8FourierTransformType, ...
                                u8UnwrapAlgorithm, this.dAnalysisRegion,this.dAnalysisRegion2,dLowPass,FilterType,....
                                NA,this.uieScanRange.get(),z2_mm,detSize);
                    end
                    
                    
                case 'Fourier Domain'
                    [this.dWx, this.dWy, this.dWxNoTilt, this.dWyNoTilt, this.dWxUnwrapped, this.dWyUnwrapped,z1,CCDrot]    = lsianalyze.utils.extractPhaseFD(...
                        this.ceInt{1}, FilterType, FilterWidth, u8UnwrapAlgorithm, this.dAnalysisRegion,this.dAnalysisRegion2,NA,this.uieScanRange.get(),z2_mm,detSize);
            end
            this.uieFDZ1.set(z1);
            this.uitGramRot.set(strcat('Gram Rot:',num2str(CCDrot),'degree'));
            %             s1=this.dWx;
            %             s2=this.dWy;
            %             s1(this.dAnalysisRegion2==0)=NaN;
            %             s2(this.dAnalysisRegion2==0)=NaN;
            %             figure(10),mesh(s1+s2')
            sig=5;
            this.dAnalysisRegion2(abs(this.dWxNoTilt)>sig*std(this.dWxNoTilt(this.dAnalysisRegion~=0&this.dAnalysisRegion2~=0)))=0;
            this.dAnalysisRegion2(abs(this.dWyNoTilt)>sig*std(this.dWyNoTilt(this.dAnalysisRegion~=0&this.dAnalysisRegion2~=0)))=0;
            this.dWxNoTilt(abs(this.dWxNoTilt)>sig*std(this.dWxNoTilt(this.dAnalysisRegion~=0&this.dAnalysisRegion2~=0)))=0;
            this.dWyNoTilt(abs(this.dWyNoTilt)>sig*std(this.dWyNoTilt(this.dAnalysisRegion~=0&this.dAnalysisRegion2~=0)))=0;
            this.dWxNoTilt=this.dWxNoTilt/2/pi;
            this.dWyNoTilt=this.dWyNoTilt/2/pi;
            % Make phase tab active:
            this.uitgAxesDisplay.selectTabByIndex(this.U8PROBEOBJECT);
            
            % Plot wavefronts on phase tab
            this.replot(this.U8PROBEOBJECT, []);
            
            % Set state:
            this.setState(this.U8STATE_PHASE_PROCESSED);
        end
        
        
        function reconstruct(this, cReconstructionType)
            if isempty(this.dWx)||isempty(this.dWy)
                msgbox('Please process phase first!', 'Error');
                return;
            end
            RimmerType  = this.uipRimmerType.getSelectedIndex();
            FittingType  = this.uipFittingType.getSelectedIndex();
            N           = this.uieRes.get();
            lambda_um   = this.uieLambda.get()/1000;
            T_um        = this.uieScanRange.get();
            %z1_um       = this.uiez1.get();
            z1_um       = this.uieFDZ1.get();
            z2_mm       = this.uiez2.get();
            alpha       = this.uieGratTilt.get() * pi/180;
            gamma       = this.uieDetTilt.get() * pi/180;
            NA          = this.uieNA.get();
            detSize     = this.uieDetSize.get();
            dZernOrder  = this.uieNZernikes.get();
            zernCouples = eval(this.uieZrn.get());
            dZrnRef     = zeros(dZernOrder,1);
            for j=1:size(zernCouples,1)
                if zernCouples(j,1)>0&&zernCouples(j,1)<=dZernOrder
                    dZrnRef(zernCouples(j,1))=zernCouples(j,2);
                end
            end
            
            
            ceAnalysisParas={lambda_um*1000,T_um,NA,z2_mm,this.uieGratTilt.get(),this.uieDetTilt.get(),...
                this.uieGlobalRot.get(),detSize,this.uieCenterObstruction.get(),this.uieBinning.get(),...
                this.uipDataType.getSelectedIndex(),eval(this.uieCCDCenter.get()),eval(this.uieObsOffset.get()), this.uipSelectMask.getSelectedIndex(),...
                this.uitgAnalysisDomain.getSelectedTabName(),this.uipUnwrapEngine.getSelectedIndex(),...
                z1_um,eval(this.uiePhaseStepsTD.get()),this.uipFTType.getSelectedIndex(),...
                this.uieFilterWidth.get(),this.uipWindowFunction.getSelectedIndex(),...
                this.uitgReconstructionType.getSelectedTabName(),this.uieNZernikes.get(),...
                this.uipRimmerType.getSelectedIndex(),this.uicbAutoLoadRim.get(),...
                this.uicbAutoLoadDB.get(),this.uicbAutoLoadFT.get(),this.uieNZernikesBasis.get(),...
                this.uicbRemoveTiltX.get(),this.uicbRemoveTiltY.get(),this.uicbRemoveDefocus.get(),...
                this.uicbScaledNW.get()};
            
            if N~=length(this.dWx)
                N=length(this.dWx);
            end
            z1r=round(z1_um);
            if z1r==0&&z1_um>0
                z1r=0.1;
            elseif z1r==0&&z1_um<0
                z1r=-0.1;
            elseif z1_um==0
                z1r=0;
            end
            
            % Check to see if bases folder exists:
            cBasisFolder = fullfile(this.cAppPath , '..','..' ,'bases');
            if  isempty(dir(cBasisFolder))
                mkdir(cBasisFolder)
            end
            
            switch cReconstructionType
                case 'Rimmer'
                    lAutoLoad = this.uicbAutoLoadRim.get();
                    if lAutoLoad
                        cBasisName = sprintf('RimmerBasis,%d,%0.1f,%0.3f,%0.3f,%0.3f,%0.3f,%0.3f,%0.2f,%0.3f', ...
                            N, lambda_um, T_um, z1r, z2_mm, alpha, gamma, NA, detSize);
                        cBasisName(cBasisName == '.') = '_';
                        cBasisName=strcat(cBasisName,'.mat');
                        % Look for a basis using these params
                        cBasisPath = fullfile(this.cAppPath ,'..', '..','bases', cBasisName);
                        ceFls = dir(cBasisPath);
                        if ~isempty(ceFls)
                            stBasis=load(cBasisPath);
                        else
                            stBasis = [];
                        end
                        [dW,dZrn,ceBasisVectors,RMSFit] = ...
                            lsianalyze.utils.reconstructByRimmer(...
                            this.dWx, this.dWy, RimmerType, this.dAnalysisRegion,this.dAnalysisRegion2, ...
                            lambda_um, T_um, NA, z1_um, z2_mm, N,dZernOrder,alpha,gamma,detSize,stBasis,this.uicbScaledNW.get());
                        save(cBasisPath,'ceBasisVectors');
                    else
                        stBasis = [];
                        [dW,dZrn,~,RMSFit] = ...
                            lsianalyze.utils.reconstructByRimmer(...
                            this.dWx, this.dWy, RimmerType, this.dAnalysisRegion,this.dAnalysisRegion2, ...
                            lambda_um, T_um, NA, z1_um, z2_mm, N,dZernOrder,alpha,gamma,detSize,stBasis,this.uicbScaledNW.get());
                    end
                    
                case 'Derivative Basis'
                    lAutoLoad = this.uicbAutoLoadDB.get();
                    
                    if lAutoLoad
                        % Generate basis name:
                        cBasisName = sprintf('BetaBasis,%d,%d,%0.1f,%0.3f,%0.3f,%0.3f,%0.3f,%0.3f,%0.2f,%0.3f', ...
                            N, dZernOrder, lambda_um, T_um, z1r, z2_mm, alpha, gamma, NA, detSize);
                        
                        cBasisName(cBasisName == '.') = '_';
                        cBasisName=strcat(cBasisName,'.mat');
                        
                        % Look for a basis using these params
                        
                        
                        cBasisPath = fullfile(this.cAppPath , '..','..' ,'bases', cBasisName);
                        ceFls = dir(cBasisPath);
                        if ~isempty(ceFls)
                            stBasis=load(cBasisPath);
                        else
                            stBasis = [];
                        end
                        
                        [dZrn,ceBasisVectors,RMSFit]    = ...
                            lsianalyze.utils.reconstructByDerivativeBasis( ...
                            this.dWx, this.dWy, dZernOrder, this.dAnalysisRegion,this.dAnalysisRegion2, ...
                            lambda_um, NA, T_um, z1_um, z2_mm, N, alpha, gamma, detSize,stBasis,this.uicbScaledNW.get());
                        save(cBasisPath,'ceBasisVectors');
                        
                    else
                        % Call vanilla DBRecon
                        stBasis = [];
                        [dZrn,~,RMSFit]   = ...
                            lsianalyze.utils.reconstructByDerivativeBasis( ...
                            this.dWx, this.dWy, dZernOrder, this.dAnalysisRegion,this.dAnalysisRegion2, ...
                            lambda_um, NA, T_um, z1_um, z2_mm, N, alpha, gamma,detSize,stBasis,this.uicbScaledNW.get() );
                    end
                case 'Fourier'
                    lAutoLoad = this.uicbAutoLoadFT.get();
                    orthogonalization = this.uicbOrthogonalization.get();
                    if lAutoLoad
                        % Generate basis name:
                        cBasisName = sprintf('FourierBasis,%d,%d,%0.1f,%0.3f,%0.3f,%0.3f,%0.3f,%0.3f,%0.2f,%0.3f', ...
                            N, dZernOrder, lambda_um, T_um, z1r, z2_mm, alpha, gamma, NA, detSize);
                        
                        cBasisName(cBasisName == '.') = '_';
                        cBasisName=strcat(cBasisName,'.mat');
                        
                        % Look for a basis using these params
                        cBasisPath = fullfile(this.cAppPath ,'..','..','bases', cBasisName);
                        ceFls = dir(cBasisPath);
                        if ~isempty(ceFls)
                            stBasis=load(cBasisPath);
                            fprintf('Found a Fourier basis for given parameters, using it!\n');
                        else
                            stBasis = [];
                            fprintf('No Fourier basis with these parameters was found, generating one now...\n');
                        end
                        
                        [FTPhase,dZrn,ceBasisVectors,RMSFit]    = ...
                            lsianalyze.utils.reconstructByFourier( ...
                            this.dWx, this.dWy, dZernOrder, this.dAnalysisRegion,this.dAnalysisRegion2, ...
                            lambda_um, NA, T_um, z1_um, z2_mm, N, alpha, gamma, detSize,stBasis,this.uicbScaledNW.get(),...
                            this.uicbRemoveTiltX.get(),this.uicbRemoveTiltY.get(),this.uicbRemoveDefocus.get(),FittingType,orthogonalization); %#ok<ASGLU>
                        
                        save(cBasisPath,'ceBasisVectors');
                        
                    else
                        % Call vanilla DBRecon
                        stBasis = [];
                        [FTPhase,dZrn,~,RMSFit]   = ...
                            lsianalyze.utils.reconstructByFourier( ...
                            this.dWx, this.dWy, dZernOrder, this.dAnalysisRegion,this.dAnalysisRegion2, ...
                            lambda_um, NA, T_um, z1_um, z2_mm, N, alpha, gamma,detSize,stBasis,this.uicbScaledNW.get(),...
                            this.uicbRemoveTiltX.get(),this.uicbRemoveTiltY.get(),this.uicbRemoveDefocus.get(),FittingType,orthogonalization);
                    end
            end
            %             dZrns=zerndecomp((abdom(pinhole(500), [0;dZrn])), dZernOrder );
            %             dZrn=dZrns(2:end);
            dZrn(1:3)  = dZrn(1:3).*[~this.uicbRemoveTiltX.get();~this.uicbRemoveTiltY.get();~this.uicbRemoveDefocus.get()];
            dZrnRef(1:3)  = dZrnRef(1:3).*[~this.uicbRemoveTiltX.get();~this.uicbRemoveTiltY.get();~this.uicbRemoveDefocus.get()];
            % set to waves for now:
            
            this.dZernike = dZrn;
            this.dZernikeResidual = dZrn(:)-dZrnRef;
            %             this.dZernike=  dZrn*lambda_um*1000;
            
            this.dReconstructed=zeros(2*round(this.dBeamWidthEstPx));
            for i=1:dZernOrder
                this.dReconstructed=this.dReconstructed+dZrn(i)*zgen(2*round(this.dBeamWidthEstPx),i);
            end
            if strcmp(cReconstructionType,'Fourier')
                this.dReconstructed=FTPhase;
                this.dResult(1)=std(this.dReconstructed(~isnan(this.dReconstructed)));
            else
                this.dReconstructed=pad2(this.dReconstructed,N,N).*(this.dAnalysisRegion).*(this.dAnalysisRegion2);
                this.dResult(1)=std(this.dReconstructed(rot90(this.dAnalysisRegion)==1));
            end
            
            this.dResult(2)=RMSFit;
            %              this.dResult(2)=RMSFit*lambda_um*1000;
            
            % Make phase tab active:
            this.uitgAxesDisplay.selectTabByIndex(this.U8RECONSTRUCTION);
            
            % Plot wavefronts on phase tab
            this.replot(this.U8RECONSTRUCTION, []);
            %analysislog
            if ~isequal(this.ceAnalysisPara,ceAnalysisParas)&&strcmp(this.uitgSelectDataSource.getSelectedTabName(),'Load P/S Series')
                this.ceAnalysisPara=ceAnalysisParas;
                this.replot(this.U8LOG, []);
                this.SaveAnalysisLog();%save analysis log
            end
            % Set state:
            this.setState(this.U8STATE_RECONSTRUCTED);
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
            z1_um       = this.uiez1.get();
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
            
            
            % Generate Wx and Wy using 2-ray method:
            N           = this.uieRes.get();
            lambda_um   = this.uieLambda.get()/1000;
            T_um        = this.uieScanRange.get();
            z1_um       = this.uiez1.get();
            z2_mm       = this.uiez2.get();
            alpha       = this.uieGratTilt.get() * pi/180;
            gamma       = this.uieDetTilt.get() * pi/180;
            zernCouples = eval(this.uieZrn.get());
            NA          = this.uieNA.get();
            Np          = this.uieNp.get();
            s2s         = this.uieShotToShot.get();
            detSize     = this.uieDetSize.get();
            CenterObstruction     = this.uieCenterObstruction.get();
            detCurve=this.uieDetCurve.get();
            FlareLevel=this.uieFlareLevel.get();
            MSFN=this.uieMSFN.get();
            airflow=this.uieAirflow.get();
            SelectMask  = this.uipSelectMask.getSelectedIndex();
            
            % Higher orders:
            d20Strength = this.uie2ndOrderStrength.get();
            d11Strength = this.uie11orderStrength.get();
            
            zfn         =   lsianalyze.utils.generateZernikeFunction(zernCouples,N,1); % generate zernike function form based on zernike couples
            dMSFN_Phase =   lsianalyze.utils.generateMSFN(MSFN,2000,10,100); % generate MSFN
            dNonlinearitys = eval(this.uieNonlinearity.get());
            dNonlinearity = dNonlinearitys(1)*randn(N);
            dNonlinearityPower = dNonlinearitys(2);
            if lComputePSStack % True if simulating "stack"
                %                 simInts = {};
                u8NumPhaseSteps = this.uieScanSteps.get();
                simInts = cell(u8NumPhaseSteps, 2);
                this.uipbExposureProgress.set(0);
                
                
                dPhaseShifts = eval(this.uiePhaseStepsSim.get());
                % X steps
                %                  phasedrift=linspace(0,pi/3,u8NumPhaseSteps^2);
                %                  phasedrift=reshape(phasedrift,u8NumPhaseSteps,u8NumPhaseSteps);
                %                  phasedrift(:,2:2:end)=flipud(phasedrift(:,2:2:end));
                
                % Drift effects:
                dNTotalSteps = u8NumPhaseSteps^2;
                
                
                dZDriftTotal = this.uieZLinearDrift.get() / 1000; % (um)
                dXDriftTotal = this.uieXLinearDrift.get() / this.uieScanRange.get() / 1000 * 2 * pi; % (rad)
                dYDriftTotal = this.uieYLinearDrift.get() / this.uieScanRange.get() / 1000 * 2 * pi; % (rad)
                
                dStageUncertainty = this.uiePhaseShiftingError.get() / this.uieScanRange.get() / 1000  * 2 * pi; %rad
                
                
                if dZDriftTotal > 0 % Then compute a drift map (in um):
                    dZDriftStep = dZDriftTotal/dNTotalSteps;
                else
                    dZDriftStep = 0;
                end
                if dXDriftTotal > 0 % Then compute a drift map (in rad):
                    dXDriftStep = dXDriftTotal/dNTotalSteps;
                else
                    dXDriftStep = 0;
                end
                if dYDriftTotal > 0 % Then compute a drift map (in rad):
                    dYDriftStep = dYDriftTotal/dNTotalSteps;
                else
                    dYDriftStep = 0;
                end
                
                % Increment cumulative drifts
                dXCumDrift = (0:dNTotalSteps-1) * dXDriftStep;
                dYCumDrift = (0:dNTotalSteps-1) * dYDriftStep;
                dZCumDrift = (0:dNTotalSteps-1) * dZDriftStep;
                dXCumDrift = reshape(dXCumDrift,[u8NumPhaseSteps,u8NumPhaseSteps]);
                dYCumDrift = reshape(dYCumDrift,[u8NumPhaseSteps,u8NumPhaseSteps]);
                dZCumDrift = reshape(dZCumDrift,[u8NumPhaseSteps,u8NumPhaseSteps]);
                
                dXPhaseShifts = dPhaseShifts(:,1);
                dYPhaseShifts = dPhaseShifts(:,2);
                
                % Aligment:
                dGratTiltEr = this.uieGratingTiltErr.get() * pi/180;
                dDetTiltEr = this.uieDetectorTiltErr.get() * pi/180;
                
                
                try
                    this.uipbExposureProgress.set(0);
                    drawnow;
                    parfor mk=1:length(dXPhaseShifts)*length(dYPhaseShifts)
                        k=mod(mk,length(dYPhaseShifts));
                        if k==0
                            k=length(dYPhaseShifts);
                        end
                        m=ceil(mk/length(dYPhaseShifts));
                        Wint = lsianalyze.utils.simulateTwoRayMethod_deferredZ(...
                            N, lambda_um, NA, T_um, z1_um + dZCumDrift(k,m), z2_mm, ...
                            alpha + dGratTiltEr, gamma + dDetTiltEr, zernCouples,zfn, Np,...
                            s2s, d11Strength, d20Strength, detCurve, FlareLevel,dMSFN_Phase,...
                            airflow,detSize,CenterObstruction,SelectMask,...
                            [dXPhaseShifts(k) + dXCumDrift(k,m) + dStageUncertainty*randn(1), ...
                            dYPhaseShifts(m) + dYCumDrift(k,m) + dStageUncertainty*randn(1)]);
                        fileName=['Interferogram',num2str(m),'_',num2str(k),'.int'];
                        fid=fopen(fileName,'w');
                        fwrite(fid,Wint,'double');
                        fclose(fid);
                    end
                    for m=1:length(dXPhaseShifts)
                        for k=1:length(dYPhaseShifts)
                            fileName=['Interferogram',num2str(m),'_',num2str(k),'.int'];
                            fid=fopen(fileName,'r');
                            Wint=reshape(fread(fid,N^2,'double'),[N,N]);
                            fclose(fid);
                            simInts{k, m} = Wint;
                        end
                    end
                    this.uipbExposureProgress.set(1);
                    delete('*.int');
                catch
                    for m = 1:length(dXPhaseShifts)
                        for k = 1:length(dYPhaseShifts)
                            %
                            Wint = lsianalyze.utils.simulateTwoRayMethod_deferredZ(...
                                N, lambda_um, NA, T_um, z1_um + dZCumDrift(k,m), z2_mm, ...
                                alpha + dGratTiltEr, gamma + dDetTiltEr, zernCouples,zfn, Np,...
                                s2s, d11Strength, d20Strength, detCurve, FlareLevel,dMSFN_Phase,...
                                airflow,detSize,CenterObstruction,SelectMask,...
                                [dXPhaseShifts(k) + dXCumDrift(k,m) + dStageUncertainty*randn(1), ...
                                dYPhaseShifts(m) + dYCumDrift(k,m) + dStageUncertainty*randn(1)]);
                            
                            simInts{k, m} = Wint;
                            this.uipbExposureProgress.set(k/length(dXPhaseShifts)/length(dYPhaseShifts)+(m-1)/length(dXPhaseShifts));
                            drawnow
                        end
                    end
                end
                % Y steps
                simInts=simInts';
                % make sure there is no negetive value in interferogram
                % plus CCD nonlinearity
                minValue=0;
                maxValue=0;
                for m = 1:length(dXPhaseShifts)
                    for k = 1:length(dYPhaseShifts)
                        if -min(simInts{k, m}(:))>minValue
                            minValue=min(simInts{k, m}(:));
                        end
                        if max(simInts{k, m}(:))>maxValue
                            maxValue=max(simInts{k, m}(:));
                        end
                    end
                end
                if minValue<0
                    for m = 1:length(dXPhaseShifts)
                        for k = 1:length(dYPhaseShifts)
                            %                             simInts{k, m}=(1-dNonlinearity).*((simInts{k, m}-minValue)/(maxValue-minValue))+dNonlinearity.*((simInts{k, m}-minValue)/(maxValue-minValue)).^dNonlinearityPower;
                            %                               simInts{k, m}=(simInts{k, m}-minValue)/(maxValue-minValue) + dNonlinearity;
                            %                               simInts{k, m}=(simInts{k, m}-minValue)/(maxValue-minValue) +dNonlinearity.*((simInts{k, m}-minValue)/(maxValue-minValue)).^dNonlinearityPower;
                            simInts{k, m}=((simInts{k, m}-minValue)/(maxValue-minValue)).^(1+dNonlinearityPower);
                        end
                    end
                else
                    for m = 1:length(dXPhaseShifts)
                        for k = 1:length(dYPhaseShifts)
                            %                             simInts{k, m}=(1-dNonlinearity).*(simInts{k, m}/maxValue)+dNonlinearity.*(simInts{k, m}/maxValue).^dNonlinearityPower;
                            %                             simInts{k, m}=simInts{k, m}/maxValue + dNonlinearity;
                            %                             simInts{k, m}=simInts{k, m}/maxValue +dNonlinearity.*(simInts{k, m}/maxValue).^dNonlinearityPower;
                            simInts{k, m}=(simInts{k, m}/maxValue).^(1+dNonlinearityPower);
                        end
                    end
                end
                
                this.handleLoadData(simInts, {'sim'});
                this.dAnalysisRegion(simInts{1,1}==0)=0;
                this.uipbExposureProgress(1);
                drawnow
            else % Simulate a singe image:
                tic
                Wint = lsianalyze.utils.simulateTwoRayMethod_deferredZ(...
                    N, lambda_um, NA, T_um, z1_um, z2_mm, ...
                    alpha, gamma, zernCouples,zfn, Np, s2s, d11Strength, d20Strength, ...
                    detCurve, FlareLevel,dMSFN_Phase,airflow,detSize,CenterObstruction, SelectMask,[0,0]);
                % make sure there is no negetive value in interferogram
                % plus CCD nonlinearity
                if min(Wint(:))<0
                    %                     Wint=(1-dNonlinearity).*((Wint-min(Wint(:)))/(max(Wint(:))-min(Wint(:))))+dNonlinearity.*((Wint-min(Wint(:)))/(max(Wint(:))-min(Wint(:)))).^dNonlinearityPower;
                    %                     Wint=(Wint-min(Wint(:)))/(max(Wint(:))-min(Wint(:))) + dNonlinearity;
                    Wint=((Wint-min(Wint(:)))/(max(Wint(:))-min(Wint(:)))).^(1+dNonlinearityPower);
                else
                    %                     Wint=(1-dNonlinearity).*(Wint/max(Wint(:)))+dNonlinearity.*(Wint/max(Wint(:))).^dNonlinearityPower;
                    %                     Wint=Wint/max(Wint(:)) + dNonlinearity;
                    Wint=(Wint/max(Wint(:))).^(1+dNonlinearityPower);
                end
                
                this.handleLoadData({Wint}, {'sim'});
                this.dAnalysisRegion(Wint==0)=0;
                fprintf('Single image took %s\n', s2f(toc));
            end
            %             dmask=zeros(N);
            %             z2_um       = z2_mm * 1000;
            %             halfD=detSize/2*1e3;
            %             NAs=NA*halfD/(tan(asin(NA))*z2_um);
            %             xIdxs        = z2_um*tan(asin(NAs*linspace(-1, 1, N)));
            %             yIdxs        = z2_um*tan(asin(NAs*linspace(-1, 1, N)));
            %             [X, Y]      = meshgrid(xIdxs,yIdxs);
            %             dmask(X.^2+Y.^2<=(z2_um*tan(asin(NA)))^2)=1;
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
            z1_mm       = this.uiez1.get()/1000;
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
                        this.dAnalysisRegion = pinhole(2*round(this.dBeamWidthEstPx - this.dShearPix), sr, sc)-...
                            circshift(pinhole(2*round(this.dObstructionWidthEstPx + this.dObstructionShearPix), sr, sc), ...
                            [ dObsOffset(2),  dObsOffset(1)]);
                    else
                        this.dAnalysisRegion = pinhole(2*round(this.dBeamWidthEstPx - this.dShearPix), sr, sc)-...
                            pinhole(2*round(this.dObstructionWidthEstPx + this.dObstructionShearPix), sr, sc);
                    end
                else
                    this.dAnalysisRegion = pinhole(2*round(this.dBeamWidthEstPx - this.dShearPix), sr, sc);
                end
            else
                if this.dObstructionWidthEstPx~=0
                    if ~isempty(dObsOffset)
                        this.dAnalysisRegion = pinhole(2*ceil(this.dBeamWidthEstPx - this.dShearPix)-1, sr, sc)-...
                            circshift(pinhole(2*floor(this.dObstructionWidthEstPx + this.dObstructionShearPix)+1, sr, sc), ...
                            [ dObsOffset(2),  dObsOffset(1)]);
                    else
                        this.dAnalysisRegion = pinhole(2*ceil(this.dBeamWidthEstPx - this.dShearPix)-1, sr, sc)-...
                            pinhole(2*floor(this.dObstructionWidthEstPx + this.dObstructionShearPix)+1, sr, sc);
                    end
                else
                    this.dAnalysisRegion = pinhole(2*ceil(this.dBeamWidthEstPx - this.dShearPix)-1, sr, sc);
                end
            end
            %remove lines
            this.dAnalysisRegion2=ones(sr,sc);
            %             LineWidth=2;
            %             HalfAngle=27;
            %             [x,y]=meshgrid(linspace(-1,1,sr));
            %             this.dAnalysisRegion2((abs(x-tan(HalfAngle/180*pi)*y)/sqrt(1+tan(HalfAngle/180*pi)^2)<LineWidth*2/sr|...
            %                 abs(x+tan(HalfAngle/180*pi)*y)/sqrt(1+tan(HalfAngle/180*pi)^2)<LineWidth*2/sr)&y<=0)=0;
        end
        
        function computeMET5Mask(this)
            lambda_um   = this.uieLambda.get()/1000;
            T_um        = this.uieScanRange.get();
            z2_mm       = this.uiez2.get();
            z1_mm       = this.uiez1.get()/1000;
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
            
            
            %
            %             this.dBeamWidthEstPx    = round(tan(asin(NA))*(z2_mm) * sr / (detSize));
            %             this.dShearPix          = ceil(this.dBeamWidthEstPx - ...
            %                                         tan(asin(NA) - lambda_um/T_um)*(z2_mm) * sr / (detSize));
            %             this.dObstructionWidthEstPx = round(CenterObstruction*tan(asin(NA))*(z2_mm) * sr / (detSize));
            %             this.dObstructionShearPix =round(tan(asin(lambda_um/T_um))*z2_mm* sr / detSize);
            %             % Set default mask:
            %             if this.dObstructionWidthEstPx~=0
            %                 this.dAnalysisRegion = pinhole(2*(this.dBeamWidthEstPx - this.dShearPix), sr, sc)-...
            %                     pinhole(2*(this.dObstructionWidthEstPx + this.dObstructionShearPix), sr, sc);
            %             else
            %                 this.dAnalysisRegion = pinhole(2*(this.dBeamWidthEstPx - this.dShearPix), sr, sc);
            %             end
            %             %remove lines
            %             this.dAnalysisRegion2=ones(sr,sc);
            %             LineWidth=4;
            %             HalfAngle=45;
            %             [x,y]=meshgrid(linspace(-1,1,sr));
            %             this.dAnalysisRegion2((abs(x-tan(HalfAngle/180*pi)*y)/sqrt(1+tan(HalfAngle/180*pi)^2)<LineWidth*2/sr|...
            %                 abs(x+tan(HalfAngle/180*pi)*y)/sqrt(1+tan(HalfAngle/180*pi)^2)<LineWidth*2/sr)&y>=0)=0;
        end
        function computeElipticalMask(this)
            lambda_um   = this.uieLambda.get()/1000;
            T_um        = this.uieScanRange.get();
            z2_mm       = this.uiez2.get();
            z1_mm       = this.uiez1.get()/1000;
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
            %             LineWidth=2;
            %             HalfAngle=27;
            %             [x,y]=meshgrid(linspace(-1,1,sr));
            %             this.dAnalysisRegion2((abs(x-tan(HalfAngle/180*pi)*y)/sqrt(1+tan(HalfAngle/180*pi)^2)<LineWidth*2/sr|...
            %                 abs(x+tan(HalfAngle/180*pi)*y)/sqrt(1+tan(HalfAngle/180*pi)^2)<LineWidth*2/sr)&y<=0)=0;
        end
        
        % Position recall getters
        function store = prGetters(this, src)
            store = [];
            switch src
                case this.prControlsSetup
                    ceFields = {'uieLambda', 'uieScanRange', 'uieNA', 'uiez2', ...
                        'uieGratTilt', 'uieDetTilt', 'uieGlobalRot', 'uieDetSize',...
                        'uieCenterObstruction','uieBinning'};
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
                        'uieCenterObstruction','uieBinning'};
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
                        % Rotate interferograms here
                        this.ceIntR{k,t} = imrotate(ceData{k,t},dGlobalRot,'crop');
                    else
                        this.ceIntR{k,t} = ceData{k,t};
                    end
                    if ~isempty(ceMeta) && strcmp(ceMeta{1},'shifted')==1
                        this.ceIntR{k,t} = this.ceInt{k,t};
                        this.ceIntMeta={};
                    end
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
            try
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
            catch
            end
            
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
                        
                        % Plot main interferogram
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
                                
                                dXCirc2 = (this.dBeamWidthEstPx - this.dShearPix) * cos(dTh) + dMidX;
                                dYCirc2 = (this.dBeamWidthEstPx - this.dShearPix) * sin(dTh)*1 + dMidY;
                                
                                this.hsaInterferogram.plot(dXCirc2, dYCirc2, 'y', 'lineWidth', 2);
                                if isempty(dObsOffset)
                                    dObsOffset=[0,0];
                                end
                                if this.uieCenterObstruction.get() > 0
                                    dXCirc3 = this.dObstructionWidthEstPx * cos(dTh) + dMidX+dObsOffset(1);
                                    dYCirc3 = this.dObstructionWidthEstPx * sin(dTh)*1 + dMidY+dObsOffset(2);
                                    
                                    this.hsaInterferogram.plot(dXCirc3, dYCirc3, 'r', 'lineWidth', 2);
                                    
                                    dXCirc4 = (this.dObstructionWidthEstPx + this.dObstructionShearPix) * cos(dTh) + dMidX+dObsOffset(1);
                                    dYCirc4 = (this.dObstructionWidthEstPx + this.dObstructionShearPix) * sin(dTh)*1 + dMidY+dObsOffset(2);
                                    
                                    this.hsaInterferogram.plot(dXCirc4, dYCirc4, 'b', 'lineWidth', 2);
                                end
                                dRedPx = (this.dBeamWidthEstPx - this.dShearPix);
                                this.hsaInterferogram.plot((dMidX + dRedPx)*[1, 1], [1, sr], 'c', 'lineWidth', 0.5);
                                this.hsaInterferogram.plot((dMidX - dRedPx)*[1, 1], [1, sr], 'c', 'lineWidth', 0.5);
                                this.hsaInterferogram.plot([1, sc], (dMidY + dRedPx)*[1, 1], 'c', 'lineWidth', 0.5);
                                this.hsaInterferogram.plot([1, sc], (dMidY - dRedPx)*[1, 1], 'c', 'lineWidth', 0.5);
                            else
                                dTh = linspace(0, 2*pi, 41);
                                dXCirc1 = this.dBeamWidthEstPx * cos(dTh) + dMidX;
                                dYCirc1 = this.dBeamWidthEstPx * sin(dTh)*0.5714 + dMidY;
                                
                                this.hsaInterferogram.plot(dXCirc1, dYCirc1, 'g', 'lineWidth', 2);
                                
                                dXCirc2 = (this.dBeamWidthEstPx - this.dShearPix) * cos(dTh) + dMidX;
                                dYCirc2 = (this.dBeamWidthEstPx - this.dShearPix) * sin(dTh)*0.5714 + dMidY;
                                
                                this.hsaInterferogram.plot(dXCirc2, dYCirc2, 'y', 'lineWidth', 2);
                                if isempty(dObsOffset)
                                    dObsOffset=[0,0];
                                end
                                if this.uieCenterObstruction.get() > 0
                                    dXCirc3 = this.dObstructionWidthEstPx * cos(dTh) + dMidX+dObsOffset(1);
                                    dYCirc3 = this.dObstructionWidthEstPx * sin(dTh)*1 + dMidY+dObsOffset(2);
                                    
                                    this.hsaInterferogram.plot(dXCirc3, dYCirc3, 'r', 'lineWidth', 2);
                                    
                                    dXCirc4 = (this.dObstructionWidthEstPx + this.dObstructionShearPix) * cos(dTh) + dMidX+dObsOffset(1);
                                    dYCirc4 = (this.dObstructionWidthEstPx + this.dObstructionShearPix) * sin(dTh)*1 + dMidY+dObsOffset(2);
                                    
                                    this.hsaInterferogram.plot(dXCirc4, dYCirc4, 'b', 'lineWidth', 2);
                                end
                                
                                
                                dRedPx = (this.dBeamWidthEstPx - this.dShearPix);
                                this.hsaInterferogram.plot((dMidX + dRedPx)*[1, 1], [0, sr], 'c', 'lineWidth', 0.5);
                                this.hsaInterferogram.plot((dMidX - dRedPx)*[1, 1], [0, sr], 'c', 'lineWidth', 0.5);
                                this.hsaInterferogram.plot([0, sc], (dMidY + dRedPx*0.5714)*[1, 1], 'c', 'lineWidth', 0.5);
                                this.hsaInterferogram.plot([0, sc], (dMidY - dRedPx*0.5714)*[1, 1], 'c', 'lineWidth', 0.5);
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
                        detSize_um  = this.uieDetSize.get()*1000; 
                        N           = this.uieRes.get();
                        dc_um       = detSize_um/N; % detector pixel pitch
                        if strcmp(propagator,'fourier')
                            do_um = lambda_um*z_um/N/dc_um;
                        else
                            do_um = dc_um; % object pixel pitch
                        end
                        [K,L] = size(this.dObject);
                        xp_mm = [1:N]*do_um/1000; % object coordinates
                        xo_mm = [1:L]*do_um/1000; % object coordinates
                        yo_mm = [1:K]*do_um/1000; % object coordinates
                        
                        imagesc(this.haProbeAmp, xp_mm,xp_mm,abs(this.dProbe));colorbar(this.haProbeAmp);axis(this.haProbeAmp,'xy');
                        this.haProbeAmp.Title.String = 'Probe amplitude';this.haProbeAmp.XLabel.String = 'mm';this.haProbeAmp.YLabel.String = 'mm';
                        imagesc(this.haProbePha, xp_mm,xp_mm,atan2(imag(this.dProbe),real(this.dProbe)));colorbar(this.haProbePha);axis(this.haProbePha,'xy');
                        this.haProbePha.Title.String = 'Probe phase';this.haProbePha.XLabel.String = 'mm';this.haProbePha.YLabel.String = 'mm';
                        imagesc(this.haObjectAmp, xo_mm,yo_mm,abs(this.dObject));colorbar(this.haObjectAmp);axis(this.haObjectAmp,'xy');
                        this.haObjectAmp.Title.String = 'Object amplitude';this.haObjectAmp.XLabel.String = 'mm';this.haObjectAmp.YLabel.String = 'mm';
                        imagesc(this.haObjectPha, xo_mm,yo_mm,atan2(imag(this.dObject),real(this.dObject)));colorbar(this.haObjectPha);axis(this.haObjectPha,'xy');
                        this.haObjectPha.Title.String = 'Object phase';this.haObjectPha.XLabel.String = 'mm';this.haObjectPha.YLabel.String = 'mm';
                    case  this.U8GUESS
                        propagator = this.uipPropagator.getOptions{this.uipPropagator.getSelectedIndex()};
                        detSize_um  = this.uieDetSize.get()*1000; 
                        N           = this.uieRes.get();
                        dc_um       = detSize_um/N; % detector pixel pitch
                        if strcmp(propagator,'fourier')
                            do_um = lambda_um*z_um/N/dc_um;
                        else
                            do_um = dc_um; % object pixel pitch
                        end
                        [K,L] = size(this.dObject);
                        xp_mm = [1:N]*do_um/1000; % object coordinates
                        xo_mm = [1:L]*do_um/1000; % object coordinates
                        yo_mm = [1:K]*do_um/1000; % object coordinates
                        
                        imagesc(this.haGuessProbeAmp, xp_mm,xp_mm,abs(this.dProbeGuess));colorbar(this.haGuessProbeAmp);axis(this.haGuessProbeAmp,'xy');
                        this.haGuessProbeAmp.Title.String = 'Probe amplitude';this.haGuessProbeAmp.XLabel.String = 'mm';this.haGuessProbeAmp.YLabel.String = 'mm';
                        imagesc(this.haGuessProbePha, xp_mm,xp_mm,atan2(imag(this.dProbeGuess),real(this.dProbeGuess)));colorbar(this.haGuessProbePha);axis(this.haGuessProbePha,'xy');
                        this.haGuessProbePha.Title.String = 'Probe phase';this.haGuessProbePha.XLabel.String = 'mm';this.haGuessProbePha.YLabel.String = 'mm';
                        imagesc(this.haGuessObjectAmp, xo_mm,yo_mm,abs(this.dObjectGuess));colorbar(this.haGuessObjectAmp);axis(this.haGuessObjectAmp,'xy');
                        this.haGuessObjectAmp.Title.String = 'Object amplitude';this.haGuessObjectAmp.XLabel.String = 'mm';this.haGuessObjectAmp.YLabel.String = 'mm';
                        imagesc(this.haGuessObjectPha, xo_mm,yo_mm,atan2(imag(this.dObjectGuess),real(this.dObjectGuess)));colorbar(this.haGuessObjectPha);axis(this.haGuessObjectPha,'xy');
                        this.haGuessObjectPha.Title.String = 'Object phase';this.haGuessObjectPha.XLabel.String = 'mm';this.haGuessObjectPha.YLabel.String = 'mm';
                    case this.U8RECONSTRUCTION
                        u8ReconUnit               = this.uipReconUnit.getSelectedIndex();
                        switch u8ReconUnit
                            case 1  %unit: Wave
                                this.uitRMS.set(strcat('RMS:',num2str(this.dResult(1)),'Waves'));
                                this.uitRMSFit.set(strcat('RMSFit:',num2str(this.dResult(2)),'Waves'));
                                if ~this.uicbResidualError.get()
                                    bar(this.haZernikeDecomp,[1:(length(this.dZernike))], this.dZernike);
                                else
                                    residaulPhase= abdom(pinhole(100), [0;this.dZernikeResidual]);
                                    RMSE = std(residaulPhase(logical(pinhole(100))));
                                    this.uitRMS.set(strcat('RMSE:',num2str(RMSE),'Waves'));
                                    bar(this.haZernikeDecomp,[1:(length(this.dZernike))], this.dZernikeResidual);
                                end
                                ylabel(this.haZernikeDecomp,'Waves');
                                xlabel(this.haZernikeDecomp,'Zernike terms');
                            case 2  %unit: mWave
                                this.dReconstructed=this.dReconstructed*1000;
                                this.dZernike=this.dZernike*1000;
                                this.dZernikeResidual=this.dZernikeResidual*1000;
                                this.dResult=this.dResult*1000;
                                this.uitRMS.set(strcat('RMS:',num2str(this.dResult(1)),'mWaves'));
                                this.uitRMSFit.set(strcat('RMSFit:',num2str(this.dResult(2)),'mWaves'));
                                if ~this.uicbResidualError.get()
                                    bar(this.haZernikeDecomp,[1:(length(this.dZernike))], this.dZernike);
                                else
                                    residaulPhase= abdom(pinhole(100), [0;this.dZernikeResidual]);
                                    RMSE = std(residaulPhase(logical(pinhole(100))));
                                    this.uitRMS.set(strcat('RMSE:',num2str(RMSE),'mWaves'));
                                    bar(this.haZernikeDecomp,[1:(length(this.dZernike))], this.dZernikeResidual);
                                end
                                ylabel(this.haZernikeDecomp,'mWaves');
                                xlabel(this.haZernikeDecomp,'Zernike terms');
                            case 3  %unit: nm
                                this.dReconstructed=this.dReconstructed*this.uieLambda.get();
                                this.dZernike=this.dZernike*this.uieLambda.get();
                                this.dZernikeResidual=this.dZernikeResidual*this.uieLambda.get();
                                this.dResult=this.dResult*this.uieLambda.get();
                                this.uitRMS.set(strcat('RMS:',num2str(this.dResult(1)),'nm'));
                                this.uitRMSFit.set(strcat('RMSFit:',num2str(this.dResult(2)),'nm'));
                                if ~this.uicbResidualError.get()
                                    bar(this.haZernikeDecomp,[1:(length(this.dZernike))], this.dZernike);
                                else
                                    residaulPhase= abdom(pinhole(100), [0;this.dZernikeResidual]);
                                    RMSE = std(residaulPhase(logical(pinhole(100))));
                                    this.uitRMS.set(strcat('RMSE:',num2str(RMSE),'nm'));
                                    bar(this.haZernikeDecomp,[1:(length(this.dZernike))], this.dZernikeResidual);
                                end
                                ylabel(this.haZernikeDecomp,'nm');
                                xlabel(this.haZernikeDecomp,'Zernike terms');
                        end
                        axis(this.haZernikeDecomp,[1 length(this.dZernike) -inf inf])
                        if ~this.uicbResidualError.get()
                            imagesc(this.haReconstructed, this.dReconstructed);colorbar(this.haReconstructed);axis(this.haReconstructed,'xy');
                        else
                            imagesc(this.haReconstructed, residaulPhase);colorbar(this.haReconstructed);axis(this.haReconstructed,'xy');
                        end
                        this.haReconstructed.Title.String = 'Reconstructed wavefront';
                        this.haZernikeDecomp.Title.String = 'Zernike decomposition';
                        this.uitRMS.setFontSize(12);
                        this.uitRMSFit.setFontSize(12);
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
            this.processPhase(this.uitgAnalysisDomain.getSelectedTabName());
            
            % Reconstruct:
            this.reconstruct(this.uitgReconstructionType.getSelectedTabName());
            
        end
        
        function setSimParams(this, varargin)
            for k = 1:2:length(varargin)
                paramName = varargin{k};
                paramValue = varargin{k+1};
                switch paramName
                    case 'nZern' % sets all zern values
                        this.uieNZernikes.set(paramValue);
                        this.uieNZernikesBasis.set(paramValue);
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
                    case 'MSFN'
                        this.uieMSFN.set(paramValue);
                    case 'airflow'
                        this.uieAirflow.set(paramValue);
                    case 'lowPass'
                        this.uieLowPass.set(paramValue);
                    case 'filterType'
                        this.uipWindowFunction.setSelectedIndex(uint8(paramValue));
                    case 'fittingType'
                        this.uipFittingType.setSelectedIndex(uint8(paramValue));
                    case 'gratingPitch'
                        this.uieScanRange.set(paramValue);
                    case 'Z_1'
                        this.uiez1.set(paramValue);
                    case 'recType'
                        switch lower(paramValue)
                            case 'fourier'
                                this.uitgReconstructionType.selectTabByName('Fourier');
                            case {'derivative basis', 'derivative-basis'}
                                this.uitgReconstructionType.selectTabByName('Derivative Basis');
                        end
                end
            end
        end
        
        function build(this, dOffsetX, dOffsetY)
            if nargin == 1
                dOffsetX = 0;
                dOffsetY = 0;
            end
            
            % build the main window
            this.hFigure = figure(...
                'name', 'PIE analysis GUI v1.20190129',...
                'Units', 'pixels',...
                'Position', [5 - dOffsetX, 5 - dOffsetY,  this.dWidth, this.dHeight],...
                'handlevisibility','off',... %out of reach gcf
                'numberTitle','off',...
                'Toolbar','none',...
                'Menubar','none');
            
            
            % Build all containers first:
            
            
            drawnow
            
            % Axes
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
            uitReconstruction = this.uitgAxesDisplay.getTabByName('Reconstruction');
            this.haReconstructed = axes('Parent', uitReconstruction, ...
                'Units', 'pixels', ...
                'Position', [200, 360, 560, 480], ...
                'XTick', [], 'YTick', []);
            this.haZernikeDecomp = axes('Parent', uitReconstruction, ...
                'Units', 'pixels', ...
                'Position', [80, 100, 780, 180], ...
                'XTick', [], 'YTick', []);
            
            this.uitRMS.build           (uitReconstruction, 80, 880, 200, 20);
            this.uitRMSFit.build           (uitReconstruction, 300, 880, 200, 20);
            % Axes:Log
            uitLog = this.uitgAxesDisplay.getTabByName('Log');
            this.htLog = uitable(uitLog,'Position',[10 10 930 900],'CellSelectionCallback',@(src, evt)this.cb(src,evt));
            % Controls main:
            this.hpControls = uipanel(...
                'Parent', this.hFigure,...
                'Units', 'pixels',...
                'Title', 'Controls',...
                'FontWeight', 'Bold',...
                'Clipping', 'on',...
                'BorderWidth',1, ...
                'Position', [970 20 620 950] ...
                );
            
            
            % Controls: Experiment setup:
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
            
            % Controls: Data panel:
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
                (this.hpData, 10, 10 + Offset1, 575, 200);drawnow
            
            this.uieCCDCenter.build     (this.hpData, 10, 215+Offset1, 120, 20);
            this.uieObsOffset.build     (this.hpData, 140, 215+Offset1, 120, 20);
            this.uibSelectCenter.build  (this.hpData, 10, 260+Offset1, 80, 20);
            this.uibResetCenter.build (this.hpData, 100, 260+Offset1, 110, 20);
            this.uicbAutoCenter.build  (this.hpData, 220, 260+Offset1, 100, 20);
            
            this.uipSelectMask.build    (this.hpData, 340, 215+Offset1, 240, 20)
            this.uibLoadMask.build      (this.hpData, 470, 260+Offset1, 100, 20);
            
            % Controls:Data:LI
            Offset2=-30;
            uitLISingle = this.uitgSelectDataSource.getTabByName('Load Interferogram');
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
            this.uibLoadPhaseStepsSim.build ...
                (uitSim, 435, 50 + 13+Offset4, 100, 20);
            this.uieZrn.build       (uitSim, 10, 90+Offset4, 410, 20);
            
            this.uibLoadZrn.build   (uitSim, 435, 103+Offset4, 100, 20);
            this.uibSimulate.build  (uitSim, 150, 140+Offset4, 100, 20);
            this.uibSimulateS.build (uitSim, 275, 140+Offset4, 100, 20);
            this.uicbMultiPhaseShifting.build (uitSim, 10, 140+Offset4, 100, 20)
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
            this.uibGenProbeObject.build    (uitProbe, 430, 130+Offsetp, 100, 20);
            this.uicbGuess.build (uitProbe, 430, 10+Offsetp, 100, 20);
            
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
            this.uieMSFN.build                  (uitCSim, 410+Offset8, 90, 90, 25);
            this.uieNonlinearity.build          (uitCSim, 20+Offset8, 130, 90, 25);
            this.uieAirflow.build               (uitCSim, 150+Offset8, 130, 90, 25);
            
            % Controls: Phase
            this.hpPhase = uipanel(...
                'Parent', this.hpControls,...
                'Units', 'pixels',...
                'Title', 'Phase',...
                'FontWeight', 'Bold',...
                'Clipping', 'on',...
                'BorderWidth',1, ...
                'Position', [10 190 600 190] ...
                );
            drawnow
            Offset5=20;
            this.uitgAnalysisDomain.build(this.hpPhase, 10, 20, 400, 160);drawnow
            
            uitFD = this.uitgAnalysisDomain.getTabByName('Fourier Domain');
            uitTD = this.uitgAnalysisDomain.getTabByName('Time Domain');
            
            this.uipUnwrapEngine.build  (this.hpPhase, 415, 90+Offset5, 170, 20);
            this.uibComputePhase.build  (this.hpPhase, 415, 140+Offset5, 160, 20);
            
            this.uieFDZ1.build          (this.hpPhase, 415, 20+Offset5, 80, 20)
            this.uitGramRot.build    (this.hpPhase, 415, 65+Offset5, 200, 20)
            
            % Controls:Phase:FD
            this.uipWindowFunction.build    (uitFD, 10, 20, 200, 20);
            this.uieFilterWidth.build    (uitFD, 10, 70, 80, 20);
            % Controls:Phase:TD
            Offset6=-40;
            this.uiePhaseStepsTD.build  (uitTD, 10, 60+Offset6, 250, 20);
            this.uieLowPass.build  (uitTD, 200, 105+Offset6, 130, 20);
            this.uipFTType.build        (uitTD, 10, 100+Offset6, 150, 20);
            this.uibTDLoadPhaseStepsTD.build ...
                (uitTD, 270, 60 + 13+Offset6, 100, 20);
            
            
            % Controls: Reconstruction
            this.hpReconstruction = uipanel(...
                'Parent', this.hpControls,...
                'Units', 'pixels',...
                'Title', 'Reconstruction',...
                'FontWeight', 'Bold',...
                'Clipping', 'on',...
                'BorderWidth', 1, ...
                'Position', [10 10 600 170] ...
                );
            drawnow
            Offset7=-30;
            this.uitgReconstructionType.build(this.hpReconstruction, 10, 25, 250, 135);
            drawnow
            
            uitRim  = this.uitgReconstructionType.getTabByName('Rimmer');
            uitDB   = this.uitgReconstructionType.getTabByName('Derivative Basis');
            uitFT   = this.uitgReconstructionType.getTabByName('Fourier');
            this.uitgReconstructionType.selectTabByName('Derivative Basis');
            
            
            
            
            this.uipRimmerType.build    (uitRim, 10, 50+Offset7, 200, 20);
            
            this.uicbAutoLoadRim.build  (uitRim, 10, 100+Offset7, 200, 20);
            
            this.uieNZernikesBasis.build(uitDB, 10, 50+Offset7, 100, 20);
            this.uicbAutoLoadDB.build   (uitDB, 10, 100+Offset7, 200, 20);
            
            this.uicbAutoLoadFT.build   (uitFT, 10, 100+Offset7, 200, 20);
            this.uicbOrthogonalization.build   (uitFT, 140, 100+Offset7, 200, 20);
            this.uipFittingType.build    (uitFT, 10, 50+Offset7, 200, 20);
            
            Offset8=20;
            this.uieNZernikes.build     (this.hpReconstruction, 415, 70+Offset8, 80, 20);
            this.uibReconstruct.build   (this.hpReconstruction, 415, 120+Offset8, 160, 20);
            this.uicbRemoveTiltX.build  (this.hpReconstruction, 280, 40+Offset8, 100, 20);
            this.uicbRemoveTiltY.build  (this.hpReconstruction, 280, 80+Offset8, 100, 20);
            this.uicbRemoveDefocus.build(this.hpReconstruction, 280, 120+Offset8, 120, 20);
            
            this.uicbScaledNW.build  (this.hpReconstruction, 280, 20, 160, 20);
            this.uipReconUnit.build  (this.hpReconstruction, 415, 45, 100, 20);
            this.uicbResidualError.build  (this.hpReconstruction, 415, 20, 100, 20);
            
            drawnow
            % Set hsa offset:
            this.hsaInterferogram.setAxesOffset([dTgPx + uitData.Position(1), dTgPy + uitData.Position(2)]);
            
            % refresh button state:
            this.setState(this.u8State);
        end
        
    end
    
end
