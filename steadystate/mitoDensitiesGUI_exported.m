classdef mitoDensitiesGUI_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                       matlab.ui.Figure
        GridLayout                     matlab.ui.container.GridLayout
        LeftPanel                      matlab.ui.container.Panel
        ReactivationratekaSlider       matlab.ui.control.Slider
        ReactivationratekaSliderLabel  matlab.ui.control.Label
        FissionratekbeachsideSlider    matlab.ui.control.Slider
        FissionratekbeachsideLabel     matlab.ui.control.Label
        FusionprobscalinggammaSlider   matlab.ui.control.Slider
        FusionprobscalinggammaSliderLabel  matlab.ui.control.Label
        FusionandFissionLabel          matlab.ui.control.Label
        ProximalfuseprobSlider         matlab.ui.control.Slider
        ProximalfuseprobSliderLabel    matlab.ui.control.Label
        PauseandRestartLabel           matlab.ui.control.Label
        RestartrateksSlider            matlab.ui.control.Slider
        RestartrateksSliderLabel       matlab.ui.control.Label
        PauseratescalinggammasSlider   matlab.ui.control.Slider
        PauseratescalinggammasSliderLabel  matlab.ui.control.Label
        ProximalpauserateSlider        matlab.ui.control.Slider
        ProximalpauserateSliderLabel   matlab.ui.control.Label
        VelocityvSlider                matlab.ui.control.Slider
        VelocityvSliderLabel           matlab.ui.control.Label
        ProductionratekpSlider         matlab.ui.control.Slider
        ProductionratekpSliderLabel    matlab.ui.control.Label
        RightPanel                     matlab.ui.container.Panel
        trunkradiusEditField           matlab.ui.control.NumericEditField
        trunkradiusEditFieldLabel      matlab.ui.control.Label
        BranchradiiButtonGroup         matlab.ui.container.ButtonGroup
        SetbalancedradiiButton         matlab.ui.control.RadioButton
        UsepresavededgevalsButton      matlab.ui.control.RadioButton
        AvgclustersizeLabel            matlab.ui.control.Label
        MitovolumefractionLabel        matlab.ui.control.Label
        DistalenrichmentLabel          matlab.ui.control.Label
        UIAxes                         matlab.ui.control.UIAxes
    end

    % Properties that correspond to apps with auto-reflow
    properties (Access = private)
        onePanelWidth = 576;
    end

    
    properties (Access = public)
        NT % Network object
        trunkedge % which edge is the trunk
        radii % radius for each network edge
        patchH % handles to patch objects
        param % parameters describing mito dynamics
        distinfo % information on mito distributions
    end
    
    methods (Access = private)
        
        function results = replotNetwork(app, fieldVals)
            % Completely redraw the network
            % using colors in fieldVals

            if isempty(app.NT)
                disp('No network set up.')
                results = NaN;
                return
            end

            % Clear axes and set it as the active target
            cla(app.UIAxes);
            plotopt = struct('Parent',app.UIAxes);
            app.patchH = app.NT.plotNetworkField(app.radii,fieldVals,plotopt);

            % turn off ticks, turn on colorbar
            set(app.UIAxes,'XTick',[],'YTick',[])
            colorbar(app.UIAxes)
        end
        
        function param = setDefaultParam(app)
            % Set up default parameters at the start

            param = struct();
            % anterograde mito flux: # units entering trunk per second
            param.kp=1/60;
            
            
            % velocity of anterograde and retrograde moving mitos (um/sec)
            param.v=0.5;
            
            % Pausing occurs at rate: As/radius^gammas (units of per sec)
            param.gammas = 0.8;
            param.As=0.1;
            
            % rate of a paused mitochondrion restarting and becoming mobile again
            param.kr = 0.003;
            
            % Rate of fission (from each side of a cluster); units of per sec
            param.kb = 0.001;
            % Probability of fusion when two mitos pass: Au/radius^gamma
            param.gamma = 2;
            param.Au = 0.1;
            
            % cut off fusion probabilities, so they cannot go above 1
            param.probmax = 1;
            
            % reactivation rate (in both directions together)
            % for a unit-size cluster left behind after fission
            param.ka = 2*param.kb;
            
            % Below are some definitions just for calculating metrics
            % distcutoff defines the distal region of the tree (for calculatin distal
            % enrichment)
            % graph distance to trunk must be more than this fraction of the max
            % graph distance
            param.distcutoff = 0.7;
            
            % volume of a mitochondrial unit (um^3)
            % (for calculating volume fraction, distal enrichment)
            param.mitovol = 0.5;

            app.param = param;
        end
        
        function updateDistrib(app)
            % update mito distribution info

            app.distinfo = getMitoDensitiesMetrics(app.NT,app.radii,app.param);
        end

        function updateNetworkColor(app, fieldVals)
            % update the colors on the network plot according to the
            % provided field values
            % if fieldVals not provided, use current distribution volume
            % densities

            if (nargin<2)
                fieldVals = app.distinfo.M1vals*app.param.mitovol./(pi*app.radii.^2);
            end

            for ec = 1:app.NT.nedge
                app.patchH(ec).CData = fieldVals(ec);
            end
        end

        function paramChangeUpdate(app)
            % run this function to recalculate and replot whenever the
            % parameters change

            app.updateDistrib()
            app.updateNetworkColor()

            app.MitovolumefractionLabel.Text = sprintf('Mito volume fraction: %g', app.distinfo.volfrac);
            app.AvgclustersizeLabel.Text = sprintf('Average cluster size: %g', app.distinfo.avgclustsize);
            app.DistalenrichmentLabel.Text = sprintf('Distal enrichment: %g', app.distinfo.distenrich);
        end
        
        function setupSliders(app)
            % set limits and logarithmic spacing as needed for sliders            

            % production rate
            ticks = [0.001,0.01,0.1,1];
            app.ProductionratekpSlider.Limits = log10([min(ticks) max(ticks)])  ;          
            app.ProductionratekpSlider.MajorTicks = log10(ticks);
            app.ProductionratekpSlider.MajorTickLabels = cellstr(string(ticks));
            app.ProductionratekpSlider.MinorTicks = log10([0.001:0.001:0.009, 0.02:0.01:0.09, 0.2:0.1:0.9]);
            
            app.ProductionratekpSlider.Value = log10(app.param.kp);
            app.VelocityvSlider.Value = app.param.v;

            % --------- Pausing and restarting --------------
            % proximal stop rate
            ticks = [0.01,0.1,1];
            app.ProximalpauserateSlider.Limits = log10([min(ticks) max(ticks)])  ;          
            app.ProximalpauserateSlider.MajorTicks = log10(ticks);
            app.ProximalpauserateSlider.MajorTickLabels = cellstr(string(ticks));
            app.ProximalpauserateSlider.MinorTicks = log10([0.01:0.01:0.09, 0.1:0.1:0.9]);
            app.ProximalpauserateSlider.Value = log10(app.param.As/app.radii(app.trunkedge)^app.param.gammas);

            % scaling
            app.PauseratescalinggammasSlider.Value = app.param.gammas;

            % restart rate
            ticks = [0.001,0.01,0.1,1];
            app.RestartrateksSlider.Limits = log10([min(ticks) max(ticks)])  ;          
            app.RestartrateksSlider.MajorTicks = log10(ticks);
            app.RestartrateksSlider.MajorTickLabels = cellstr(string(ticks));
            app.RestartrateksSlider.MinorTicks = log10([0.001:0.001:0.009, 0.01:0.01:0.09, 0.1:0.1:0.9]);
            app.RestartrateksSlider.Value = log10(app.param.kr);

            % -------------- Fusion and fission --------------
            % proximal fusion probability
            ticks = [0.001,0.01,0.1,1];
            app.ProximalfuseprobSlider.Limits = log10([min(ticks) max(ticks)])  ;          
            app.ProximalfuseprobSlider.MajorTicks = log10(ticks);
            app.ProximalfuseprobSlider.MajorTickLabels = cellstr(string(ticks));
            app.ProximalfuseprobSlider.MinorTicks = log10([0.001:0.001:0.009, 0.01:0.01:0.09, 0.1:0.1:0.9]);
            app.ProximalfuseprobSlider.Value = log10(app.param.Au/app.radii(app.trunkedge)^app.param.gamma);

            app.FusionprobscalinggammaSlider.Value = app.param.gamma;

            % fission rate (each side)
            ticks = [0.0001,0.001,0.01,0.1];
            app.FissionratekbeachsideSlider.Limits = log10([min(ticks) max(ticks)])  ;          
            app.FissionratekbeachsideSlider.MajorTicks = log10(ticks);
            app.FissionratekbeachsideSlider.MajorTickLabels = cellstr(string(ticks));
            app.FissionratekbeachsideSlider.MinorTicks = log10([1e-4:1e-4:9e-4, 0.001:0.001:0.009, 0.01:0.01:0.09]);
            app.FissionratekbeachsideSlider.Value = log10(app.param.kb);

             % reactivation rate (both sides together)
            ticks = [0.0001,0.001,0.01,0.1];
            app.ReactivationratekaSlider.Limits = log10([min(ticks) max(ticks)])  ;          
            app.ReactivationratekaSlider.MajorTicks = log10(ticks);
            app.ReactivationratekaSlider.MajorTickLabels = cellstr(string(ticks));
            app.ReactivationratekaSlider.MinorTicks = log10([1e-4:1e-4:9e-4, 0.001:0.001:0.009, 0.01:0.01:0.09]);
            app.ReactivationratekaSlider.Value = log10(app.param.ka);

        end
        
        function setRadii(app)
            % set the branch radii field, according to the GUI radio
            % buttons

            if (app.SetbalancedradiiButton.Value)
                % Set balanced radii with whatever trunk radius is entered
                % in the textbox

                % Da Vinci Law: r1^2 + r2^2 = r0^2
                a = 2;
                rm = 0; % minimal allowed radius (in um)

                % calculate length, depth for all subtrees
                [stL,~,stD,~] = setSubtreeInfo(app.NT,app.trunkedge,a,'L/D');

                rtrunk = app.trunkradiusEditField.Value;

                % Split branch areas in proportion to bushiness: r1^2/r2^2 = (L1/D1)/(L2/D2)    
                app.radii= setRadiiWithRm(app.NT,app.trunkedge,2,rm,rtrunk,stL./stD);

            elseif (app.UsepresavededgevalsButton.Value)
                % assume edgevals contains radii^2
                app.radii = sqrt(app.NT.edgevals(:,2))';

                % update the trunk radius textbox
                app.trunkradiusEditField.Value = app.radii(app.trunkedge);
            else
                error('this should never happen')
            end
        end
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app, varargin)
            % set defaults
            app.NT = [];
            app.radii = [];
            app.patchH = [];
            
            % default parameters
            app.setDefaultParam()                        
          
            % optional arguments
            % Set up a network if one was passed to the GUI
            inputradii = false;
            for index = 1:2:length(varargin)
                switch(lower(varargin{index}))
                    case 'nt'
                        app.NT = varargin{index+1};   
                    case 'radii'
                        app.radii = varargin{index+1};
                        inputradii = true; % input radii are provided
                end          
            end

            if (~isempty(app.NT))

                % set trunkedge
                app.trunkedge = app.NT.nodeedges(app.NT.rootnode,1);

                if (~inputradii)
                    % calculate radii according to the radio buttons
                    app.setRadii()
                end
                
                % Compute mitochondrial distribution
                app.updateDistrib()

                % get volume density on each branch
                voldens = app.distinfo.M1vals*app.param.mitovol./(pi*app.radii.^2);

                % plot with identical field values                
                app.replotNetwork(voldens);

                % also set all the metric labels
                app.paramChangeUpdate();
            end

            % set up the slider values and limits
            app.setupSliders()


            % try to fix axes appearance
            app.UIAxes.PositionConstraint = 'innerposition';
        end

        % Value changed function: FissionratekbeachsideSlider, 
        % ...and 8 other components
        function paramsFromSliders(app, event)
                       
            app.param.v = app.VelocityvSlider.Value;
            app.param.kp = 10.^app.ProductionratekpSlider.Value;            
            app.param.gammas = app.PauseratescalinggammasSlider.Value;
            proxstop = 10.^app.ProximalpauserateSlider.Value;
            app.param.As = proxstop*app.radii(app.trunkedge)^app.param.gammas;
            app.param.kr = 10^app.RestartrateksSlider.Value;
            
            app.param.gamma = app.FusionprobscalinggammaSlider.Value;
            Pu0 = 10.^app.ProximalfuseprobSlider.Value;
            app.param.Au = Pu0*app.radii(app.trunkedge)^app.param.gamma;
            app.param.kb = 10^app.FissionratekbeachsideSlider.Value;
            app.param.ka = 10^app.ReactivationratekaSlider.Value;
            
            app.paramChangeUpdate()
        end

        % Callback function: BranchradiiButtonGroup, trunkradiusEditField
        function radiiChangeUpdate(app, event)
            selectedButton = app.BranchradiiButtonGroup.SelectedObject;
            
            fieldVals= zeros(1,app.NT.nedge);            
            app.setRadii();
            app.replotNetwork(fieldVals)
            app.paramChangeUpdate();            

            
        end

        % Changes arrangement of the app based on UIFigure width
        function updateAppLayout(app, event)
            currentFigureWidth = app.UIFigure.Position(3);
            if(currentFigureWidth <= app.onePanelWidth)
                % Change to a 2x1 grid
                app.GridLayout.RowHeight = {651, 651};
                app.GridLayout.ColumnWidth = {'1x'};
                app.RightPanel.Layout.Row = 2;
                app.RightPanel.Layout.Column = 1;
            else
                % Change to a 1x2 grid
                app.GridLayout.RowHeight = {'1x'};
                app.GridLayout.ColumnWidth = {367, '1x'};
                app.RightPanel.Layout.Row = 1;
                app.RightPanel.Layout.Column = 2;
            end
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.AutoResizeChildren = 'off';
            app.UIFigure.Position = [100 100 1014 651];
            app.UIFigure.Name = 'MATLAB App';
            app.UIFigure.SizeChangedFcn = createCallbackFcn(app, @updateAppLayout, true);

            % Create GridLayout
            app.GridLayout = uigridlayout(app.UIFigure);
            app.GridLayout.ColumnWidth = {367, '1x'};
            app.GridLayout.RowHeight = {'1x'};
            app.GridLayout.ColumnSpacing = 0;
            app.GridLayout.RowSpacing = 0;
            app.GridLayout.Padding = [0 0 0 0];
            app.GridLayout.Scrollable = 'on';

            % Create LeftPanel
            app.LeftPanel = uipanel(app.GridLayout);
            app.LeftPanel.Layout.Row = 1;
            app.LeftPanel.Layout.Column = 1;

            % Create ProductionratekpSliderLabel
            app.ProductionratekpSliderLabel = uilabel(app.LeftPanel);
            app.ProductionratekpSliderLabel.HorizontalAlignment = 'right';
            app.ProductionratekpSliderLabel.Position = [27 579 114 22];
            app.ProductionratekpSliderLabel.Text = 'Production rate (kp):';

            % Create ProductionratekpSlider
            app.ProductionratekpSlider = uislider(app.LeftPanel);
            app.ProductionratekpSlider.Limits = [0 1];
            app.ProductionratekpSlider.ValueChangedFcn = createCallbackFcn(app, @paramsFromSliders, true);
            app.ProductionratekpSlider.Position = [162 588 150 3];
            app.ProductionratekpSlider.Value = 0.01;

            % Create VelocityvSliderLabel
            app.VelocityvSliderLabel = uilabel(app.LeftPanel);
            app.VelocityvSliderLabel.HorizontalAlignment = 'right';
            app.VelocityvSliderLabel.Position = [45 510 67 22];
            app.VelocityvSliderLabel.Text = 'Velocity (v):';

            % Create VelocityvSlider
            app.VelocityvSlider = uislider(app.LeftPanel);
            app.VelocityvSlider.Limits = [0 1];
            app.VelocityvSlider.ValueChangedFcn = createCallbackFcn(app, @paramsFromSliders, true);
            app.VelocityvSlider.Position = [173 520 136 3];
            app.VelocityvSlider.Value = 0.5;

            % Create ProximalpauserateSliderLabel
            app.ProximalpauserateSliderLabel = uilabel(app.LeftPanel);
            app.ProximalpauserateSliderLabel.HorizontalAlignment = 'right';
            app.ProximalpauserateSliderLabel.Position = [28 417 118 22];
            app.ProximalpauserateSliderLabel.Text = 'Proximal pause rate: ';

            % Create ProximalpauserateSlider
            app.ProximalpauserateSlider = uislider(app.LeftPanel);
            app.ProximalpauserateSlider.ValueChangedFcn = createCallbackFcn(app, @paramsFromSliders, true);
            app.ProximalpauserateSlider.Position = [162 427 144 3];

            % Create PauseratescalinggammasSliderLabel
            app.PauseratescalinggammasSliderLabel = uilabel(app.LeftPanel);
            app.PauseratescalinggammasSliderLabel.HorizontalAlignment = 'right';
            app.PauseratescalinggammasSliderLabel.Position = [26 355 168 22];
            app.PauseratescalinggammasSliderLabel.Text = 'Pause rate scaling (gammas): ';

            % Create PauseratescalinggammasSlider
            app.PauseratescalinggammasSlider = uislider(app.LeftPanel);
            app.PauseratescalinggammasSlider.Limits = [0 2];
            app.PauseratescalinggammasSlider.ValueChangedFcn = createCallbackFcn(app, @paramsFromSliders, true);
            app.PauseratescalinggammasSlider.Position = [213 364 93 3];

            % Create RestartrateksSliderLabel
            app.RestartrateksSliderLabel = uilabel(app.LeftPanel);
            app.RestartrateksSliderLabel.HorizontalAlignment = 'right';
            app.RestartrateksSliderLabel.Position = [27 290 94 22];
            app.RestartrateksSliderLabel.Text = 'Restart rate (ks):';

            % Create RestartrateksSlider
            app.RestartrateksSlider = uislider(app.LeftPanel);
            app.RestartrateksSlider.ValueChangedFcn = createCallbackFcn(app, @paramsFromSliders, true);
            app.RestartrateksSlider.Position = [162 300 139 3];

            % Create PauseandRestartLabel
            app.PauseandRestartLabel = uilabel(app.LeftPanel);
            app.PauseandRestartLabel.FontWeight = 'bold';
            app.PauseandRestartLabel.Position = [130 438 110 30];
            app.PauseandRestartLabel.Text = 'Pause and Restart';

            % Create ProximalfuseprobSliderLabel
            app.ProximalfuseprobSliderLabel = uilabel(app.LeftPanel);
            app.ProximalfuseprobSliderLabel.HorizontalAlignment = 'right';
            app.ProximalfuseprobSliderLabel.Position = [28 194 108 22];
            app.ProximalfuseprobSliderLabel.Text = 'Proximal fuse prob:';

            % Create ProximalfuseprobSlider
            app.ProximalfuseprobSlider = uislider(app.LeftPanel);
            app.ProximalfuseprobSlider.ValueChangedFcn = createCallbackFcn(app, @paramsFromSliders, true);
            app.ProximalfuseprobSlider.Position = [162 204 139 3];

            % Create FusionandFissionLabel
            app.FusionandFissionLabel = uilabel(app.LeftPanel);
            app.FusionandFissionLabel.FontWeight = 'bold';
            app.FusionandFissionLabel.Position = [131 215 114 29];
            app.FusionandFissionLabel.Text = 'Fusion and Fission';

            % Create FusionprobscalinggammaSliderLabel
            app.FusionprobscalinggammaSliderLabel = uilabel(app.LeftPanel);
            app.FusionprobscalinggammaSliderLabel.HorizontalAlignment = 'right';
            app.FusionprobscalinggammaSliderLabel.Position = [22 140 164 22];
            app.FusionprobscalinggammaSliderLabel.Text = 'Fusion prob scaling (gamma):';

            % Create FusionprobscalinggammaSlider
            app.FusionprobscalinggammaSlider = uislider(app.LeftPanel);
            app.FusionprobscalinggammaSlider.Limits = [0 2];
            app.FusionprobscalinggammaSlider.ValueChangedFcn = createCallbackFcn(app, @paramsFromSliders, true);
            app.FusionprobscalinggammaSlider.Position = [220 150 93 3];

            % Create FissionratekbeachsideLabel
            app.FissionratekbeachsideLabel = uilabel(app.LeftPanel);
            app.FissionratekbeachsideLabel.HorizontalAlignment = 'right';
            app.FissionratekbeachsideLabel.Position = [41 80 94 30];
            app.FissionratekbeachsideLabel.Text = {'Fission rate (kb):'; '(each side)'};

            % Create FissionratekbeachsideSlider
            app.FissionratekbeachsideSlider = uislider(app.LeftPanel);
            app.FissionratekbeachsideSlider.ValueChangedFcn = createCallbackFcn(app, @paramsFromSliders, true);
            app.FissionratekbeachsideSlider.Position = [170 97 140 3];

            % Create ReactivationratekaSliderLabel
            app.ReactivationratekaSliderLabel = uilabel(app.LeftPanel);
            app.ReactivationratekaSliderLabel.HorizontalAlignment = 'right';
            app.ReactivationratekaSliderLabel.Position = [29 36 126 22];
            app.ReactivationratekaSliderLabel.Text = 'Reactivation rate (ka): ';

            % Create ReactivationratekaSlider
            app.ReactivationratekaSlider = uislider(app.LeftPanel);
            app.ReactivationratekaSlider.ValueChangedFcn = createCallbackFcn(app, @paramsFromSliders, true);
            app.ReactivationratekaSlider.Position = [171 45 140 3];

            % Create RightPanel
            app.RightPanel = uipanel(app.GridLayout);
            app.RightPanel.TitlePosition = 'centertop';
            app.RightPanel.Layout.Row = 1;
            app.RightPanel.Layout.Column = 2;

            % Create UIAxes
            app.UIAxes = uiaxes(app.RightPanel);
            title(app.UIAxes, 'Mito volume fraction on each branch')
            app.UIAxes.Position = [178 106 381 385];

            % Create DistalenrichmentLabel
            app.DistalenrichmentLabel = uilabel(app.RightPanel);
            app.DistalenrichmentLabel.FontSize = 14;
            app.DistalenrichmentLabel.Position = [188 514 252 34];
            app.DistalenrichmentLabel.Text = 'Distal enrichment:';

            % Create MitovolumefractionLabel
            app.MitovolumefractionLabel = uilabel(app.RightPanel);
            app.MitovolumefractionLabel.FontSize = 14;
            app.MitovolumefractionLabel.Position = [188 590 252 34];
            app.MitovolumefractionLabel.Text = 'Mito volume fraction:';

            % Create AvgclustersizeLabel
            app.AvgclustersizeLabel = uilabel(app.RightPanel);
            app.AvgclustersizeLabel.FontSize = 14;
            app.AvgclustersizeLabel.Position = [188 549 235 38];
            app.AvgclustersizeLabel.Text = 'Avg cluster size: ';

            % Create BranchradiiButtonGroup
            app.BranchradiiButtonGroup = uibuttongroup(app.RightPanel);
            app.BranchradiiButtonGroup.SelectionChangedFcn = createCallbackFcn(app, @radiiChangeUpdate, true);
            app.BranchradiiButtonGroup.Title = 'Branch radii:';
            app.BranchradiiButtonGroup.Position = [202 14 175 83];

            % Create UsepresavededgevalsButton
            app.UsepresavededgevalsButton = uiradiobutton(app.BranchradiiButtonGroup);
            app.UsepresavededgevalsButton.Text = 'Use presaved edgevals';
            app.UsepresavededgevalsButton.Position = [17 30 147 22];

            % Create SetbalancedradiiButton
            app.SetbalancedradiiButton = uiradiobutton(app.BranchradiiButtonGroup);
            app.SetbalancedradiiButton.Text = 'Set balanced radii';
            app.SetbalancedradiiButton.Position = [17 8 118 22];
            app.SetbalancedradiiButton.Value = true;

            % Create trunkradiusEditFieldLabel
            app.trunkradiusEditFieldLabel = uilabel(app.RightPanel);
            app.trunkradiusEditFieldLabel.HorizontalAlignment = 'right';
            app.trunkradiusEditFieldLabel.Position = [386 19 71 22];
            app.trunkradiusEditFieldLabel.Text = 'trunk radius:';

            % Create trunkradiusEditField
            app.trunkradiusEditField = uieditfield(app.RightPanel, 'numeric');
            app.trunkradiusEditField.ValueChangedFcn = createCallbackFcn(app, @radiiChangeUpdate, true);
            app.trunkradiusEditField.Position = [472 17 35 25];
            app.trunkradiusEditField.Value = 2.7;

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = mitoDensitiesGUI_exported(varargin)

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            % Execute the startup function
            runStartupFcn(app, @(app)startupFcn(app, varargin{:}))

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end