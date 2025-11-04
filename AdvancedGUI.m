classdef AdvancedGUI < handle
    % Advanced MATLAB GUI with multiple interactive components
    % Features: Sliders, Text boxes, Notifications, Buttons, Timers, Listeners, Dropdown menus
    
    properties
        % Main figure handle
        fig
        
        % UI Components
        amplitudeSlider
        frequencySlider
        phaseSlider
        amplitudeEdit
        frequencyEdit
        phaseEdit
        startButton
        stopButton
        resetButton
        signalTypeDropdown
        plotTypeDropdown
        notificationPanel
        dataList
        saveButton
        clearButton
        
        % Plot handles
        signalAxes
        spectrumAxes
        
        % Data properties
        time
        signal
        samplingRate = 1000
        duration = 2
        isRunning = false
        
        % Timer
        updateTimer
        
        % Listeners
        sliderListeners
    end
    
    methods
        function obj = AdvancedGUI()
            % Constructor - creates and initializes the GUI
            obj.createGUI();
            obj.setupTimer();
            obj.setupListeners();
            obj.initializeData();
        end
        
        function createGUI(obj)
            % Create main figure
            obj.fig = figure('Name', 'Advanced Signal Generator & Monitor', ...
                'NumberTitle', 'off', ...
                'Position', [100, 100, 1200, 800], ...
                'Resize', 'on', ...
                'Color', [0.95, 0.95, 0.95], ...
                'CloseRequestFcn', @obj.closeGUI);
            
            % Create control panel
            controlPanel = uipanel('Parent', obj.fig, ...
                'Title', 'Signal Controls', ...
                'Position', [0.02, 0.02, 0.25, 0.96], ...
                'BackgroundColor', [0.9, 0.9, 0.9], ...
                'FontSize', 11, ...
                'FontWeight', 'bold');
            
            % Create display panel
            displayPanel = uipanel('Parent', obj.fig, ...
                'Title', 'Signal Display', ...
                'Position', [0.28, 0.52, 0.70, 0.46], ...
                'BackgroundColor', [0.9, 0.9, 0.9]);
            
            % Create spectrum panel
            spectrumPanel = uipanel('Parent', obj.fig, ...
                'Title', 'Frequency Spectrum', ...
                'Position', [0.28, 0.02, 0.70, 0.46], ...
                'BackgroundColor', [0.9, 0.9, 0.9]);
            
            % Create axes for signal plot
            obj.signalAxes = axes('Parent', displayPanel, ...
                'Position', [0.1, 0.15, 0.85, 0.75]);
            title(obj.signalAxes, 'Real-time Signal', 'FontSize', 12, 'FontWeight', 'bold');
            xlabel(obj.signalAxes, 'Time (s)');
            ylabel(obj.signalAxes, 'Amplitude');
            grid(obj.signalAxes, 'on');
            
            % Create axes for spectrum plot
            obj.spectrumAxes = axes('Parent', spectrumPanel, ...
                'Position', [0.1, 0.15, 0.85, 0.75]);
            title(obj.spectrumAxes, 'Frequency Spectrum', 'FontSize', 12, 'FontWeight', 'bold');
            xlabel(obj.spectrumAxes, 'Frequency (Hz)');
            ylabel(obj.spectrumAxes, 'Magnitude');
            grid(obj.spectrumAxes, 'on');
            
            % Create dropdown menus
            uicontrol('Parent', controlPanel, ...
                'Style', 'text', ...
                'String', 'Signal Type:', ...
                'Position', [20, 650, 100, 20], ...
                'BackgroundColor', [0.9, 0.9, 0.9], ...
                'HorizontalAlignment', 'left');
            
            obj.signalTypeDropdown = uicontrol('Parent', controlPanel, ...
                'Style', 'popupmenu', ...
                'String', {'Sine Wave', 'Square Wave', 'Sawtooth Wave', 'Chirp Signal'}, ...
                'Position', [120, 650, 130, 25], ...
                'Callback', @obj.signalTypeChanged, ...
                'BackgroundColor', 'white');
            
            uicontrol('Parent', controlPanel, ...
                'Style', 'text', ...
                'String', 'Plot Type:', ...
                'Position', [20, 610, 100, 20], ...
                'BackgroundColor', [0.9, 0.9, 0.9], ...
                'HorizontalAlignment', 'left');
            
            obj.plotTypeDropdown = uicontrol('Parent', controlPanel, ...
                'Style', 'popupmenu', ...
                'String', {'Line Plot', 'Stem Plot', 'Stairs Plot', 'Area Plot'}, ...
                'Position', [120, 610, 130, 25], ...
                'Callback', @obj.plotTypeChanged, ...
                'BackgroundColor', 'white');
            
            % Create amplitude controls
            uicontrol('Parent', controlPanel, ...
                'Style', 'text', ...
                'String', 'Amplitude:', ...
                'Position', [20, 550, 100, 20], ...
                'BackgroundColor', [0.9, 0.9, 0.9], ...
                'HorizontalAlignment', 'left');
            
            obj.amplitudeSlider = uicontrol('Parent', controlPanel, ...
                'Style', 'slider', ...
                'Min', 0.1, ...
                'Max', 10, ...
                'Value', 1, ...
                'Position', [120, 550, 130, 20], ...
                'BackgroundColor', [0.8, 0.8, 1]);
            
            obj.amplitudeEdit = uicontrol('Parent', controlPanel, ...
                'Style', 'edit', ...
                'String', '1.0', ...
                'Position', [260, 550, 50, 25], ...
                'BackgroundColor', 'white', ...
                'Callback', @obj.amplitudeEditChanged);
            
            % Create frequency controls
            uicontrol('Parent', controlPanel, ...
                'Style', 'text', ...
                'String', 'Frequency (Hz):', ...
                'Position', [20, 500, 100, 20], ...
                'BackgroundColor', [0.9, 0.9, 0.9], ...
                'HorizontalAlignment', 'left');
            
            obj.frequencySlider = uicontrol('Parent', controlPanel, ...
                'Style', 'slider', ...
                'Min', 0.1, ...
                'Max', 50, ...
                'Value', 5, ...
                'Position', [120, 500, 130, 20], ...
                'BackgroundColor', [0.8, 1, 0.8]);
            
            obj.frequencyEdit = uicontrol('Parent', controlPanel, ...
                'Style', 'edit', ...
                'String', '5.0', ...
                'Position', [260, 500, 50, 25], ...
                'BackgroundColor', 'white', ...
                'Callback', @obj.frequencyEditChanged);
            
            % Create phase controls
            uicontrol('Parent', controlPanel, ...
                'Style', 'text', ...
                'String', 'Phase (rad):', ...
                'Position', [20, 450, 100, 20], ...
                'BackgroundColor', [0.9, 0.9, 0.9], ...
                'HorizontalAlignment', 'left');
            
            obj.phaseSlider = uicontrol('Parent', controlPanel, ...
                'Style', 'slider', ...
                'Min', 0, ...
                'Max', 2*pi, ...
                'Value', 0, ...
                'Position', [120, 450, 130, 20], ...
                'BackgroundColor', [1, 0.8, 0.8]);
            
            obj.phaseEdit = uicontrol('Parent', controlPanel, ...
                'Style', 'edit', ...
                'String', '0.0', ...
                'Position', [260, 450, 50, 25], ...
                'BackgroundColor', 'white', ...
                'Callback', @obj.phaseEditChanged);
            
            % Create control buttons
            obj.startButton = uicontrol('Parent', controlPanel, ...
                'Style', 'pushbutton', ...
                'String', 'Start Generation', ...
                'Position', [20, 350, 140, 40], ...
                'BackgroundColor', [0.4, 0.8, 0.4], ...
                'FontWeight', 'bold', ...
                'Callback', @obj.startGeneration);
            
            obj.stopButton = uicontrol('Parent', controlPanel, ...
                'Style', 'pushbutton', ...
                'String', 'Stop Generation', ...
                'Position', [170, 350, 140, 40], ...
                'BackgroundColor', [0.8, 0.4, 0.4], ...
                'FontWeight', 'bold', ...
                'Callback', @obj.stopGeneration, ...
                'Enable', 'off');
            
            obj.resetButton = uicontrol('Parent', controlPanel, ...
                'Style', 'pushbutton', ...
                'String', 'Reset Parameters', ...
                'Position', [20, 300, 290, 30], ...
                'BackgroundColor', [0.8, 0.8, 0.4], ...
                'FontWeight', 'bold', ...
                'Callback', @obj.resetParameters);
            
            % Create data management section
            uicontrol('Parent', controlPanel, ...
                'Style', 'text', ...
                'String', 'Data Management:', ...
                'Position', [20, 250, 150, 20], ...
                'BackgroundColor', [0.9, 0.9, 0.9], ...
                'FontWeight', 'bold', ...
                'HorizontalAlignment', 'left');
            
            obj.saveButton = uicontrol('Parent', controlPanel, ...
                'Style', 'pushbutton', ...
                'String', 'Save Data', ...
                'Position', [20, 210, 140, 30], ...
                'BackgroundColor', [0.4, 0.6, 0.9], ...
                'Callback', @obj.saveData);
            
            obj.clearButton = uicontrol('Parent', controlPanel, ...
                'Style', 'pushbutton', ...
                'String', 'Clear Data', ...
                'Position', [170, 210, 140, 30], ...
                'BackgroundColor', [0.9, 0.6, 0.4], ...
                'Callback', @obj.clearData);
            
            % Create notification panel
            obj.notificationPanel = uipanel('Parent', controlPanel, ...
                'Title', 'Notifications & Status', ...
                'Position', [0.05, 0.02, 0.9, 0.18], ...
                'BackgroundColor', [1, 1, 0.9]);
            
            % Create data listbox
            uicontrol('Parent', obj.notificationPanel, ...
                'Style', 'text', ...
                'String', 'System Messages:', ...
                'Position', [10, 150, 200, 20], ...
                'BackgroundColor', [1, 1, 0.9], ...
                'HorizontalAlignment', 'left');
            
            obj.dataList = uicontrol('Parent', obj.notificationPanel, ...
                'Style', 'listbox', ...
                'String', {'Welcome to Advanced Signal Generator!', 'System initialized successfully.'}, ...
                'Position', [10, 10, 280, 140], ...
                'BackgroundColor', 'white');
            
            % Add status indicators
            obj.createStatusIndicators(controlPanel);
        end
        
        function createStatusIndicators(obj, parent)
            % Create status indicator panel
            statusPanel = uipanel('Parent', parent, ...
                'Title', 'System Status', ...
                'Position', [0.05, 0.21, 0.9, 0.12], ...
                'BackgroundColor', [0.95, 0.95, 0.95]);
            
            % Status indicators
            uicontrol('Parent', statusPanel, ...
                'Style', 'text', ...
                'String', 'Generator:', ...
                'Position', [10, 45, 80, 20], ...
                'BackgroundColor', [0.95, 0.95, 0.95], ...
                'HorizontalAlignment', 'left');
            
            obj.createStatusLight(statusPanel, [100, 45], [0.8, 0.2, 0.2]); % Red - stopped
            
            uicontrol('Parent', statusPanel, ...
                'Style', 'text', ...
                'String', 'Timer:', ...
                'Position', [10, 20, 80, 20], ...
                'BackgroundColor', [0.95, 0.95, 0.95], ...
                'HorizontalAlignment', 'left');
            
            obj.createStatusLight(statusPanel, [100, 20], [0.8, 0.2, 0.2]); % Red - stopped
        end
        
        function createStatusLight(obj, parent, position, color)
            % Create a status light indicator
            uicontrol('Parent', parent, ...
                'Style', 'pushbutton', ...
                'String', '', ...
                'Position', [position(1), position(2), 15, 15], ...
                'BackgroundColor', color, ...
                'Enable', 'off');
        end
        
        function setupTimer(obj)
            % Setup timer for real-time updates
            obj.updateTimer = timer(...
                'ExecutionMode', 'fixedRate', ...
                'Period', 0.1, ... % Update every 100ms
                'TimerFcn', @obj.updateDisplay, ...
                'Name', 'GUIUpdateTimer');
        end
        
        function setupListeners(obj)
            % Setup listeners for slider changes
            obj.sliderListeners = [...
                addlistener(obj.amplitudeSlider, 'ContinuousValueChange', @obj.amplitudeSliderChanged); ...
                addlistener(obj.frequencySlider, 'ContinuousValueChange', @obj.frequencySliderChanged); ...
                addlistener(obj.phaseSlider, 'ContinuousValueChange', @obj.phaseSliderChanged)];
        end
        
        function initializeData(obj)
            % Initialize time and signal data
            obj.time = 0:1/obj.samplingRate:obj.duration;
            obj.signal = zeros(size(obj.time));
            obj.updatePlots();
            obj.addNotification('Data initialized. Ready to start generation.');
        end
        
        function updatePlots(obj)
            % Update both signal and spectrum plots
            cla(obj.signalAxes);
            cla(obj.spectrumAxes);
            
            % Plot signal based on selected type
            plotType = obj.plotTypeDropdown.Value;
            switch plotType
                case 1 % Line plot
                    plot(obj.signalAxes, obj.time, obj.signal, 'b-', 'LineWidth', 2);
                case 2 % Stem plot
                    stem(obj.signalAxes, obj.time, obj.signal, 'filled', 'MarkerSize', 2);
                case 3 % Stairs plot
                    stairs(obj.signalAxes, obj.time, obj.signal, 'LineWidth', 2);
                case 4 % Area plot
                    area(obj.signalAxes, obj.time, obj.signal, 'FaceColor', [0.2, 0.6, 1]);
            end
            
            grid(obj.signalAxes, 'on');
            title(obj.signalAxes, 'Real-time Signal', 'FontSize', 12, 'FontWeight', 'bold');
            xlabel(obj.signalAxes, 'Time (s)');
            ylabel(obj.signalAxes, 'Amplitude');
            
            % Plot spectrum
            if length(obj.signal) > 1
                L = length(obj.signal);
                Y = fft(obj.signal);
                P2 = abs(Y/L);
                P1 = P2(1:floor(L/2)+1);
                P1(2:end-1) = 2*P1(2:end-1);
                f = obj.samplingRate*(0:(L/2))/L;
                
                plot(obj.spectrumAxes, f, P1, 'r-', 'LineWidth', 2);
                grid(obj.spectrumAxes, 'on');
                title(obj.spectrumAxes, 'Frequency Spectrum', 'FontSize', 12, 'FontWeight', 'bold');
                xlabel(obj.spectrumAxes, 'Frequency (Hz)');
                ylabel(obj.spectrumAxes, 'Magnitude');
                xlim(obj.spectrumAxes, [0, 100]);
            end
        end
        
        function generateSignal(obj)
            % Generate signal based on current parameters
            amplitude = obj.amplitudeSlider.Value;
            frequency = obj.frequencySlider.Value;
            phase = obj.phaseSlider.Value;
            signalType = obj.signalTypeDropdown.Value;
            
            switch signalType
                case 1 % Sine wave
                    obj.signal = amplitude * sin(2*pi*frequency*obj.time + phase);
                case 2 % Square wave
                    obj.signal = amplitude * square(2*pi*frequency*obj.time + phase);
                case 3 % Sawtooth wave
                    obj.signal = amplitude * sawtooth(2*pi*frequency*obj.time + phase);
                case 4 % Chirp signal
                    obj.signal = amplitude * chirp(obj.time, 1, obj.duration, frequency*10);
            end
        end
        
        % Callback functions
        function amplitudeSliderChanged(obj, ~, ~)
            value = obj.amplitudeSlider.Value;
            obj.amplitudeEdit.String = sprintf('%.2f', value);
            obj.generateSignal();
            obj.updatePlots();
        end
        
        function frequencySliderChanged(obj, ~, ~)
            value = obj.frequencySlider.Value;
            obj.frequencyEdit.String = sprintf('%.2f', value);
            obj.generateSignal();
            obj.updatePlots();
        end
        
        function phaseSliderChanged(obj, ~, ~)
            value = obj.phaseSlider.Value;
            obj.phaseEdit.String = sprintf('%.2f', value);
            obj.generateSignal();
            obj.updatePlots();
        end
        
        function amplitudeEditChanged(obj, ~, ~)
            value = str2double(obj.amplitudeEdit.String);
            if ~isnan(value) && value >= obj.amplitudeSlider.Min && value <= obj.amplitudeSlider.Max
                obj.amplitudeSlider.Value = value;
                obj.generateSignal();
                obj.updatePlots();
            else
                obj.amplitudeEdit.String = sprintf('%.2f', obj.amplitudeSlider.Value);
            end
        end
        
        function frequencyEditChanged(obj, ~, ~)
            value = str2double(obj.frequencyEdit.String);
            if ~isnan(value) && value >= obj.frequencySlider.Min && value <= obj.frequencySlider.Max
                obj.frequencySlider.Value = value;
                obj.generateSignal();
                obj.updatePlots();
            else
                obj.frequencyEdit.String = sprintf('%.2f', obj.frequencySlider.Value);
            end
        end
        
        function phaseEditChanged(obj, ~, ~)
            value = str2double(obj.phaseEdit.String);
            if ~isnan(value) && value >= obj.phaseSlider.Min && value <= obj.phaseSlider.Max
                obj.phaseSlider.Value = value;
                obj.generateSignal();
                obj.updatePlots();
            else
                obj.phaseEdit.String = sprintf('%.2f', obj.phaseSlider.Value);
            end
        end
        
        function signalTypeChanged(obj, ~, ~)
            obj.generateSignal();
            obj.updatePlots();
            signalTypes = get(obj.signalTypeDropdown, 'String');
            currentType = signalTypes{obj.signalTypeDropdown.Value};
            obj.addNotification(sprintf('Signal type changed to: %s', currentType));
        end
        
        function plotTypeChanged(obj, ~, ~)
            obj.updatePlots();
            plotTypes = get(obj.plotTypeDropdown, 'String');
            currentType = plotTypes{obj.plotTypeDropdown.Value};
            obj.addNotification(sprintf('Plot type changed to: %s', currentType));
        end
        
        function startGeneration(obj, ~, ~)
            if ~obj.isRunning
                obj.isRunning = true;
                start(obj.updateTimer);
                obj.startButton.Enable = 'off';
                obj.stopButton.Enable = 'on';
                obj.addNotification('Signal generation started.');
            end
        end
        
        function stopGeneration(obj, ~, ~)
            if obj.isRunning
                obj.isRunning = false;
                stop(obj.updateTimer);
                obj.startButton.Enable = 'on';
                obj.stopButton.Enable = 'off';
                obj.addNotification('Signal generation stopped.');
            end
        end
        
        function resetParameters(obj, ~, ~)
            % Reset all parameters to default values
            obj.amplitudeSlider.Value = 1;
            obj.frequencySlider.Value = 5;
            obj.phaseSlider.Value = 0;
            obj.amplitudeEdit.String = '1.0';
            obj.frequencyEdit.String = '5.0';
            obj.phaseEdit.String = '0.0';
            obj.signalTypeDropdown.Value = 1;
            obj.plotTypeDropdown.Value = 1;
            
            obj.generateSignal();
            obj.updatePlots();
            obj.addNotification('All parameters reset to default values.');
        end
        
        function saveData(obj, ~, ~)
            % Save current signal data to workspace
            assignin('base', 'generatedSignal', obj.signal);
            assignin('base', 'signalTime', obj.time);
            assignin('base', 'signalParameters', struct(...
                'amplitude', obj.amplitudeSlider.Value, ...
                'frequency', obj.frequencySlider.Value, ...
                'phase', obj.phaseSlider.Value));
            
            obj.addNotification('Signal data saved to workspace.');
        end
        
        function clearData(obj, ~, ~)
            % Clear the notification list
            obj.dataList.String = {};
            obj.addNotification('Notification list cleared.');
        end
        
        function updateDisplay(obj, ~, ~)
            % Update display with new data (called by timer)
            if obj.isRunning
                % Add some random variation to simulate real-time data
                variation = 0.1 * randn(size(obj.signal));
                tempSignal = obj.signal + variation;
                
                % Update plot with new data
                plot(obj.signalAxes, obj.time, tempSignal, 'b-', 'LineWidth', 1.5);
                grid(obj.signalAxes, 'on');
                title(obj.signalAxes, 'Real-time Signal (Live)', 'FontSize', 12, 'FontWeight', 'bold');
                xlabel(obj.signalAxes, 'Time (s)');
                ylabel(obj.signalAxes, 'Amplitude');
                
                drawnow;
            end
        end
        
        function addNotification(obj, message)
            % Add a notification message to the list
            timestamp = datestr(now, 'HH:MM:SS');
            fullMessage = sprintf('[%s] %s', timestamp, message);
            currentList = obj.dataList.String;
            if isempty(currentList)
                currentList = {fullMessage};
            else
                currentList = [currentList; {fullMessage}];
            end
            obj.dataList.String = currentList;
            obj.dataList.Value = length(currentList); % Auto-scroll to bottom
        end
        
        function closeGUI(obj, ~, ~)
            % Clean up when closing the GUI
            if obj.isRunning
                stop(obj.updateTimer);
            end
            delete(obj.updateTimer);
            delete(obj.fig);
        end
    end
end

% Helper function to create and display the GUI
function launchAdvancedGUI()
    % Launch the advanced GUI
    gui = AdvancedGUI();
end

