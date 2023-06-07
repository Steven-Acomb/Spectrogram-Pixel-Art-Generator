classdef AdditiveSynthPatch
    properties
        partials
        xRng
        yRng
        fs
        td
        logf
        fMax
        fMin
        dBpk
        pp
        y
        tt
        ff
        aa
    end 
    methods

        % Class Constructor
        function obj = AdditiveSynthPatch(timbrePartials,xRng,yRng, ...
                fs,td,logf,fMax,fMin,dBpk)
            if nargin == 0
                obj.partials = [];
                obj.xRng = [0 10];
                obj.yRng = [0 10];
                obj.fs = 4e4;
                obj.td = 1;
                obj.logf = false;
                obj.fMax = 0.25*obj.fs;
                obj.fMax = 0;
                obj.dBpk = -1;
            else
                obj.partials = timbrePartials;
                obj.xRng = xRng;
                obj.yRng = yRng;
                obj.fs = fs;
                obj.td = td;
                obj.logf = logf;
                obj.fMax = fMax;
                obj.fMin = fMin;
                obj.dBpk = dBpk;
            end

            % Construct frequency and amplitude trajectories
                tScale = obj.td/diff(obj.xRng);
                obj.tt = linspace(obj.xRng(1),obj.xRng(2)-1/obj.fs,obj.td*obj.fs)'*tScale;
                obj.ff = zeros(length(obj.tt),length(obj.partials));
                obj.aa = obj.ff;
                for k = 1:length(obj.partials)
                    if obj.logf
                        obj.ff(:,k) = interp1(obj.partials(k).fx*tScale, ...
                        yRef_Hz*2.^(obj.partials(k).fy*yInt), obj.tt, obj.partials(k).fim);
                    else
                        obj.ff(:,k) = interp1(obj.partials(k).fx*tScale, ...
                        ((obj.partials(k).fy-obj.yRng(1))/diff(obj.yRng)) * ...
                        (obj.fMax-obj.fMin)+obj.fMin,obj.tt,obj.partials(k).fim);
                    end
                    % Convert to linear amplitude using A = 2^(AdB/6)
                    obj.aa(:,k) = 2.^(interp1(obj.partials(k).ax*tScale,obj.partials(k).ay,obj.tt,obj.partials(k).aim)/6);
                    %     obj.aa(isnan(obj.ff(:,k)),k) = 0;
                    %     obj.ff(isnan(obj.ff(:,k)),k) = obj.fMin;
                end

            % Look for NaNs in frequency or amplitude trajectory, mute those in
            % amplitude trajectory, and then replace them with DC in the frequency
            % trajectory
            obj.aa(isnan(obj.ff) | isnan(obj.aa)) = 0;
            obj.ff(isnan(obj.ff)) = 0;

            % Create the audio partials, sum to output, and normalize level
            % Integrate frequency trajectory to obtain phase trajectory for sine
            % function and scale by amplitude trajectory
            obj.pp = obj.aa .* sin(2*pi*cumsum(obj.ff)/obj.fs);

            % Sum the partials and normalize level to desired peak level
            obj.y = sum(obj.pp,2);
            obj.y = 2^(obj.dBpk/6)*obj.y/max(abs(obj.y(:)));
        end

        % Plot the frequency & amplitude trajectories and the spectrogram
        function show(obj,specFloor_dB,specColormap)
            tiledlayout(1,3)
            xx = linspace(obj.xRng(1),obj.xRng(2),500);
            for k = 1:length(obj.partials)
                nexttile(1)
                plot(xx,interp1(obj.partials(k).fx,obj.partials(k).fy,xx,obj.partials(k).fim));
                grid on, axis padded, hold on
                plot(obj.partials(k).fx,obj.partials(k).fy,'o')
                text(-0.05*diff(obj.xRng)+obj.partials(k).fx(1),obj.partials(k).fy(1),num2str(k))
                nexttile(2)
                plot(xx,interp1(obj.partials(k).ax,obj.partials(k).ay,xx,obj.partials(k).aim));
                grid on, axis padded, hold on
                plot(obj.partials(k).ax,obj.partials(k).ay,'o')
                text(-0.05*diff(obj.xRng)+obj.partials(k).ax(1),obj.partials(k).ay(1),num2str(k))
            end
            nexttile(1), hold off
            title('Frequency Trajectories'), xlabel('x'), ylabel('y')
            nexttile(2), hold off
            title('Amplitude Trajectories'), xlabel('x'), ylabel('y (dB)')
            
            nexttile
            winLength = floor(length(obj.y)/64);
            spectrogram(obj.y,kaiser(winLength,5),winLength-1,[], ... 
                obj.fs,"yaxis", MinThreshold = specFloor_dB)
            colormap(specColormap)
        end

        % Plot the frequency trajectories
        function show_freq_traj(obj)
            xx = linspace(obj.xRng(1),obj.xRng(2),500);
            for k = 1:length(obj.partials)
                plot(xx,interp1(obj.partials(k).fx,obj.partials(k).fy,xx,obj.partials(k).fim));
                grid on, axis padded, hold on
                plot(obj.partials(k).fx,obj.partials(k).fy,'o')
                text(-0.05*diff(obj.xRng)+obj.partials(k).fx(1),obj.partials(k).fy(1),num2str(k))
            end
            hold off
            title('Frequency Trajectories'), xlabel('x'), ylabel('y')
        end

        % Plot the amplitude trajectories
        function show_ampl_traj(obj)
            xx = linspace(obj.xRng(1),obj.xRng(2),500);
            for k = 1:length(obj.partials)
                plot(xx,interp1(obj.partials(k).ax,obj.partials(k).ay,xx,obj.partials(k).aim));
                grid on, axis padded, hold on
                plot(obj.partials(k).ax,obj.partials(k).ay,'o')
                text(-0.05*diff(obj.xRng)+obj.partials(k).ax(1),obj.partials(k).ay(1),num2str(k))
            end
            hold off
            title('Amplitude Trajectories'), xlabel('x'), ylabel('y (dB)')
        end

        % Plot the spectrogram
        function show_spec(obj,specFloor_dB,specColormap)
            winLength = floor(length(obj.y)/64);
%             spectrogram(obj.y,kaiser(winLength,5),winLength-1,[], ... 
%                 obj.fs,"yaxis", MinThreshold = specFloor_dB)
            spectrogram(obj.y,kaiser(winLength,5),0,[], ... 
                obj.fs,"yaxis", MinThreshold = specFloor_dB)
            colormap(specColormap)
        end

        % Listen to each partial individually or the complete sound
        function listen(obj,specOneAtATime)
            if (specOneAtATime)
%                 disp('Individual partials...')
                for k = 1:size(obj.pp,2)
                    fprintf('Partial #%g\n',k);
                    sound(obj.pp(:,k),obj.fs)
                    pause(obj.td+1)
                end
            else
%                 disp('Complete sound')
                sound(obj.y,obj.fs)
            end
            
        end

        % Write audio to a file
        function write(obj,audioFilename)
            audiowrite(audioFilename,obj.y,obj.fs);
        end
    end
end