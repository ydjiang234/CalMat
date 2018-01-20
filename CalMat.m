classdef CalMat
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        E; K; A; I; L; N; backbone; targetX; targetY; targetEnergy;%input parameters, to descripe the target member and target data
        K0; K0amp; as, thetap,thetapc; thetay; Fy; thetac; Fc; thetau; Res; backboneShifted;%describe the backbone of the CP model
        c_S = 1.0; c_C = 1.0; c_A = 1.0; c_K = 1.0;
        D = 1.0;
        revK; revKamp;%describe the additional reversed material
        ampFactor;
        turning; DispList; Energy; BBcyclic; targetBB;
        template;
        considerNum = 10;
        d_incr;
        working_path = 'Working';
    end
    
    methods
        function obj = CalMat(A, I, L, N, revK, backbone, targetData, ampFactor, d_incr)
            obj.A = A;
            obj.I = I;
            obj.L = L;
            obj.N = N;
            obj.revK =revK;
            obj.backbone = backbone;
            obj.targetX = targetData(:,1);
            obj.targetY = targetData(:,2);
            obj.ampFactor = ampFactor;
            obj.d_incr = d_incr;
            obj = obj.initialize();
        end
        
        function obj = initialize(obj)
            %stiffness of the elastic member
            obj.K = obj.backbone(1,2) / obj.backbone(1,1);
            obj.E = obj.K * obj.L / obj.I / 3.0;
            %To get the property of the CP backbone curve
            obj = obj.shiftBackbone();
            obj.thetay = obj.backboneShifted(1,1);
            obj.Fy = obj.backboneShifted(1,2);
            obj.thetac = obj.backboneShifted(2,1);
            obj.Fc = obj.backboneShifted(2,2);
            obj.thetau = obj.backbone(3,1);%Be careful to use backbone instead of backboneShifted
            obj.K0 = obj.Fy / obj.thetay;
            obj.K0amp = obj.K0 * (obj.ampFactor + 1.0);
            obj.thetap = obj.thetac - obj.thetay;
            obj.thetapc = obj.backboneShifted(3,1) - obj.thetac;
            obj.as = (obj.Fc - obj.Fy) / obj.thetap / obj.K0;
            %for the reverse material
            obj.revKamp = obj.revK;
            %other
            obj.turning = obj.findTurning(horzcat(obj.targetX, obj.targetY));
            obj = obj.convertDisp();
            obj = obj.preRender();
            
            %Initial targetEnergy
            %obj.Energy = obj.getEnergy(obj.targetX, obj.targetY);
            %obj.Energy = obj.monoData(obj.Energy);
            %xRange = linspace(obj.Energy(1,1), obj.Energy(length(obj.Energy),1), obj.considerNum);
            %obj.targetEnergy = obj.simplifiedEnergy(obj.Energy, xRange);
            
            %Initial backbone
            obj.BBcyclic = obj.monoData(obj.turning(:,2:3));
            xRange = linspace(obj.BBcyclic(1,1), obj.BBcyclic(length(obj.BBcyclic),1), obj.considerNum);
            obj.targetBB = obj.simplifiedData(obj.BBcyclic, xRange);
        end
        
        function obj = preRender(obj)
            %load template;
            filePath = 'Tcl_Template/ModifiedCPTemplate.tcl';
            f = fopen(filePath, 'r');
            str = textscan(f, '%s','delimiter','\n');
            fclose(f);
            %str = strjoin(str{1,1}, '\n');
            str = str{1,1};
            str = strrep(str, '{{L}}', sprintf('%f', obj.L));
            str = strrep(str, '{{E}}', sprintf('%f', obj.E));
            str = strrep(str, '{{A}}', sprintf('%f', obj.A));
            str = strrep(str, '{{I}}', sprintf('%f', obj.I));
            str = strrep(str, '{{revE}}', sprintf('%f', obj.revKamp));
            str = strrep(str, '{{ampFactor}}', sprintf('%f', obj.ampFactor));
            str = strrep(str, '{{d_incr}}', sprintf('%f', obj.d_incr));
            str = strrep(str, '{{Naxial}}', sprintf('%f', obj.N));
            str = strrep(str, '{{DispList}}', obj.DispList);
            obj.template = str;
        end
        
        
        function obj = shiftBackbone(obj)
            obj.backboneShifted(:,1) = obj.backbone(:,1);
            obj.backboneShifted(:,2) = obj.backbone(:,1) .* obj.revK + obj.backbone(:,2);
            obj.Res = obj.backboneShifted(3,2) / obj.backboneShifted(2,2);
            x1 = obj.backboneShifted(2,1);
            y1 = obj.backboneShifted(2,2);
            x2 = obj.backboneShifted(3,1);
            y2 = obj.backboneShifted(3,2);
            x3 = x1 - y1 / (y1 - y2) * (x1 - x2);
            obj.backboneShifted(3,1) = x3;
            obj.backboneShifted(3,2) = 0.0;
        end
        
        function cmdLine = CP_cmdLine(obj, vector)
            lambda_S = vector(1);
            lambda_C = vector(2);
            lambda_A = vector(3);
            lambda_K = vector(4);
            cmdLine = sprintf('%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f', obj.K0, obj.as, obj.as, obj.Fy, -obj.Fy, lambda_S, lambda_C, lambda_A, lambda_K, obj.c_S, obj.c_C, obj.c_A, obj.c_K, obj.thetap, obj.thetap, obj.thetapc, obj.thetapc, obj.Res, obj.Res, obj.thetau, obj.thetau, obj.D, obj.D);
        end
        
        function renderTemplate(obj, outPath, vector)
            str = obj.template; 
            %replace parameters
            
            str = strrep(str, '{{CP_CMDLine}}', obj.CP_cmdLine(vector));
            
            str = strrep(str, '{{outname}}', outPath);
            %str = strrep(str, '{{}}', num2str());
            %str = strrep(str, '{{}}', num2str());
            f = fopen(sprintf('%s.out', outPath), 'w');
            fprintf(f, '%s\n', str{:});
            fclose(f);
            %dlmwrite(sprintf('%s.out', outPath),char(str),'delimiter','');
        end
        
        function output = runOpenSees(obj, filePath)
            [x, y] = system(sprintf('OpenSees %s.out', filePath));
            dataX = load(sprintf('%s_rotation.out', filePath));
            dataY = load(sprintf('%s_moment.out', filePath));
            delete *.out
            output = horzcat(dataX*-1, dataY);
        end
        
        function [output, fitness] = Analyze(obj, vector)
            
            outPath = sprintf('%s/rendered', obj.working_path);
            obj.renderTemplate(outPath, vector);
            output = runOpenSees(obj, outPath);
            %energy = obj.getEnergy(output(:,1), output(:,2));
            fitness = obj.Fitness(output);
        end
        
        function fitness = fit_fun(obj, vector)
            [output, fitness] = obj.Analyze(vector);
        end
        
        function fitness = Fitness(obj, data)
            data = obj.findTurning(data);
            data = obj.monoData(data(:,2:3));
            data = obj.simplifiedData(data, obj.targetBB(:,1));
            
            fitness = -1 * sum(abs(data(:,2) - obj.targetBB(:,2)));
        end
        
        function obj = convertDisp(obj)
            %find turning point
            newX = obj.turning(:,2);
            newX = vertcat(obj.targetX(1),newX);
            obj.DispList = sprintf('%s ', newX);

%             out = cell(1, length(obj.targetX));
%             for i = 1 : length(obj.targetX)
%                 newX(i) = obj.targetX(i);
%                 out{i} = sprintf('%f', newX(i));
%             end
%             obj.DispList = strjoin(out);
        end
        
        function energy = getEnergy(obj, dataX, dataY)
            data_d = dataX(2:length(dataX)) - dataX(1:length(dataX)-1);
            energy(:,1) = cumsum(abs(data_d));
            energy(:,2) = cumsum(data_d .* (dataY(1:length(dataY)-1) + dataY(2:length(dataY))) ./ 2.0);
        end
        
        function data_simplifed = simplifiedData(obj, data, xRange)
            data_simplifed(:,1) = xRange;
            data_simplifed(:,2) = interp1(data(:,1), data(:,2), xRange, 'linear','extrap');
        end
        
        function output = monoData(obj, data)
            output(1,:) = data(1,:);
            pre_x = output(1,1);
            for i=2:length(data)
                if data(i,1) > pre_x
                    output = vertcat(output, data(i,:));
                    pre_x = data(i,1);
                end
            end
        end
        
        function output = findTurning(obj, data)
            %find turning point
            len = length(data);
            x_temp = data(2:len,1) - data(1:len-1,1);
            pre_dx = x_temp(1);
            i = 2;
            output = [];
            while i<len-1
                if pre_dx ~= 0.0
                    cur_dx = x_temp(i);
                    if pre_dx * cur_dx < 0
                        output = vertcat(output, [i, data(i,1), data(i,2)]);
                        pre_dx = cur_dx;
                    end
                    i = i + 1;
                else
                    pre_dx = x_temp(i);
                end
            end
            output = vertcat(output, [len, data(len,1), data(len,1)]);
        end
    end
    
end

