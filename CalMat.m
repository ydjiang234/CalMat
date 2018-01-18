classdef CalMat
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        E; K; A; I; L; N; backbone; targetX; targetY;%input parameters, to descripe the target member and target data
        K0; K0amp; as, thetap,thetapc; thetay; Fy; thetac; Fc; thetau; Res; backboneShifted;%describe the backbone of the CP model
        c_S = 1.0; c_C = 1.0; c_A = 1.0; c_K = 1.0;
        D = 1.0;
        revK; revKamp;%describe the additional reversed material
        ampFactor;
        DispList;
    end
    
    methods
        function obj = CalMat(A, I, L, N, revK, backbone, targetData, ampFactor)
            obj.A = A;
            obj.I = I;
            obj.L = L;
            obj.N = N;
            obj.revK =revK;
            obj.backbone = backbone;
            obj.targetX = targetData(:,1);
            obj.targetY = targetData(:,2);
            obj.ampFactor = ampFactor;
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
            obj.revKamp = obj.revK * (obj.ampFactor + 1.0);
            %other
            obj = obj.convertDisp();
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
        
        function cmdLine = CP_cmdLine(obj, lambda_S, lambda_C, lambda_A, lambda_K)
            cmdLine = sprintf('%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f', obj.K0amp, obj.as, obj.as, obj.Fy, -obj.Fy, lambda_S, lambda_C, lambda_A, lambda_K, obj.c_S, obj.c_C, obj.c_A, obj.c_K, obj.thetap, obj.thetap, obj.thetapc, obj.thetapc, obj.Res, obj.Res, obj.thetau, obj.thetau, obj.D, obj.D);
        end
        
        function str = renderTemplate(obj, filePath, outPath, lambda_S, lambda_C, lambda_A, lambda_K )
            f = fopen(filePath, 'r');
            str = textscan(f, '%s','delimiter','\n');
            fclose(f);
            %str = strjoin(str{1,1}, '\n');
            str = str{1,1};
            
            %replace parameters
            str = strrep(str, '{{L}}', num2str(obj.L));
            str = strrep(str, '{{E}}', num2str(obj.E));
            str = strrep(str, '{{A}}', num2str(obj.A));
            str = strrep(str, '{{I}}', num2str(obj.I));
            str = strrep(str, '{{revE}}', num2str(obj.revKamp));
            str = strrep(str, '{{ampFactor}}', num2str(obj.ampFactor));
            str = strrep(str, '{{Naxial}}', num2str(obj.N));
            str = strrep(str, '{{CP_CMDLine}}', obj.CP_cmdLine(lambda_S, lambda_C, lambda_A, lambda_K));
            str = strrep(str, '{{DispList}}', obj.DispList);
            %str = strrep(str, '{{}}', num2str());
            %str = strrep(str, '{{}}', num2str());
            
            %dlmwrite(outPath,str,'delimiter','');
        end
        
        function runOpenSees(obj, filePath)
            system(sprintf('OpenSees %s', filePath));
        end
        
        function obj = convertDisp(obj)
            out = cell(1, length(obj.targetX));
            for i=1:length(out)
                out{i} = num2str(obj.targetX(i));
            end
            obj.DispList = strjoin(out);
        end
    end
    
end

