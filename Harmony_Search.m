classdef Harmony_Search
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        ins_num; pit_range; hms; fw_ratio; fw; hmcr; par;
        fit_fun;
        HM;
        fitness;
    end
    
    methods
        function obj = Harmony_Search(pit_range, hms, fit_fun, hmcr, par, fw_ratio)
            obj.pit_range =pit_range;
            obj.hms = hms;
            obj.fit_fun = fit_fun;
            obj.hmcr = hmcr;
            obj.par = par;
            obj.fw_ratio = fw_ratio;
            obj.ins_num = length(pit_range);
            obj = obj.covert_default();
            obj = obj.generate_HM();
        end
        
        function obj = covert_default(obj)
            if length(obj.hmcr) == 1
                obj.hmcr = ones(1, obj.ins_num) * obj.hmcr;
            end
            if length(obj.par) == 1
                obj.par = ones(1, obj.ins_num) * obj.par;
            end
            obj.fw = obj.fw_ratio .* (obj.pit_range(:,2) - obj.pit_range(:,1));
        end
        
        function bool = bool_probability(obj, probability)
            temp  = rand;
            if temp <= probability
                bool = 1;
            else
                bool = 0;
            end
        end
        
        function out = random_pit_from_list(obj, target_list)
            ind = floor(length(target_list) * rand) + 1;
            out =  target_list(ind);
        end
        
        function new_pit = select_pitch(obj, ins_ind)
            if obj.bool_probability(obj.hmcr(ins_ind))
                new_pit = obj.random_pit_from_list(obj.HM(:,ins_ind));
                if obj.bool_probability(obj.par(ins_ind))
                    new_pit = new_pit + obj.fw(ins_ind) * obj.randNum(-1.0, 1.0);
                end
            else
                new_pit = obj.randNum(obj.pit_range(ins_ind,1), obj.pit_range(ins_ind,2));
            end
        end
        
        function out = randNum(obj, a, b)
            out = a + (b-a)*rand;
        end
        
        function new_vector = new_harmony_vector(obj)
            new_vector = zeros(1, obj.ins_num);
            for ins_ind =1 : obj.ins_num
                new_vector(ins_ind) = obj.select_pitch(ins_ind);
            end
        end
        
        function obj = generate_HM(obj)
            obj.HM = zeros(obj.hms, obj.ins_num);
            obj.fitness = zeros(1, obj.hms);
            for i = 1:obj.hms
                for ins_ind = 1: obj.ins_num
                    obj.HM(i, ins_ind) = obj.randNum(obj.pit_range(ins_ind,1), obj.pit_range(ins_ind,2));
                end
                obj.fitness(i) = obj.fit_fun(obj.HM(i,:));
            end
        end
        
        function obj = update_HM(obj, harmony_vector)
            new_vector = harmony_vector;
            new_f = obj.fit_fun(new_vector);
            [fitness_min, ind_min] = min(obj.fitness);
            if new_f > fitness_min
                obj.fitness(ind_min) = new_f;
                obj.HM(ind_min, :) = new_vector;
            end
        end
        
        function obj = next(obj)
            obj = obj.update_HM(obj.new_harmony_vector());
        end
        
        function [vector, fitness_max] = Optimized(obj)
            [fitness_max, ind_max] = max(obj.fitness);
            vector = obj.HM(ind_max,:);
        end
        
    end
    
end

