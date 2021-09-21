classdef Analysis2P
    %ANALYSIS2P analyze 2p data by summarizing individual results according
    %to Session2P or other function
    %   Detailed explanation goes here
    
    properties
        name
        path
    end
    
    methods
        function obj = Analysis2P(inputArg1,inputArg2)
            %ANALYSIS2P Construct an instance of this class
            %   Detailed explanation goes here
            obj.Property1 = inputArg1 + inputArg2;
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

