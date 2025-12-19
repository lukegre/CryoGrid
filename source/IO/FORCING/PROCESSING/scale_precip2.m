%========================================================================
% CryoGrid FORCING processing class
%
%
% Authors:
% S. Westermann, December 2022
% L. Gregor, Nov 2025 - total_fraction parameter added to scale total precipitation - run before splitting rain/snow
%========================================================================

classdef scale_precip2< matlab.mixin.Copyable 
    
    properties
        CONST
        PARA
        STATVAR
        TEMP
    end
    
    methods
        function proc = provide_PARA(proc)
            proc.PARA.total_fraction = [];  %total precipitation is scaled - must be run before splitting rain/snow
            proc.PARA.rain_fraction = [];  %rainfall fraction assumed in sumulations (rainfall from the forcing data file is multiplied by this parameter)
            proc.PARA.snow_fraction = [];  %snowfall fraction assumed in sumulations (snowfall from the forcing data file is multiplied by this parameter)
        end
        
        
        function proc = provide_CONST(proc)
        end
        
        
        function proc = provide_STATVAR(proc)
        end
        
        
        function proc = finalize_init(proc, tile)
        end
        
        
        function forcing = process(proc, forcing, tile)
            if length(proc.PARA.total_fraction) == 1
                forcing.DATA.precip = forcing.DATA.precip .* proc.PARA.total_fraction;
            elseif or(length(proc.PARA.rain_fraction) == 1, length(proc.PARA.snow_fraction) == 1)
                forcing.DATA.rainfall = forcing.DATA.rainfall .* proc.PARA.rain_fraction;
                forcing.DATA.snowfall = forcing.DATA.snowfall .* proc.PARA.snow_fraction;
            end

        end
        
    end
    
end

