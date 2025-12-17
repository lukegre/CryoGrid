%========================================================================
% CryoGrid RUN_INFO class RUN_SPATIAL_SPINUP
% RUN_INFO class for spatially distributed runs (using an appropriate 
% SPATIAL_REFERENCE class, DATA_PROVIDER classes and FORCING classes)
% which can run several TILE classes per point sequentially for model spin-up 
%
% S. westermann, Dec 2022
% L. Gregor, Dec 2025 - Removed MULTITILE functionality, refactored, added worker stagger start
%========================================================================

classdef RUN_SPATIAL_SPINUP_stagger < matlab.mixin.Copyable
    
    properties
        PPROVIDER
        PARA
        CONST
        STATVAR
        SPATIAL
        TILE
        CLUSTER
    end

    methods
        
        function run_info = provide_PARA(run_info)
            
            run_info.PARA.parallel_pool_open = 0;

            run_info.PARA.number_of_cores = [];
            
            run_info.PARA.tile_class = [];
            run_info.PARA.tile_class_index = [];
            run_info.PARA.number_of_runs_per_tile = []; %vector
            
            run_info.PARA.projection_class = [];
            run_info.PARA.projection_class_index = [];
            run_info.PARA.stagger_interval = 0; %seconds to stagger worker start times
                        
        end
        
        function run_info = provide_CONST(run_info)
        end
        
        function run_info = provide_STATVAR(run_info)
        end
        
        function run_info = finalize_init(run_info)
        
            if ~isempty(run_info.PARA.projection_class) && ~(sum(isnan(run_info.PARA.projection_class))>0)
                disp('get spatial data')
                spatial_class = copy(run_info.PPROVIDER.CLASSES.(run_info.PARA.projection_class){run_info.PARA.projection_class_index,1});
                if ~spatial_class.PARA.new_reference
                    spatial_class.STATVAR = run_info.SPATIAL.STATVAR;
                end
                run_info.SPATIAL = spatial_class;
                run_info.SPATIAL.RUN_INFO = run_info;
                run_info.SPATIAL = finalize_init(run_info.SPATIAL);
            end
        end

        function [run_info, tile] = run_model(run_info)
            if run_info.PARA.number_of_cores > 1
                poolobj = gcp('nocreate'); % If no pool, do not create new one.
                if isempty(poolobj)
                    parpool(run_info.PARA.number_of_cores)
                end

                spmd
                    [run_info, tile] = run_model_parallel(run_info);
                end
                
                delete(poolobj);
            else
                [run_info, tile] = run_model_sequential(run_info);
            end
        end

        function [run_info, tile] = run_model_parallel(run_info)
            tile = 0;
            worker_number = spmdIndex();

            pause_duration = (worker_number-1) * run_info.PARA.stagger_interval;
            fprintf('Worker %d pausing for %d seconds to stagger start times...\n', worker_number, pause_duration);
            pause(pause_duration); %stagger worker start times to reduce file access conflicts

            number_spatial_points = size(run_info.SPATIAL.STATVAR.key,1);
            max_number_of_gridcells = min(run_info.PARA.number_of_cores, size(run_info.SPATIAL.STATVAR.key,1));
            number_of_runs = number_spatial_points ./ max_number_of_gridcells;
            run_raster = [0; round([number_of_runs : number_of_runs : number_spatial_points]')];
            run_raster = [run_raster(1:end-1,1)+1 run_raster(2:end,1)];

            if worker_number <= size(run_raster,1)
                for run_number = run_raster(worker_number,1):run_raster(worker_number,2)
                    disp(['running grid cell ' num2str(run_number)])
                    [run_info, tile] = run_TILE(run_info, worker_number, run_number);
                end
            end
        end

        function [run_info, tile] = run_model_sequential(run_info)
            tile = 0;
            for run_number = 1:size(run_info.SPATIAL.STATVAR.key,1)

                disp(['running grid cell ' num2str(run_number)])
                [run_info, tile] = run_TILE(run_info, 1, run_number);
            end
        end

        function run_info = param_file_info(run_info)
            run_info = provide_PARA(run_info);

            run_info.PARA.STATVAR = [];
            run_info.PARA.class_category = 'RUN_INFO';
            run_info.PARA.default_value = [];
            run_info.PARA.comment = [];
            
            run_info.PARA.comment.number_of_cores = {'number of cores to be used for calculation'};
            run_info.PARA.default_value.number_of_cores = {2};
            
            run_info.PARA.options.tile_class.name =  'H_LIST';
            run_info.PARA.options.tile_class.entries_x = {'TILE_1D_standard' 'TILE_1D_standard'};
            
            run_info.PARA.options.tile_class_index.name =  'H_LIST'; 
            run_info.PARA.options.tile_class_index.entries_x = {1 2};
            
            run_info.PARA.options.number_of_runs_per_tile.name =  'H_LIST'; % 
            run_info.PARA.options.number_of_runs_per_tile.entries_x = {1 1};
            
            run_info.PARA.comment.projection_class = {'projection class providing providing information on the locations and additinal data for each target point'};
            
        end
    end

    methods (Access = private)
        function [run_info, tile] = run_TILE(run_info, worker_number, run_number)
            % Shared spin-up sequence across TILE parallel/sequential.
            tile = 0;
            for i=1:size(run_info.PARA.tile_class,1)
                disp(['running tile number ' num2str(i)])
                for j=1:run_info.PARA.number_of_runs_per_tile(i,1)
                    if j > 1
                        disp(['running round ' num2str(j)])
                    end

                    for ai=1:size(run_info.SPATIAL.ACTION,1)
                        run_info.SPATIAL.ACTION{ai,1} = assign_tile_properties(run_info.SPATIAL.ACTION{ai,1}, run_number); %writes the provider class
                    end

                    new_tile = copy(run_info.PPROVIDER.CLASSES.(run_info.PARA.tile_class{i,1}){run_info.PARA.tile_class_index(i,1),1});
                    new_tile.RUN_INFO = run_info;
                    new_tile = finalize_init(new_tile);
                    tile = new_tile;
                    run_info.TILE = tile;

                    tile.PARA.worker_number = worker_number;
                    tile.PARA.range = run_number;

                    tile = run_model(tile);  %time integration
                
                end
            end
        end  % end of function
    end

end



