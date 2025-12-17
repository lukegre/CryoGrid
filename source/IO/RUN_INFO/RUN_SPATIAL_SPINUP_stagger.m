%========================================================================
% CryoGrid RUN_INFO class RUN_SPATIAL_SPINUP
% RUN_INFO class for spatially distributed runs (using an appropriate 
% SPATIAL_REFERENCE class, DATA_PROVIDER classes and FORCING classes)
% which can run several TILE classes per point sequentially for model spin-up 
%
% S. westermann, Dec 2022
% L. Gregor, Dec 2025 - Removed MULTITILE functionality, switched to parfor, added worker stagger start
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

    properties (Constant)
        stagger_interval = 10; %seconds to stagger worker start times
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
            tile = 0;
            if run_info.PARA.number_of_cores > 1
                parpool(run_info.PARA.number_of_cores)
                spmd
                    [run_info, tile] = run_model_TILE_parallel_parfor(run_info);
                end
            else
                [run_info, tile] = run_model_TILE_sequential(run_info);
            end
        end

        function [run_info, tile] = run_model_TILE_parallel_parfor(run_info)
            tile = 0;
            number_spatial_points = size(run_info.SPATIAL.STATVAR.key, 1);
            stagger = run_info.stagger_interval;

            parfor run_number = 1:number_spatial_points
                % Get the unique ID of the worker currently running this iteration
                t = getCurrentTask(); 
                worker_id = t.ID;

                % Use persistent var that lives on the worker's process
                persistent worker_initialized;  
                if isempty(worker_initialized)  % block ONLY runs 1st time a worker executes code 
                    pause_time = (worker_id - 1) * stagger;
                    
                    fprintf('Worker %d: First time initialization. Staggering %d sec.\n', worker_id, pause_time);
                    pause(pause_time);
                    
                    worker_initialized = true;
                end
                run_spinup_for_range(run_info, worker_id, run_number, false);
            end
        end

        function [run_info, tile] = run_model_TILE_parallel_spmd(run_info)
            tile = 0;
            worker_number = spmdIndex();

            number_spatial_points = size(run_info.SPATIAL.STATVAR.key,1);
            max_number_of_gridcells = min(run_info.PARA.number_of_cores, size(run_info.SPATIAL.STATVAR.key,1));
            number_of_runs = number_spatial_points ./ max_number_of_gridcells;
            run_raster = [0; round([number_of_runs : number_of_runs : number_spatial_points]')];
            run_raster = [run_raster(1:end-1,1)+1 run_raster(2:end,1)];

            if worker_number <= size(run_raster,1)
                for run_number = run_raster(worker_number,1):run_raster(worker_number,2)
                    disp(['running grid cell ' num2str(run_number)])
                    [run_info, tile] = run_spinup_for_range(run_info, worker_number, run_number, false);
                end
            end
        end

        function [run_info, tile] = run_model_TILE_sequential(run_info)
            tile = 0;
            for run_number = 1:size(run_info.SPATIAL.STATVAR.key,1)

                disp(['running grid cell ' num2str(run_number)])
                [run_info, tile] = run_spinup_for_range(run_info, 1, run_number, true);
            end
        end

        %-------------param file generation-----
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
        
        function [run_info, tile] = run_spinup_for_range(run_info, worker_number, range, display_latlon)
            % Shared spin-up sequence across TILE parallel/sequential.
            tile = 0;
            for i=1:size(run_info.PARA.tile_class,1)
                disp(['running tile number ' num2str(i)])
                for j=1:run_info.PARA.number_of_runs_per_tile(i,1)
                    disp(['running round ' num2str(j)])

                    for ai=1:size(run_info.SPATIAL.ACTION,1)
                        run_info.SPATIAL.ACTION{ai,1} = assign_tile_properties(run_info.SPATIAL.ACTION{ai,1}, range); %writes the provider class
                    end

                    new_tile = copy(run_info.PPROVIDER.CLASSES.(run_info.PARA.tile_class{i,1}){run_info.PARA.tile_class_index(i,1),1});
                    new_tile.RUN_INFO = run_info;
                    new_tile = finalize_init(new_tile);
                    tile = new_tile;
                    run_info.TILE = tile;

                    tile.PARA.worker_number = worker_number;
                    tile.PARA.range = range;

                    if display_latlon
                        fprintf('LAT: %.3f, LON: %.3f\n', tile.PARA.latitude, tile.PARA.longitude);
                    end

                    tile = run_model(tile);  %time integration
                end
            end
        end  % end of function
    end
end



