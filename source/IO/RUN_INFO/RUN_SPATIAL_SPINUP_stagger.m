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
            unfinished_runs = distribute_remaining_runs(run_info);
            fprintf("Remaining runs: %d", sum(unfinished_runs > 0))
            
            if run_info.PARA.number_of_cores > 1
                poolobj = gcp('nocreate'); % If no pool, do not create new one.
                if isempty(poolobj)
                    parpool("Processes", run_info.PARA.number_of_cores)
                end

                spmd
                    [run_info, tile] = run_model_parallel(run_info, unfinished_runs);
                end
                
                delete(poolobj);
            else
                [run_info, tile] = run_model_sequential(run_info, unfinished_runs);
            end
        end

        function [run_info, tile] = run_model_parallel(run_info, unfinished_runs)
            tile = 0;
            worker_number = spmdIndex();
            
            pause_duration = (worker_number-1) * run_info.PARA.stagger_interval;
            fprintf('Worker %d paused for %d seconds to stagger start times...\n', worker_number, pause_duration);
            pause(pause_duration); %stagger worker start times to reduce file access conflicts

            runs = unfinished_runs(worker_number, :);
            runs = runs(runs > 0);  % some runs may be assigned a 0 index as padding

            if size(runs, 1) > 0
                for run_number = runs
                        [run_info, tile] = run_TILE(run_info, worker_number, run_number);
                end
            end
        end

        function [run_info, tile] = run_model_sequential(run_info, unfinished_runs)
            tile = 0;
            unfinished_runs = unfinished_runs(unfinished_runs > 0);
            for run_number = unfinished_runs
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
        function distributed_run_numbers = distribute_remaining_runs(run_info)
            unfinished_runs = get_unfinished_runs(run_info);
            ncpus = min(run_info.PARA.number_of_cores, size(unfinished_runs,1));

            padsize = ncpus- mod(size(unfinished_runs, 1), ncpus);
            runs_padded = padarray(unfinished_runs, padsize, 'post');
            distributed_run_numbers = reshape(runs_padded, ncpus, []);
        end

        function unfinished_runs = get_unfinished_runs(run_info)
            
            num_spatial_points = size(run_info.SPATIAL.STATVAR.key,1);
            run_finished = zeros(num_spatial_points, 1);
            for i = 1:num_spatial_points
                fname = make_final_tile_output_fname(run_info, i);
                run_finished(i, 1) = isfile(fname);
            end
            
            unfinished_runs = find(~run_finished);

        end
        
        function name = make_final_tile_output_fname(run_info, run_number)
            % This is quite an ugly function in that it doesn't generalise
            % to other configs, but it does the trick - creates the
            % filename of the final output
            
            run_name = run_info.PPROVIDER.PARA.run_name;
            classes = run_info.PPROVIDER.CLASSES;
            
            start_end_time_forcing_classes = classes.set_start_end_time;
            num_classes = length(start_end_time_forcing_classes);
            start_end_time = start_end_time_forcing_classes{num_classes};
            end_time_vector = start_end_time.PARA.end_time;
            end_time_str = join(string(end_time_vector), "");

            stem = join(string([run_name; run_number; end_time_str]), "_");
            name = join(string([stem; ".mat"]), "");
        end

        function [run_info, tile] = run_TILE(run_info, worker_number, run_number)
            % Shared spin-up sequence across TILE parallel/sequential.
            
            fname = make_final_tile_output_fname(run_info, run_number);
            
            if isfile(fname)
                fprintf("Final output found, skipping run [run_number=%d]\n", run_number)
                tile = 0;
                return 
            end

            disp(['running grid cell ' num2str(run_number)])
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



