%========================================================================
% CryoGrid SPATIAL_REFERENCE class LIST_OF_POINTS
% S. Westermann, Oct 2025
% L. Gregor, Dec 2025 - import points from CSV (parameter_list_csv_fname) - must match excel point columns
%========================================================================


classdef LIST_OF_POINTS_csv < matlab.mixin.Copyable

    properties
        RUN_INFO
        PARA
        CONST
        STATVAR
        TEMP
        ACTION
    end
    
    methods
        function proj = provide_PARA(proj)
            %optoional parameters, can be used to set constant values for
            %all points, overwritten if they are contained in
            %parameter_list
            proj.PARA.latitude = [];
            proj.PARA.longitude = [];
            proj.PARA.altitude = [];
            proj.PARA.slope_angle = [];     %
            proj.PARA.aspect = [];     %
            proj.PARA.area = [];
            
            proj.PARA.parameter_list_csv_fname = [];
            proj.PARA.parameter_list = []; %first row contains variable names

            proj.PARA.mask_class = []; %acts on the entire 2d matirx
            proj.PARA.mask_class_index = [];
            
            proj.PARA.data_class = [];
            proj.PARA.data_class_index = [];
            
            proj.PARA.data_mask_class = []; %
            proj.PARA.data_mask_class_index = [];
            
            proj.PARA.assign_tile_properties_class = [];
            proj.PARA.assign_tile_properties_class_index = [];

            proj.PARA.new_reference = 1;
        end
        
        function proj = provide_STATVAR(proj)

        end
        
        function proj = provide_CONST(proj)
            
        end
        
        function proj = finalize_init(proj)
            
            proj = get_parameter_list(proj);
            
            vars = fieldnames(proj.PARA.parameter_list);

            for i=1:size(vars,1)
                if ~strcmp(vars{i,1}, 'depth')
                    proj.STATVAR.(vars{i,1}) = proj.PARA.parameter_list.(vars{i,1});
                    number_of_points = size(proj.STATVAR.(vars{i,1}),1);
                end
            end
            vars = {'latitude'; 'longitude'; 'altitude'; 'slope_angle'; 'aspect'; 'area'};
            for i=1:size(vars,1)
                if ~any(strcmp(vars{i,1}, fieldnames(proj.STATVAR)))
                    proj.STATVAR.(vars{i,1}) = repmat(proj.PARA.(vars{i,1}), number_of_points,1);
                end
            end
           
            proj.STATVAR.key = [1:size(proj.STATVAR.latitude,1)]';
            
            %apply masks before data sets
            proj.STATVAR.mask = logical(proj.STATVAR.longitude.*1);
            for i=1:size(proj.PARA.mask_class_index,1)
                mask_class = copy(proj.RUN_INFO.PPROVIDER.CLASSES.(proj.PARA.mask_class{i,1}){proj.PARA.mask_class_index(i,1),1});
                mask_class.PARENT = proj;
                mask_class = finalize_init(mask_class);
                mask_class = apply_mask(mask_class); %can be additive or subtractive
            end

            %reduce the list to the ones inside the masks
            mask = proj.STATVAR.mask;
            fn = fieldnames(proj.STATVAR);
            for i=1:size(fn,1)  
                    proj.STATVAR.(fn{i,1})(~mask) = [];
            end


            %load data sets
            for i=1:size(proj.PARA.data_class,1)
                data_class = copy(proj.RUN_INFO.PPROVIDER.CLASSES.(proj.PARA.data_class{i,1}){proj.PARA.data_class_index(i,1),1});
                data_class.PARENT = proj;
                data_class = finalize_init(data_class);
                data_class = load_data(data_class); %can be additive or subtractive
            end

            for i=1:size(proj.PARA.data_mask_class_index,1)
                mask_class = copy(proj.RUN_INFO.PPROVIDER.CLASSES.(proj.PARA.data_mask_class{i,1}){proj.PARA.data_mask_class_index(i,1),1});
                mask_class.PARENT = proj;
                mask_class = finalize_init(mask_class);
                mask_class = apply_mask(mask_class); %can be additive or subtractive
            end

            mask = proj.STATVAR.mask;
            fn = fieldnames(proj.STATVAR);
            for i=1:size(fn,1)  
                    proj.STATVAR.(fn{i,1})(~mask) = [];
            end
            
            for i=1:size(proj.PARA.assign_tile_properties_class,1)
                proj.ACTION{i,1} = copy(proj.RUN_INFO.PPROVIDER.CLASSES.(proj.PARA.assign_tile_properties_class{i,1}){proj.PARA.assign_tile_properties_class_index(i,1),1});
                proj.ACTION{i,1} = finalize_init(proj.ACTION{i,1});
                proj.ACTION{i,1}.PROJ = proj;
            end

        end

    end

    methods (Access = private)
        function proj = get_parameter_list(proj)
            % if PARA does not have parameter_list but has parameter_csv_fname
            file_provided = ~isempty(proj.PARA.parameter_list_csv_fname);
            list_provided = ~isempty(proj.PARA.parameter_list);

            if and(file_provided, list_provided)
                error("You have provided both parameter_list_csv_fname and parameter_list to LIST_OF_POINTS2")
            elseif ~or(file_provided, list_provided)
                error("You must provide at least parameter_list_csv_fname or parameter_list to LIST_OF_POINTS2")
            elseif list_provided
                disp("Using list given in the excel config file for LIST_OF_POINTS2")
                return 
            elseif file_provided
                disp(['Getting parameter_list from the given filename: ' proj.PARA.parameter_list_csv_fname])
                table = readtable(proj.PARA.parameter_list_csv_fname);
                data = table2struct(table, 'ToScalar',true);
                proj.PARA.parameter_list = data;
            
            end
        end
    end
end

