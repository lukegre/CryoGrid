%========================================================================
% CryoGrid FORCING processing class
%%uses "normal" interpolation of Sin, not the superior algroithm relying on
%S_TOA used in process_topoScale2
%
% Authors:
% S. Westermann, December 2022
%
%========================================================================

classdef process_topoScale3 < process_BASE
    
    properties (Constant)
        %constants
        gas_constant_for_dry_air = 287.05;  % Gas constant for dry air [JK^-1kg^-1]
        universal_gas_constant = 8.3144598
        molar_mass_of_air = 0.0289644;
        epsilon0 = 0.622; % Ratio of molecular weight of water and dry air [-]
        stefan_bolzman_constant = 5.67e-8;
        
        gravity = 9.81; % Acceleration of gravity [ms^-1]
        
        solar_constant = 1370; % Solar constat (total TOA solar irradiance) [Wm^-2] used in ECMWF's IFS

        precipitation_adjustment_factor = 0.27;  % Mean adjustment factor (following Fiddes and Gruber, 2014)
    end

    methods
        function proc = provide_PARA(proc)
        end
        
        function proc = provide_CONST(proc)
            proc.CONST.Tmfw = [];
            proc.CONST.sigma = [];
            proc.CONST.limit_above_orography = 1;
        end
        
        function proc = provide_STATVAR(proc)
        end
        
        function proc = finalize_init(proc, tile)
        end
        
        function forcing = process(proc, forcing, tile)
            
            disp('applying downscaling with TopoScale')
            era = forcing.TEMP.era;

            lat_point = forcing.SPATIAL.STATVAR.latitude;
            lon_point = forcing.SPATIAL.STATVAR.longitude;
            alt_point = forcing.SPATIAL.STATVAR.altitude;
            point = struct('lat', lat_point, 'lon', lon_point, 'alt', alt_point);

            lat_grid = era.lat;
            lon_grid = era.lon;
            
            n_timesteps = size(era.t, 2);
            n_levels = size(era.p, 2);
            n_coords = 4;
            

            % PREPARING COORDINATE WEIGHTS FOR WEIGHTED AVERAGES
            [weights_lat, weights_lon, ind_lat, ind_lon] = proc.calc_horz_weights(point.lat, point.lon, lat_grid, lon_grid, n_timesteps);
            inds = struct('lat', ind_lat, 'lon', ind_lon);
            
            era_alt    = reshape(era.Z (inds.lon, inds.lat, :, :), n_coords, n_levels, n_timesteps);
            era_alt_sl = reshape(era.Zs(inds.lon, inds.lat), n_coords, 1);
            era_alt_sl = repmat(era_alt_sl, 1, 1, n_timesteps);  % lat_lon, level, timestep
            
            limit_above_orography = proc.CONST.limit_above_orography;
            [weights_Z_below, weights_Z_above, factor] = proc.calc_vert_weights(point.alt, era_alt, era_alt_sl, n_levels, limit_above_orography);
            weights = struct('lat', weights_lat, 'lon', weights_lon, 'Z_above', weights_Z_above, 'Z_below', weights_Z_below, 'factor', factor);
            

            % FUNCTIONS FOR SELECTING DATA
            select_surf = @(var) double(var(inds.lon, inds.lat, :));
            % select_surf -> (2 x 2 x t)
            select_3D = @(var) double(var(inds.lon, inds.lat, :, :));
            % select_3d -> (2 x 2 x levels x t)
            
            
            % DOWNSCALING TEMPERATURE 
            era_T = select_3D(era.T) .* era.T_sf;
            era_T_sl = select_surf(era.T2) .* era.T_sf;
            % downscale_3D first does vertical downscaling, then horizontal using downscale_2D
            % downscale_3D -> (t x 1)
            T_topoScale = proc.downscale_3D_var(era_T, era_T_sl, weights, point.alt, era_alt_sl);
            
            
            % DOWNSCALING WIND SPEED
            era_u = select_3D(era.u) .* era.wind_sf;
            era_v = select_3D(era.v) .* era.wind_sf;
            era_u10 = select_surf(era.u10) .* era.wind_sf;
            era_v10 = select_surf(era.v10) .* era.wind_sf;
            era_wind = proc.calc_wind_speed(era_u, era_v);
            era_wind_sl = proc.calc_wind_speed(era_u10, era_v10);
            wind_topoScale = proc.downscale_3D_var(era_wind, era_wind_sl, weights, point.alt, era_alt_sl);
            

            % DOWNSCALING PRESSURE
            % Since pressure does not exist in ERA5, we compute it from Z-pressure levels
            era_p = repmat(era.p, n_coords, 1, n_timesteps);  % dims = 4, 6, t
            era_p_sl  = select_surf(era.ps)  .* era.ps_sf;  % surface pressure
            p_topoScale = sum(era_p .* (weights.Z_above + weights.Z_below), 2);
            p_topoScale = proc.barometric_formula(p_topoScale, era_p_sl, point.alt, era_alt_sl, era_T_sl);
            p_topoScale = proc.downscale_2D_var(p_topoScale, weights);


            % DOWNSCALING HUMIDITY
            % era_q = select_3D(era.q) .* era.q_sf;  % humidity
            era_Td_sl = select_surf(era.Td2) .* era.T_sf;  % dew point temp
            % era_q_sl = proc.calc_humidity(era_Td_sl, era_p_sl);  % surface humidity
            era_RH_sl = proc.magnus_formula(era_Td_sl) ./ proc.magnus_formula(era_T_sl);
            RH_topoScale_sl = proc.downscale_2D_var(era_RH_sl, weights);
            % in the original code, era_q was downscaled, and later q_topoScale recomputed with the code below. 
            % The latter is more accurate, so I use this directly
            q_topoScale = proc.epsilon0 .* RH_topoScale_sl .* proc.magnus_formula(T_topoScale) ./ p_topoScale;
            % q_topoScale = proc.downscale_3D_var(era_q, era_q_sl, weights, point.alt, era_alt_sl);
            

            % DOWNSCALING LONGWAVE RADIATION
            % in the original process_topoScale function, humidity used for downscaling is not the final version of the humidity that is
            % later stored - in this version, we use the final humidity, meaning that results are slightly different to the original. 
            sbc = proc.stefan_bolzman_constant;
            era_vapour_press_sl = proc.magnus_formula(era_Td_sl);

            era_Lin_sl = select_surf(era.LW) .* era.rad_sf;
            era_Lin_sl = proc.downscale_2D_var(era_Lin_sl, weights);
            era_TK_sl = proc.downscale_2D_var(era_T_sl, weights) + 273.15;
            era_vp_sl= proc.downscale_2D_var(era_vapour_press_sl, weights);
            TK_topoScale = T_topoScale + 273.15;
            
            mask = era_Lin_sl <= 0;
            era_Lin_sl(mask) = sbc .* era_TK_sl(mask).^4;
            vp_topoScale = proc.calc_vapour_pressure(q_topoScale, p_topoScale);
            
            % Use the vapor pressure and temperature to calculate clear sky emssivity at grid and subgrid.
            cef = proc.calc_clear_sky_emissivity(vp_topoScale, T_topoScale);
            cec = proc.calc_clear_sky_emissivity(era_vp_sl, era_TK_sl - 273.15);  % expects C
            aec = era_Lin_sl ./ (sbc .* era_TK_sl.^4);  % Diagnose the all sky emissivity at grid.
            deltae = aec - cec;  % Calculate the "cloud" emissivity at grid, assume this is the same at subgrid.
            aef = cef + deltae;  % Use the former cloud emissivity to compute the all sky emissivity at subgrid.
            Lin_topoScale = aef .* sbc .* TK_topoScale.^4;
            

            % DOWNSCALING SHORTWAVE RADIATION
            era_Sin_sl = select_surf(era.SW) .* era.rad_sf;
            era_Sin_sl = proc.downscale_2D_var(era_Sin_sl, weights);
            era_S_TOA_sl = select_surf(era.S_TOA) .* era.rad_sf;
            era_S_TOA_sl = proc.downscale_2D_var(era_S_TOA_sl, weights);
            era_p_sl2 = proc.downscale_2D_var(era_p_sl, weights);
            % solar zenith angle
            [~, solar_zen] = solargeom(proc, era.t, point.lat, point.lon);
            % determines the intensity of the sun rays (1 = overhead, 0 = horizon)
            [mu0, is_sunset] = proc.calc_solar_geometry(solar_zen');
            era_S_TOA_point = proc.solar_constant .* mu0; % The theoretical maximum radiation if there were no atmosphere
            kt = proc.calc_clearness_index(era_Sin_sl, era_S_TOA_sl);
            kd = proc.calc_diffuse_fraction(kt);

            % May also want to consider treating values mu0<0 for prominent topography when the horizon angles are less than 0.
            Sin_diff_topoScale = kd .* era_Sin_sl .* double(~is_sunset); % Diffuse shortwave radiation.
            Sin_dir_topoScale = era_Sin_sl - Sin_diff_topoScale; % Direct shortwave radiation
            Sin_dir_topoScale = era_S_TOA_point .* (Sin_dir_topoScale ./ (max(1e-10, era_S_TOA_point))).^(p_topoScale ./ era_p_sl2);


            % DOWNSCALING PRECIPITATION
            adjf = proc.precipitation_adjustment_factor; 
            era_precip_sl = select_surf(era.P) .* era.P_sf;
            era_precip_sl = proc.downscale_2D_var(era_precip_sl, weights);
            era_alt_sl2 = proc.downscale_2D_var(era_alt_sl, weights);

            %  Apply Liston & Elder (MicroMet) elevation-based precip adjustment
            dZ = point.alt - era_alt_sl2; % m
            dZ = dZ ./ 1e3; % km
            dZ = min(dZ, 2); % No larger that 2 km=3.3 adjustment
            dZ = max(dZ, -2);% For symmetry
            adj = (1 + adjf .* dZ) ./ (1 - adjf .* dZ);
            precip_topoScale = era_precip_sl .* adj .* 24; %in mm/day, check if timestep must be taken into account when not using 1h input data
            

            % HOUSEKEEPING AND SETTING UP STRUCTS
            forcing.DATA.Tair = double(T_topoScale);
            forcing.DATA.q = double(q_topoScale);
            forcing.DATA.wind = double(wind_topoScale); 
            forcing.DATA.Sin_dir = double(Sin_dir_topoScale);
            forcing.DATA.Sin_dif = double(Sin_diff_topoScale);
            forcing.DATA.Sin =  forcing.DATA.Sin_dir + forcing.DATA.Sin_dif;
            forcing.DATA.S_TOA = era_S_TOA_point;
            forcing.DATA.Lin = double(Lin_topoScale);
            forcing.DATA.p = double(p_topoScale);
            forcing.DATA.precip = double(precip_topoScale);
            forcing.DATA.timeForcing = era.t';

            forcing.TEMP.era = [];

        end
        
        function p = satPresIce(proc, T)
            Tmfw = proc.CONST.Tmfw;
            p = 6.112.* 100.* exp(22.46.*(T-Tmfw)./(272.61+T-Tmfw));
        end
        
        function p = satPresWater(proc, T)
            Tmfw = proc.CONST.Tmfw;
            p = 6.112 .* 100 .* exp(17.62.*(T-Tmfw)./(243.12+T-Tmfw));
        end
        
    end

    methods(Static)

        function var = update_variable(var, var_sl, merge_w_sl, use_sl, factor)
            var = double(~merge_w_sl) .* var + double(merge_w_sl) .* factor .* var + double(merge_w_sl) .* (1-factor) .* double(var_sl);
            var = double(~use_sl) .* var + double(use_sl) .* double(var_sl);
        end

        function [weights_lat, weights_lon, ind_lat, ind_lon] = calc_horz_weights(lat_point, lon_point, lat_grid, lon_grid, n_timesteps)

            dist_lat = abs(lat_point - lat_grid);
            dist_lon = abs(lon_point - lon_grid);

            [dist_lat, ind_lat] = sort(dist_lat);
            [dist_lon, ind_lon] = sort(dist_lon);
            
            dist_lat=dist_lat(1:2);
            dist_lon=dist_lon(1:2);
            
            ind_lat = ind_lat(1:2);
            ind_lon = ind_lon(1:2);

            weights_lat = 1 - dist_lat./sum(dist_lat);
            weights_lon = 1 - dist_lon./sum(dist_lon);

            weights_lat = repmat(weights_lat', 2, 1, 1, n_timesteps);
            weights_lat = reshape(weights_lat, 4, 1 , n_timesteps);
            
            weights_lon = repmat(weights_lon, 1, 2, 1, n_timesteps);
            weights_lon = reshape(weights_lon, 4, 1 , n_timesteps);
            
        end
        
        function [weights_above, weights_below, factor] = calc_vert_weights(alt_point, era_alt, era_alt_sl, n_levels, limit_above_orography)
            
            if nargin <= 4
                limit_above_orography = 1;
            elseif nargin == 5
                assert(limit_above_orography == 0 | limit_above_orography == 1, "limit_above_orography must be [0,1]");
            end
            layer_below = int16(era_alt .* 0);
            layer_above = int16(era_alt .* 0);
            
            %do another one to get the lowermost pl above the orography
            for i = 2:n_levels
                level_a = era_alt(:,i,:);
                level_b = era_alt(:,i-1,:);

                point_lt_grid1 = level_a < alt_point;
                grid0_gt_point = level_b >= alt_point;
                
                % this is added so that valleys can also be included for very mountainous areas
                if limit_above_orography
                    a_above_orography = level_a > era_alt_sl;
                    b_above_orography = level_b > era_alt_sl;
                else
                    a_above_orography = 1;
                    b_above_orography = 1;
                end

                layer_below(:,  i,:) = double(point_lt_grid1 & grid0_gt_point & a_above_orography);
                layer_above(:,i-1,:) = double(point_lt_grid1 & grid0_gt_point & b_above_orography);
            end
            
            distance_Z_above = abs(sum(era_alt .* layer_above, 2) - alt_point) .* double(sum(layer_above, 2) > 0);
            distance_Z_below = abs(sum(era_alt .* layer_below, 2) - alt_point) .* double(sum(layer_below, 2) > 0);
            
            weights_Z_above = 1-distance_Z_above ./ max(1e-10, distance_Z_above + distance_Z_below);
            weights_Z_below = 1-distance_Z_below ./ max(1e-10, distance_Z_above + distance_Z_below);
            
            weights_above = repmat(weights_Z_above, 1, size(layer_above,2),1) .* double(layer_above);
            weights_below = repmat(weights_Z_below, 1, size(layer_below,2),1) .* double(layer_below);

            factor = min(1, max(0, distance_Z_above./100));
            
        end
        
        function var = downscale_3D_var(var4d, var_surf, weights, alt_point, era_alt_sl)
            class = process_topoScale3;

            % Function can be applied to temperature (T), wind, and humidity (q)
            % takes var3d where variables have already been selected
            
            % VARIABLE PREPARATION
            [n_lon, n_lat, n_levels, n_timesteps] = size(var4d);
            var3d = reshape(var4d, n_lon * n_lat, n_levels, n_timesteps);
            var_surf = double(reshape(var_surf, n_lon * n_lat, 1, n_timesteps));  % check that 1 is correct
            
            % WEIGHT PREPARATION
            weights_Z = weights.Z_above + weights.Z_below;
            
            % VERTICAL DOWNSCALING
            % T_topoScale = sum(double(era_T) .*  (weights_Z_above + weights_Z_below), 2);
            var3d = sum(double(var3d) .* weights_Z, 2);
            
            % This part needs to change if we don't want the bottom layer to be used
            use_sl = sum(weights_Z, 2) < 1-1e-9   |   alt_point < era_alt_sl;
            merge_w_sl = sum(weights.Z_below, 2) == 0;
            
            % T_topoScale = double(~merge_w_sl) .* T_topoScale + double(merge_w_sl) .* factor .* T_topoScale + double(merge_w_sl) .* (1-factor) .* double(era_T_sl);
            % T_topoScale = double(~use_sl) .* T_topoScale + double(use_sl) .* double(era_T_sl);
            term1 = double(~merge_w_sl) .* var3d;
            term2 = double( merge_w_sl) .* var3d .* weights.factor;
            term3 = double( merge_w_sl) .* var_surf .* (1 - weights.factor);
            
            var3d = double(~use_sl) .* (term1 + term2 + term3) + double(use_sl) .* var_surf;
            
            % HORIZONTAL DOWNSCALING
            var = class.downscale_2D_var(var3d, weights);
            
        end

        function var1d = downscale_2D_var(var2d, weights)
            n_timesteps = size(var2d, ndims(var2d));
            var2d = reshape(var2d, 4, 1, n_timesteps);
            % dimensions are lon x lat, time, where lon x lat are collapsed to n_coords=4
            weights_lon = double(weights.lon(1:2,:,:) + weights.lon(3:4,:,:)) ./ 2;
            weights_lat = double(weights.lat);

            % T_topoScale = T_topoScale .* weights_lat;
            var2d = var2d .* weights_lat;
            % T_topoScale = (T_topoScale(1:2,:,:) +T_topoScale(3:4,:,:));
            var2d = var2d(1:2,:,:) + var2d(3:4,:,:);
            % T_topoScale = squeeze(sum(T_topoScale .* weights_lon,1)) .* era.T_sf;    
            var1d = sum(var2d .* weights_lon,1);
            var1d = squeeze(var1d);

        end

        function wind_speed = calc_wind_speed(u, v)
            wind_speed = sqrt(single(u).^2 + single(v).^2);
        end
    
        function humidity = calc_humidity(temp_C, press)
            class = process_topoScale3;
            
            eps0 = class.epsilon0;
            vpsl = class.magnus_formula(temp_C);
            humidity = eps0 * (vpsl ./ press);
        end
    
        function vap_press_surf = magnus_formula(temp_C)
            tc = double(temp_C);

            check = any(tc > 120);
            if check
                error("Temperatures are not all in degC")
            end

            % AERK from Alduchov and Eskridge (1996).
            A1=17.62; 
            B1=243.12; 
            C1=611.2;  

            vap_press_surf = double(tc>=0) .* C1.*exp(A1.*tc./(B1+tc)) + double(tc<0) .* C1.*exp(22.46.*tc./(272.61+tc)); 
        end
    
        function press_adjusted = barometric_formula(press, press_surf, alt_point, alt_ref0, tempC_ref0)
            class = process_topoScale3;
            g = class.gravity; % Acceleration of gravity [ms^-1]
            R = class.universal_gas_constant;  % Universal Gas Constant (J/(mol*K))
            M = class.molar_mass_of_air; % Molar mass of Earth's air (kg/mol)

            % create a copy before masking
            press_adjusted = press .* 1;
            
            % masking 
            mask = press == 0;
            % creating inputs for barometric formula
            press_surf = press_surf(mask);
            tempK_ref0 = double(tempC_ref0(mask) + 273.15);
            delta_height = alt_point - alt_ref0(mask);
            
            % barometric formula
            press_adjusted(mask) = press_surf .* exp(-g .* M .* delta_height ./ R ./ tempK_ref0);

        end

        function vapour_pressure = calc_vapour_pressure(specific_humidity, atm_pressure)
            class = process_topoScale3;

            epsilon0 = class.epsilon0;

            q2w = @(q) 0.5 .* (1 - sqrt(1 - 4 .* q)); % Mixing ratio from specific humidity based on 2nd order Taylor series expansion.
            wp2e = @(w,p) 0.5 .* p .* (-1 + sqrt(1 + 4 .* w ./ epsilon0)); % Vapor pressure from mixing ratio based on 2nd order Taylor series expansion.
            
            mixing_ratio = q2w(specific_humidity);
            vapour_pressure = wp2e(mixing_ratio, atm_pressure);

        end
        
        function emmissivity = calc_clear_sky_emissivity(vapour_pressure, tempC)
             x1 = 0.43; 
             x2 = 5.7;

             tempK = tempC + 273.15;

             emmissivity = real(0.23 + x1 .* (vapour_pressure ./ tempK).^(1 / x2));
        end
        
        function [mu0, is_sunset] = calc_solar_geometry(solar_zen)
            % Calculates cosine of zenith and determines sunset flag
            mu0 = max(cos(solar_zen), 0); % Ensure non-negative
            is_sunset = mu0 < cosd(89);
        end
        
        function kt = calc_clearness_index(S_surface, S_toa)
            % Ratio of surface to TOA radiation (Clamped to avoid div/0)
            kt = S_surface ./ max(1e-10, S_toa);
        end
        
        function kd = calc_diffuse_fraction(kt)
            % Empirical Model (Likely Skartveit & Olseth)
            % Converts Clearness Index (kt) to Diffuse Fraction (kd)
            
            % The logistic function
            kd = 0.952 - 1.041 .* exp(-1 .* exp(2.3 - 4.702 .* kt));
            
            % Physical Constraints
            kd = max(kd, 0); 
            % kd = min(kd, 1); % Optional safety: fraction can't be > 1
        end
    
    end
    
end

