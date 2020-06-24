
%%% An example to extract wetland area for a site from global simulated wetland map
%%% Create date: 2020/06/24
%%% by yixi@pku.edu.cn

%% 1: Generate the global wetland map with the optimized parameters

YOURdir = 'F:\7WetlandChange\pro\GLOBAL_201901127\pro_share\';
% please replace it with your directory

WTDdir = [YOURdir, 'water_table_depth\'];
ParaOPTdir = [YOURdir, 'parameter\'];
FWETdir = [YOURdir, 'fwet\'];


% 1: load water table depth data
load([WTDdir, 'rcp85_FGOALS-g2_200608.mat'], 'wt');


% 2: compute wetland fraction with optimal parameters
load([ParaOPTdir, 'GLDAS-Noahv2.0_reso0.5_RFW_max_optPara_RMSE.mat'], 'vkq_fmax_opt');
vp = vkq_fmax_opt(:, :, 1);
kp = vkq_fmax_opt(:, :, 2);
qp = vkq_fmax_opt(:, :, 3);
fmax = vkq_fmax_opt(:, :, 4);
x = (1 + vp.*exp(-kp.*(wt - qp)));
x(x<0) = nan;
x = x.^ (-1 ./ vp);
fmax(isnan(x)) = nan;
fwet_opt = min(x, fmax);


% 3: save
if ~exist(FWETdir, 'dir')
    mkdir(FWETdir)
end
save([FWETdir, 'rcp85_FGOALS-g2_200608.mat'], 'fwet_opt');



%% 2: Extract the wetland area for a site from global wetland map

YOURdir = 'F:\7WetlandChange\pro\GLOBAL_201901127\pro_share\';
% please replace it with your directory


% 1: read shapefile of a site
RAMSARdir = [YOURdir, 'ramsar_shapefile\'];
wet_polygon = shaperead([RAMSARdir, 'ShengjinhLake_ramsarid_2248.shp']);


% 2: create a buffer
wet_polygon_buffer = wet_polygon;
buffersize = 5; % times of site area
% vertices of initial site
poly = polyshape(wet_polygon.X, wet_polygon.Y);
x = wet_polygon.X; % logitude
y = wet_polygon.Y; % latitude
ind_nan = find(isnan(x)|isnan(y));
x(ind_nan) = [];
y(ind_nan) = [];
% area of the initial site
s_site = polyarea(x, y);
% create buffer
buffer_distance = (sqrt(buffersize)-1)*sqrt(s_site/pi);
poly_buffer = polybuffer(poly, buffer_distance);
% save
wet_polygon_buffer.X = poly_buffer.Vertices(:, 1);
wet_polygon_buffer.Y = poly_buffer.Vertices(:, 2);
clear wet_polygon x y s_site buffer_distance poly poly_buffer ind_nan


% 3: extract grids within the buffered site boundary
reso = 0.5;
[nb_lat, nb_lon] = deal(180/reso, 360/reso);
lat = -90+reso/2:reso:90;
lon = -180+reso/2:reso:180;

rg = wet_polygon_buffer.BoundingBox; % boundary
x = wet_polygon_buffer.X;
y = wet_polygon_buffer.Y;

poly = polyshape(x, y);

Ln = nan(1, 1);
Col = nan(1, 1);
n = 0;
for la = lat(find(lat<=(rg(1, 2)), 1, 'last')):reso:lat(find(lat>=(rg(2, 2)), 1))
    for lo = lon(find(lon<=(rg(1, 1)), 1, 'last')):reso:lon(find(lon>=(rg(2, 1)), 1))
        
        x_grid = [lo-reso/2, lo+reso/2, lo+reso/2, lo-reso/2, lo-reso/2];
        y_grid = [la-reso/2, la-reso/2, la+reso/2, la+reso/2, la-reso/2];
        
        poly2 = polyshape(x_grid, y_grid);
        % the intersect of site boundary and grid
        poly_intersect = intersect(poly, poly2);
        
        if isempty(poly_intersect.Vertices)
            continue;
        else
            ln = round((lat(1, end)-la)/reso)+1;
            col = round((lo-lon(1, 1))/reso)+1;
            n = n+1;
            Ln(n, 1) = ln;
            Col(n, 1) = col;
        end
        
    end
end
clear la* lo* rg x* y* poly* n ln col


% 4: extract wetland area for the site
FWETdir = [YOURdir, 'fwet\'];
load([FWETdir, 'rcp85_FGOALS-g2_200608.mat'], 'fwet_opt');

s = area_weighted(nb_lat, nb_lon)/1e6; % km2

CTIfile = 'F:\DATA\TOPMODEL\ga2\ga2.nc';
% Please replace it with your directory
% The CTI file can be available from Marthews et al. (2015)
% https://catalogue.ceh.ac.uk/documents/6b0c4358-2bf3-4924-aa8f-793d468b92be

reso_CTI = 1/240;
lat = -90+reso_CTI/2:reso_CTI:90;
lon = -180+reso_CTI/2:reso_CTI:180;
no_lat = [8076; 42264]; % start and end

tt = reso/reso_CTI;

nb_gd = length(Ln);
fwet_ramsar = nan(nb_gd, 1);
for igd = 1:nb_gd
    ln = Ln(igd);
    col = Col(igd);
    if ~isnan(fwet_opt(ln, col))
        % extract CTI value for the grid
        CTIx = ncread(CTIfile, 'Band1', [(col-1)*tt+1 length(lat)-ln*tt+1-no_lat(1)+1], [tt tt]);
        CTIx = rot90(CTIx, 1);
        % sort CTI of all subgrids
        [~, I] = sort(CTIx(:), 'descend', 'MissingPlacement', 'last');
        CTIx_no = nan(size(CTIx));
        CTIx_no(I) = 1:(tt*tt);
        CTIx_no(isnan(CTIx)) = nan;
        
        lat_sub = lat(length(lat)-(ln-1)*tt):-reso_CTI:lat(length(lat)-ln*tt+1);
        lon_sub = lon((col-1)*tt+1):reso_CTI:lon(col*tt);
        la_ln = 0;
        n = 0;
        CTIsub = nan(1, 1);
        for la = lat_sub
            la_ln = la_ln+1;
            lo_ln = 0;
            for lo = lon_sub
                lo_ln = lo_ln+1;
                % determine if the subgrid is in/on the site boundary
                [isINCLUDE, isONEDGE] = inpolygon(la, lo, ...
                    wet_polygon_buffer.Y, wet_polygon_buffer.X);
                if isINCLUDE==1 || isONEDGE==1
                    n = n+1;
                    CTIsub(n) = CTIx_no(la_ln, lo_ln);
                end
            end
        end

        CTIsub = CTIsub';
        
        % determine the inundated subgrids with CTI
        fwet_ramsar(igd, 1) = ...
            length(find(CTIsub<=(tt*tt)*fwet_opt(ln, col)))/(tt*tt)*s(ln, col);
    end
end
clear I CTI* la* lo* ln col


