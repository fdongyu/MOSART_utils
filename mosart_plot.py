import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import matplotlib
from matplotlib.path import Path
from matplotlib.patches import PathPatch
from netCDF4 import Dataset
from osgeo import osr
import ogr

import sys
sys.path.append('.') 
from myearth import findNearset1D

import pdb

class mosart_plot(object):
    """
    general class for visualizing MOSART model result
    """
    
    figname = None
    scale_factor = 5 
    bounds = None 
    cmap = 'bone_r'
    mask_on = True 
    river_network_on = False 
    topolines_on = False 
    clip_on = False 
    outlet_pts_on = False 
    subset = False
    usgs_on = False
    flood_fraction_on = False
    multi_river_basin_on = False
    xlim = None
    ylim = None
    vmin = None
    vmax = None
    

    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)
        

    def readgrid(self,gridfile):
        """
        function reads the model grid information
        """

        nc = Dataset(gridfile, 'r')
        self.lon = nc.variables['longxy'][:]
        self.lat = nc.variables['latixy'][:]
        self.area = nc.variables['area'][:]
        self.dnID = nc.variables['dnID'][:]
        self.ID = nc.variables['ID'][:]
    
    def river_network(self, ax, line_color='gray', line_width=0.5):
        """
        This function provides the option to overlay the river network on the plot 
        """
        def show_river_2d(dnID,ID,longxy,latixy):
            [m,n] = dnID.shape
            for i in range(m):
                for j in range(n):
                    if dnID[i,j] == -9999:
                        continue
                    else:
                        [i2,j2] = np.concatenate(np.where(ID == dnID[i,j]))
                        ax.plot([longxy[i,j],longxy[i2,j2]], [latixy[i,i], latixy[i2,j2]], linestyle='-', color=line_color, linewidth=line_width)

        def show_river_1d(dnID,ID,longxy,latixy):
            m = len(dnID);
            for i in range(m):
                if dnID[i] == -9999:
                    continue
                else:
                    i2 = np.concatenate(np.where(ID == dnID[i]))
                    ax.plot([longxy[i], longxy[i2]], [latixy[i], latixy[i2]], linestyle='-', color=line_color, linewidth=line_width)

        if self.dnID.ndim == 1:
            show_river_1d(self.dnID, self.ID, self.lon, self.lat)
        elif self.dnID.ndim == 2:
            show_river_2d(self.dnID, self.ID, self.lon, self.lat)

    def contour_map(self, var, varname, MASK):
        """
        function generates the contour map
        """
        

        if self.cmap == 'bwr':
            vmin = (-1)*max(abs(var.min()), abs(var.max()))
            vmax = max(abs(var.min()), abs(var.max()))
            levels = np.linspace(vmin/self.scale_factor, vmax/self.scale_factor, 100)
        else:
            vmin = var.min()
            vmax = var.max()
            levels = np.linspace(vmin/self.scale_factor, vmax/self.scale_factor, 100)            
        if self.vmin != None and self.vmax != None:
            vmin = self.vmin
            vmax = self.vmax
        
        if type(self.cmap) == str:   ## not self defined colormap
            cmap = plt.set_cmap(self.cmap)
        else:
            plt.register_cmap(None, self.cmap)
            cmap = self.cmap
        
        plt.rcParams.update({'font.size': 18}) 
        fig = plt.figure(figsize=(12,10))
        ax = fig.add_subplot(111)
        
        if self.mask_on:
            #cs = self.mask_plot(ax, var, MASK, cmap, levels)
            cs = self.mask_plot_patch(ax, var, MASK, cmap, vmax/self.scale_factor, vmin/self.scale_factor)
        else:
            #cs = self.unmask_plot(ax, var, cmap, levels)
            cs = self.unmask_plot_patch(ax, var, cmap, vmax/self.scale_factor, vmin/self.scale_factor)


        from mpl_toolkits.axes_grid1 import make_axes_locatable
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="3%", pad=0.05)
        if type(self.bounds) == type(None):
            cb = fig.colorbar(cs, cax=cax, orientation='vertical')
        else:
            cb = fig.colorbar(cs, cax=cax, ticks=self.bounds, orientation='vertical')
        cb.ax.tick_params(labelsize=12)
        cb.ax.yaxis.offsetText.set_fontsize(12)
        cb.set_label(varname, fontsize=14)


        if self.river_network_on:
            self.river_network(ax, 'r') 

        if self.topolines_on:
            self.topo_lines(ax)
            
        if self.multi_river_basin_on:
            self.multi_river_basin(ax)
            
        if self.clip_on:
            #clip_shp = '/Users/feng779/OneDrive - PNNL/Documents/DATA/SHP/ICoM_domain/srb_reprojected/srb_reprojected.shp'
            clip_shp = '/Users/feng779/OneDrive - PNNL/Documents/DATA/SHP/ICoM_domain/Chesapeake_Bay_Watershed_Boundary/Chesapeake_Bay_Watershed_Boundary.shp'
            clip_shp2 = '/Users/feng779/OneDrive - PNNL/Documents/DATA/SHP/ICoM_domain/drbbnd_reprojected/drb_bnd_polygon_reproj.shp'
            #clip = self.shp_clip(ax, clip_shp)
            clip = self.shp_clip2(ax, clip_shp, clip_shp2)
            if self.flood_fraction_on:
                ff = self.ff_calc(ax, clip_shp, clip_shp2, var)
            #pdb.set_trace()
            try: 
                for contour in cs.collections:
                    contour.set_clip_path(clip)
            except:
                cs.set_clip_path(clip)

        if self.outlet_pts_on:
            self.outlet_pts(ax, MASK, self.subset)
            
        if self.usgs_on:
            self.USGS_stations(ax)

        ax.set_xlabel('Longitude')
        ax.set_ylabel('Latitude')
        if self.xlim == None or self.ylim == None:
            ax.set_xlim([self.lon.min(), self.lon.max()])
            ax.set_ylim([self.lat.min(), self.lat.max()])
        else:
            ax.set_xlim(self.xlim)
            ax.set_ylim(self.ylim)
        
        ax.set_aspect('equal')
        fig.tight_layout()
        #plt.savefig(figname)
        #plt.close()
        if self.figname != None:
            plt.savefig(self.figname)
            plt.close()
        else:
            plt.show()
        
    def mask_plot(self, ax, var, MASK, cmap, levels):
        """
        plot with mask enabled 
        """
        if var.ndim == 1:
            triang = tri.Triangulation(self.lon[MASK==1], self.lat[MASK==1])            
            cs = ax.tricontourf(triang, var[MASK==1], cmap=cmap, levels=levels, extend='max')
        elif var.ndim == 2:
            var_masked = np.ma.array(var, mask=MASK==2)
            if type(self.bounds) == type(None):
                cs = ax.contourf(self.lon, self.lat, var_masked, cmap=self.cmap, levels=levels, extend='max')
            else:
                #c_nodata = '#ffffff'
                c1 = '#4679fa'
                c20 = '#46befa'
                c40 = '#61fa46'
                c60 = '#e2fa46'
                c80 = '#faa346'
                c100 = '#fa4c46'
                cmap = matplotlib.colors.ListedColormap([c1, c20, c40, c60, c80, c100])
                cs = ax.contourf(self.lon, self.lat, var_masked, self.bounds, cmap=cmap, vmin=0, vmax=self.bounds[-1])
        return cs
    
    def unmask_plot(self, ax, var, cmap, levels):
        """
        plot without mask enabled 
        """
        if var.ndim == 1:
            triang = tri.Triangulation(self.lon, self.lat)            
            cs = ax.tricontourf(triang, var, cmap=cmap, levels=levels, extend='both')
        elif var.ndim == 2:
            if type(self.bounds) == type(None):
                cs = ax.contourf(self.lon, self.lat, var, cmap=self.cmap, levels=levels, extend='both')
            else:
                #c_nodata = '#ffffff'
                c1 = '#4679fa'
                c20 = '#46befa'
                c40 = '#61fa46'
                c60 = '#e2fa46'
                c80 = '#faa346'
                c100 = '#fa4c46'
                cmap = matplotlib.colors.ListedColormap([c1, c20, c40, c60, c80, c100])
                cs = ax.contourf(self.lon, self.lat, var, self.bounds, cmap=cmap, vmin=0, vmax=self.bounds[-1])
        return cs
    
    def mask_plot_patch(self, ax, var, MASK, cmap, vmax, vmin):
        """
        plot using pcolor function
        """
        
        #lonp = np.arange(-180, 180+0.5, 0.5) # This is only for half degree grid
        #latp = np.arange(-90, 90+0.5, 0.5)
        
        dlon = self.lon[0,1] - self.lon[0,0]
        dlat = self.lat[1,0] - self.lat[0,0]
        
        lonp = np.arange(self.lon.min()-dlon/2., self.lon.max()+dlon/2.+dlon, dlon)
        latp = np.arange(self.lat.min()-dlat/2., self.lat.max()+dlat/2.+dlat, dlat)
        
        var_masked = np.ma.array(var, mask=MASK==2)
        cs = ax.pcolor(lonp, latp, var_masked, cmap=cmap, edgecolors=None, linewidths=0.1, \
                       vmin=vmin, vmax=vmax)
        return cs
        
    
    def unmask_plot_patch(self, ax, var, cmap, vmax, vmin):
        """
        plot using pcolor function
        """
        from matplotlib.patches import Polygon
        from matplotlib.collections import PatchCollection
        
        #lonp = np.arange(-180, 180+0.5, 0.5) # This is only for half degree grid
        #latp = np.arange(-90, 90+0.5, 0.5)
        
        dlon = self.lon[0,1] - self.lon[0,0]
        dlat = self.lat[1,0] - self.lat[0,0]
        
        lonp = np.arange(self.lon.min()-dlon/2., self.lon.max()+dlon/2.+dlon, dlon)
        latp = np.arange(self.lat.min()-dlat/2., self.lat.max()+dlat/2.+dlat, dlat)
        
        # dlon = 0.25
        # dlat = 0.25
        # patches = []
        # for i, llon in enumerate(self.lon):
        #     for j, llat in enumerate(self.lat):
        #         lon_w = llon - dlon
        #         lon_e = llon + dlon
        #         lat_n = llat + dlat
        #         lat_s = llat - dlat          
        #         ## counter-clockwise vertices for each cell
        #         xp = np.asarray([lon_w, lon_e, lon_e, lon_w])
        #         yp = np.asarray([lat_s, lat_n, lat_n, lat_s])
        #         patches.append(Polygon(np.vstack([xp, yp]).T))
        
        # cs = PatchCollection(patches, cmap=cmap)
        # cs.set_array(var.flatten())
        # cs.set_lw(0.1)
        # ax.add_collection(cs)
        cs = ax.pcolor(lonp, latp, var, cmap=cmap, edgecolors=None, linewidths=0.1, \
                       vmin=vmin, vmax=vmax)
        return cs
        
    
    def ff_calc(self, ax, clip_shp, clip_shp2, var):
        """
        calculate flood fraction
        """
        
        XY1, field1  = readShpPoly(clip_shp)
        XY2, field2  = readShpPoly(clip_shp2)
        XY = [XY1[0], XY2[0]]
        
        lon, lat, FF = self.subset_domain(self.lon, self.lat, var)
        lon = lon.flatten()
        lat = lat.flatten()
        FF = FF.flatten()
        
        
        ff = []
        for i in range(2):
            vertices = []
            codes = []
            for j in range(len(XY[i])):
                vertices.append((XY[i][j][0], XY[i][j][1]))
            codes += [Path.MOVETO]
            codes += [Path.LINETO] * (len(XY[i]) -2)
            codes += [Path.CLOSEPOLY]
        
        
            clip = Path(vertices, codes)
            flags = clip.contains_points(np.vstack((lon, lat)).T)
            
            FF_inpoly = FF[flags]
            ff_tem = len(FF_inpoly[FF_inpoly>1e-3]) / len(FF_inpoly) 
            ff.append(ff_tem)
        
        ## write flood fraction
        #ax.text(0.01, 0.92, 'Susquehanna: {:.2f}'.format(ff[0]), fontsize=15, transform=ax.transAxes)
        #ax.text(0.01, 0.92, 'Chesapeake: {:.2f}'.format(ff[0]), fontsize=15, transform=ax.transAxes)
        ax.text(0.01, 0.88, 'Delaware:    {:.2f}'.format(ff[1]), fontsize=15, transform=ax.transAxes)
        
    

    def topo_lines(self, ax):
        """
        function that reads topology shapefiles and overlay on the figure
        """
#        import cartopy.io.shapereader as shpreader                   
#        import sys
#        sys.path.append('../utils') #
#        from maptools import readShpPointLine
        shpfile = '/qfs/people/feng779/DATA/SHP/boundary_lines/ne_10m_coastline/ne_10m_coastline.shp'
        XY, field = readShpPointLine(shpfile)
        
        for line in XY:
            X = line[:,0]
            Y = line[:,1]
            ax.plot(X, Y, '-k', linewidth=0.1)
            
    
    def multi_river_basin(self, ax):
        """
        multiple river basin boundaries
        """
        
        shp_delaware = '/Users/feng779/OneDrive - PNNL/Documents/DATA/SHP/ICoM_domain/drbbnd_reprojected/drb_bnd_polygon_reproj.shp'
        shp_susquehanna = '/Users/feng779/OneDrive - PNNL/Documents/DATA/SHP/ICoM_domain/srb_reprojected/srb_reprojected.shp'
        shp_potomac = '/Users/feng779/OneDrive - PNNL/Documents/DATA/SHP/chesapeake-bay-data/Potomac_river_basin/Potomac_river_basin.shp'
        shp_patuxent = '/Users/feng779/OneDrive - PNNL/Documents/DATA/SHP/chesapeake-bay-data/Patuxent_river_basin/Patuxent_river_basin.shp'
        shp_choptank = '/Users/feng779/OneDrive - PNNL/Documents/DATA/SHP/chesapeake-bay-data/Choptank_river_basin/Choptank_river_basin.shp'
        shp_rappahannock = '/Users/feng779/OneDrive - PNNL/Documents/DATA/SHP/chesapeake-bay-data/Rappahannock_river_basin/Rappahannock_river_basin.shp'
        shp_mattaponi_pamunkey = '/Users/feng779/OneDrive - PNNL/Documents/DATA/SHP/chesapeake-bay-data/Mattaponi_Pamunkey_river_basin/Mattaponi_Pamunkey_river_basin.shp'
        shp_james_appomattox = '/Users/feng779/OneDrive - PNNL/Documents/DATA/SHP/chesapeake-bay-data/James_river_basin/James_river_basin.shp'
        
        shp_susquehanna_bay = '/Users/feng779/OneDrive - PNNL/Documents/DATA/SHP/ICoM_domain/clipped_bay_boundary/Chesapeake_bay_poly.shp'
        shp_delaware_bay = '/Users/feng779/OneDrive - PNNL/Documents/DATA/SHP/ICoM_domain/clipped_bay_boundary/Delaware_bay_poly.shp'
    
        self.basin_bound(ax, shp_susquehanna)
        self.basin_bound(ax, shp_potomac)
        self.basin_bound(ax, shp_patuxent)
        self.basin_bound(ax, shp_choptank)
        self.basin_bound(ax, shp_rappahannock)
        self.basin_bound(ax, shp_mattaponi_pamunkey)
        self.basin_bound(ax, shp_james_appomattox)
        
        self.basin_bound(ax, shp_susquehanna_bay)
        self.basin_bound(ax, shp_delaware_bay)
        
    
    def basin_bound(self, ax, shpfile):
        
        XY, field = readShpPoly(shpfile)
        
        for line in XY:
            X = line[:,0]
            Y = line[:,1]
            ax.plot(X, Y, '-k', linewidth=0.3)

    def outlet_pts(self, ax, mask, subset):
        """
        function that plot basin outlet points
        """

        if subset:
            ind_lon0 = findNearset1D(-77.5, self.lon[0,:])[0][0]
            ind_lon1 = findNearset1D(-75, self.lon[0,:])[0][0]
            ind_lat0 = findNearset1D(35, self.lat[:,0])[0][0]
            ind_lat1 = findNearset1D(40, self.lat[:,0])[0][0]
            lon_sub = self.lon[ind_lat0:ind_lat1+1, ind_lon0:ind_lon1+1]
            lat_sub = self.lat[ind_lat0:ind_lat1+1, ind_lon0:ind_lon1+1]
            mask_sub = mask[ind_lat0:ind_lat1+1, ind_lon0:ind_lon1+1]
            pdb.set_trace()
            ax.plot(lon_sub[mask_sub[:,:]==3], lat_sub[mask_sub[:,:]==3], '.k', markersize=2)
        else:
           ax.plot(self.lon[mask[:,:]==3], self.lat[mask[:,:]==3], '.k', markersize=1) 
           
    def shp_clip(self, ax, clip_shp):
        """
        create clip based on a shapefile 
        Ex:
            https://basemaptutorial.readthedocs.io/en/latest/clip.html
            https://stackoverflow.com/questions/55483158/how-to-clip-contours-inside-of-shape-file-using-matplotlib
        """
        #shp = '/Users/feng779/OneDrive - PNNL/Documents/DATA/SHP/ICoM_domain/srb/srb.shp'
        ## reprojected in QGIS, set coordinate system to WGS84 EPSG(4326)
        ## use reproject layer tool
        ## https://gis.stackexchange.com/questions/35590/reprojecting-vector-layer-in-qgis
        
        XY, field = readShpPoly(clip_shp)
        
        # for line in XY[0]:
        #     X = line[:,0]
        #     Y = line[:,1]
        #     ax.plot(X, Y, '-r', linewidth=0.5)
        X = XY[0][:,0]
        Y = XY[0][:,1]
        ax.plot(X, Y, '-k', linewidth=1.0)
        
        vertices = np.asarray(XY[0])
        codes = [Path.MOVETO] + [Path.LINETO]*(vertices.shape[0]-2) + [Path.CLOSEPOLY]
        
        clip = Path(vertices, codes)
        clip = PathPatch(clip, transform=ax.transData)
        #clip = PathPatch(clip, facecolor = 'white')
        #pdb.set_trace()
        return clip
    
        
    def shp_clip2(self, ax, clip_shp, clip_shp2):
        """
        create clip based on multi shapefile 
        """
        XY1, field1  = readShpPoly(clip_shp)
        XY2, field2 = readShpPoly(clip_shp2)
        #pdb.set_trace()
        XY = [XY1[0], XY2[0]]
        
        X = XY1[0][:,0]
        Y = XY1[0][:,1]
        ax.plot(X, Y, '-k', linewidth=1.0)
            
        X = XY2[0][:,0]
        Y = XY2[0][:,1]
        ax.plot(X, Y, '-k', linewidth=1.0)
        
        vertices = []
        codes = []
        for i in range(2):
            for j in range(len(XY[i])):
                vertices.append((XY[i][j][0], XY[i][j][1]))
                #pdb.set_trace()
            codes += [Path.MOVETO]
            codes += [Path.LINETO] * (len(XY[i]) -2)
            codes += [Path.CLOSEPOLY]
        
        #vertices = np.asarray(XY[0])
        #codes = [Path.MOVETO] + [Path.LINETO]*(vertices.shape[0]-2) + [Path.CLOSEPOLY]
        
        clip = Path(vertices, codes)
        clip = PathPatch(clip, transform=ax.transData)
        #clip = PathPatch(clip, facecolor = 'white')
        #pdb.set_trace()
        return clip
    
    def USGS_stations(self, ax):
        """
        plot USGS gauges
        """
        
        ## USGS streamflow
        filename = '/Users/feng779/OneDrive - PNNL/Documents/DATA/USGS/USGS_streamflow.nc'
        nc = Dataset(filename, 'r')
        for g in nc.groups:
            tmp_grp = nc.groups[g]
            tmp_Q = tmp_grp.variables['discharge']
            stationid = tmp_Q.getncattr('StationID')
            lon_gauge = tmp_grp.variables['longitude'][:]
            lat_gauge = tmp_grp.variables['latitude'][:]
            ax.plot(lon_gauge, lat_gauge, 'ok', markersize=5)
            ax.annotate(stationid, (lon_gauge, lat_gauge), color='k', fontsize=12)
        
        ## USGS gauge height
        filename = '/Users/feng779/OneDrive - PNNL/Documents/DATA/USGS/USGS_gauge_height.nc'
        nc = Dataset(filename, 'r')
        for g in nc.groups:
            tmp_grp = nc.groups[g]
            tmp_Q = tmp_grp.variables['gauge_height']
            stationid = tmp_Q.getncattr('StationID')
            lon_gauge = tmp_grp.variables['longitude'][:]
            lat_gauge = tmp_grp.variables['latitude'][:]
            ax.plot(lon_gauge, lat_gauge, 'ok', markersize=5)
            ax.annotate(stationid, (lon_gauge, lat_gauge), color='k', fontsize=12)
            
    def subset_domain(self, lonc, latc, data):
        """
        subset a tile for fast processing
        """
        ## bbox of the subset region
        bbox = [-79, -74.3, 38.6, 43]
        if self.lon.ndim > 1:
            ind_lon0 = findNearset1D(bbox[0], lonc[0,:])
            ind_lon1 = findNearset1D(bbox[1], lonc[0,:])
            ind_lat0 = findNearset1D(bbox[2], latc[:,0])
            ind_lat1 = findNearset1D(bbox[3], latc[:,0])
            
            return lonc[ind_lat0:ind_lat1+1, ind_lon0:ind_lon1+1], latc[ind_lat0:ind_lat1+1, ind_lon0:ind_lon1+1], \
                data[ind_lat0:ind_lat1+1, ind_lon0:ind_lon1+1]
            
        else:
            ind_lon0 = findNearset1D(bbox[0], lonc)
            ind_lon1 = findNearset1D(bbox[1], lonc)
            ind_lat0 = findNearset1D(bbox[2], latc)
            ind_lat1 = findNearset1D(bbox[3], latc)
        
            return lonc[ind_lon0:ind_lon1+1], latc[ind_lat0:ind_lat1+1], \
                data[ind_lat0:ind_lat1+1, ind_lon0:ind_lon1+1]
                
                
        

def readShpPointLine(shpfile,FIELDNAME=None):
    """ Reads a shapefile with line or point geometry and returns x,y,z
    
    See this tutorial:
        http://www.gis.usu.edu/~chrisg/python/2009/lectures/ospy_slides1.pdf
    """
    

    # Open the shapefile
    driver = ogr.GetDriverByName('ESRI Shapefile')
    shp = driver.Open(shpfile, 0)
    #pdb.set_trace()
    #shp = ogr.Open(shpfile)
    
    lyr = shp.GetLayer()
    
    lyr.ResetReading()
    XY=[]
    field=[]
    for feat in lyr:
        feat_defn = lyr.GetLayerDefn()
        for i in range(feat_defn.GetFieldCount()):
            field_defn = feat_defn.GetFieldDefn(i)

            if FIELDNAME==None:
                geom = feat.GetGeometryRef()
                #ztmp = float(feat.GetField(i))
                if geom.GetGeometryType() == ogr.wkbPoint: # point
                    #field.append(ztmp)
                    XY.append([geom.getX(),geom.getY()])
                elif geom.GetGeometryType() == 2:  # line
                    xyall=geom.GetPoints()  
                    XY.append(np.asarray(xyall))
                    #field.append(ztmp)
                        
                elif geom.GetGeometryType() == 5:  # multiline
                    for ii in range(0,geom.GetGeometryCount()):
                        geom2 = geom.GetGeometryRef(ii)
                        xyall=geom2.GetPoints()
                        XY.append(np.asarray(xyall))
                        #field.append(ztmp)

            elif field_defn.GetName() == FIELDNAME:
                geom = feat.GetGeometryRef()
                #ztmp = float(feat.GetField(i))
                ztmp = feat.GetField(i)
                if geom.GetGeometryType() == ogr.wkbPoint: # point
                    field.append(ztmp)
                    #XY.append([geom.getX(),geom.getY()])
                    XY.append(geom.GetPoints())
                elif geom.GetGeometryType() == 2:  # line
                    xyall=geom.GetPoints()
                    XY.append(np.asarray(xyall))
                    field.append(ztmp)
                        
                elif geom.GetGeometryType() == 5:  # multiline
                    for ii in range(0,geom.GetGeometryCount()):
                        geom2 = geom.GetGeometryRef(ii)
                        xyall=geom2.GetPoints()
                        XY.append(np.asarray(xyall))
                        field.append(ztmp)

    shp=None
    
    return XY,field    

def readShpPoly(shpfile,FIELDNAME = None):
    """ Reads a shapefile with polygon geometry and returns x,y and FIELDNAME value
    
    See this tutorial:
        http://www.gis.usu.edu/~chrisg/python/2009/lectures/ospy_slides1.pdf
    """
    # Open the shapefile
    driver = ogr.GetDriverByName('ESRI Shapefile')
    
    shp = driver.Open(shpfile, 0)
    
    lyr = shp.GetLayer()
    
    lyr.ResetReading()
    XY=[]
    field=[]
    for feat in lyr:
        feat_defn = lyr.GetLayerDefn()
        for i in range(feat_defn.GetFieldCount()):
            field_defn = feat_defn.GetFieldDefn(i)
            if FIELDNAME==None:
                # Get all of the polygons
                geom = feat.GetGeometryRef()
                ztmp = feat.GetField(i)
               
                if geom.GetGeometryType() == ogr.wkbPolygon:  # Polygon
                
                    for ii in range(0,geom.GetGeometryCount()):
                        geom2 = geom.GetGeometryRef(ii)
                        xyall=geom2.GetPoints()
                        
                        XY.append(np.asarray(xyall))
                        
                if geom.GetGeometryType() == ogr.wkbMultiPolygon:  # Multi Polygon
                    for ii in range(0,geom.GetGeometryCount()):
                        geom2 = geom.GetGeometryRef(ii)
                        for jj in range(0,geom2.GetGeometryCount()):
                            geom3 = geom2.GetGeometryRef(jj)
                            xyall=geom3.GetPoints()
                            
                            XY.append(np.asarray(xyall))
                
                
            
            elif field_defn.GetName() == FIELDNAME:
                geom = feat.GetGeometryRef()
                ztmp = feat.GetField(i)
               
                if geom.GetGeometryType() == ogr.wkbPolygon:  # Polygon
                
                    for ii in range(0,geom.GetGeometryCount()):
                        geom2 = geom.GetGeometryRef(ii)
                        xyall=geom2.GetPoints()
                        
                        XY.append(np.asarray(xyall))
                        field.append(ztmp)
                        
                if geom.GetGeometryType() == ogr.wkbMultiPolygon:  # Multi Polygon
                    for ii in range(0,geom.GetGeometryCount()):
                        geom2 = geom.GetGeometryRef(ii)
                        for jj in range(0,geom2.GetGeometryCount()):
                            geom3 = geom2.GetGeometryRef(jj)
                            xyall=geom3.GetPoints()
                            
                            XY.append(np.asarray(xyall))
                            field.append(ztmp)  
    
    shp=None
    return XY,field
