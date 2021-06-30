#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 23 14:37:36 2021

@author: feng779
"""

from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.collections import PolyCollection

import sys
sys.path.append('.') 
from latlon_coordinate_transforms import fix_periodicity
from myearth import readShpPointLine

import pdb

class mosart_mpaso_grid(object):
    
    rof_grid_on = True
    ocn_grid_on = False
    topolines_on = False 
    usgs_on = False
    xlim = None
    ylim = None
    
    def __init__(self, rof_gridfile, rof_maskfile, ocn_gridfile=None, **kwargs):
        self.__dict__.update(kwargs)
        
        self.rof_gridfile = rof_gridfile
        self.rof_maskfile = rof_maskfile
        self.ocn_gridfile = ocn_gridfile
        
    def read_rof_grid(self):
        """
        function reads the model grid information
        """

        nc = Dataset(self.rof_gridfile, 'r')
        self.lon = nc.variables['longxy'][:]
        self.lat = nc.variables['latixy'][:]
        self.area = nc.variables['area'][:]
        self.dnID = nc.variables['dnID'][:]
        self.ID = nc.variables['ID'][:]
        
        dlon = self.lon[0,1] - self.lon[0,0]
        dlat = self.lat[1,0] - self.lat[0,0]
        
        self.lonp = np.arange(self.lon.min()-dlon/2., self.lon.max()+dlon/2.+dlon, dlon)
        self.latp = np.arange(self.lat.min()-dlat/2., self.lat.max()+dlat/2.+dlat, dlat)
        
    def read_rof_mask(self):
        """
        function reads the model output to get the mask information 
        """
        nc = Dataset(self.rof_maskfile, 'r')
        self.mask = nc.variables['MASK'][0,:,:]
        
    def read_ocn_grid(self):
        """
        function reads the MPASO grid information
        """

        nc = Dataset(self.ocn_gridfile, 'r')
        latCell = nc.variables['latCell'][:]
        lonCell = nc.variables['lonCell'][:]
        xCell = nc.variables['xCell'][:]
        yCell = nc.variables['yCell'][:]
        latVertex = nc.variables['latVertex'][:]
        lonVertex = nc.variables['lonVertex'][:]
        xVertex = nc.variables['xVertex'][:]
        yVertex = nc.variables['yVertex'][:]
        nCells = len(xCell)
        
        nEdgesOnCell = nc.variables['nEdgesOnCell'][:]
        verticesOnCell = nc.variables['verticesOnCell'][:]
        
        vert_lats = np.array([])
        vert_lons = np.array([])
        
        for lat in latVertex:
            vert_lats = np.append(vert_lats, lat * (180 / np.pi)) 
        for lon in lonVertex:
            vert_lons = np.append(vert_lons, lon * (180 / np.pi))
        
        
        # Normalize latitude and longitude
        vert_lons[vert_lons > 180] = vert_lons[vert_lons > 180] - 360
        
        cell_lats = np.array([])
        cell_lons = np.array([])
        
        for lat in latCell:
            cell_lats = np.append(cell_lats, lat * (180 / np.pi)) 
        for lon in lonCell:
            cell_lons = np.append(cell_lons, lon * (180 / np.pi))
        
        # Normalize latitude and longitude
        cell_lons[cell_lons>180] = cell_lons[cell_lons>180] - 360
        
        self.cellxy = []
        for cell in range(nCells):
            vertices = verticesOnCell[cell,:nEdgesOnCell[cell]]
            vertices -= 1
            
            xp = vert_lons[vertices]
            yp = vert_lats[vertices]
            
            xp = fix_periodicity(xp, cell_lons[cell], 360)
            yp = fix_periodicity(yp, cell_lats[cell], 0)
        
            verts = np.vstack([xp, yp]).T
            self.cellxy.append(verts)

        
        
    def grid_plot(self):
        """
        plot grid
        """
        plt.rcParams.update({'font.size': 18}) 
        fig = plt.figure(figsize=(12,10))
        ax = fig.add_subplot(111)
        
        if self.rof_grid_on:
            self.rof_grid_plot(ax)
            
        if self.ocn_grid_on:
            self.ocn_grid_plot(ax)
            
        if self.topolines_on:
            self.topo_lines(ax)
            
        if self.usgs_on:
            self.USGS_stations(ax)
        
        if self.xlim == None or self.ylim == None:
            ax.set_xlim([self.lonp.min(), self.lonp.max()])
            ax.set_ylim([self.latp.min(), self.latp.max()])
        else:
            ax.set_xlim(self.xlim)
            ax.set_ylim(self.ylim)
            
        ax.set_xlabel('Longitude')
        ax.set_ylabel('Latitude')
        ax.set_aspect('equal')
        fig.tight_layout()
        plt.show()
        
        
    def rof_grid_plot(self, ax):
        """
        plot mosart structure grid
        """
        # w = self.lonp[1:] - self.lonp[0:-1]
        # h = self.latp[1:] - self.latp[0:-1]   
        
        # for i, x in enumerate(self.lonp[:-1]):
        #     print(i)
        #     for j, y in enumerate(self.latp[:-1]):
        #         if i % 2 == j % 2: # racing flag style
        #             ax.add_patch( Rectangle((x, y), w[i], h[j], fill=True) )
            
        var = np.ones_like(self.mask)        
        var_masked = np.ma.array(var, mask=self.mask==2)
        cs = ax.pcolor(self.lonp, self.latp, var_masked, facecolor='none', edgecolors='b', linewidths=0.4)
        
    def ocn_grid_plot(self, ax):
        """
        plot mpaso unstructure grid
        """
        self.read_ocn_grid(self.ocn_gridfile)
        
        collection = PolyCollection(self.cellxy,facecolors='none', closed=True)
        collection.set_edgecolors('r')
        collection.set_linewidths(0.4)
        ax.add_collection(collection)
    
    
    def topo_lines(self, ax):
        """
        function that reads topology shapefiles and overlay on the figure
        """
        
        shpfile = '/Users/feng779/OneDrive - PNNL/Documents/DATA/SHP/boundary_lines/ne_10m_coastline/ne_10m_coastline.shp'
        XY, field = readShpPointLine(shpfile)
        
        for line in XY:
            X = line[:,0]
            Y = line[:,1]
            ax.plot(X, Y, '-k', linewidth=0.1)
            
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
            ax.plot(lon_gauge, lat_gauge, 'or', markersize=5)
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
            ax.plot(lon_gauge, lat_gauge, 'ob', markersize=5)
            ax.annotate(stationid, (lon_gauge, lat_gauge), color='k', fontsize=12)
            
    
    
        
if __name__ == "__main__":

    # NLDAS     
    rof_mask_file = '/Users/feng779/OneDrive - PNNL/Documents/CODE/Coupling/NLDAS/files/IELM_NLDAS_NLDAS_yr2002_benchmark.2021-03-30-21.mosart.h0.2004-11-18-00000.nc'
    rof_grid_file = '/Users/feng779/OneDrive - PNNL/Documents/CODE/Coupling/Coupling_input_file/files/MOSART_NLDAS_8th_20210331_coupling.nc'
    ocn_grid_file = None

    # global half-degree
    #rof_mask_file = '/Users/feng779/OneDrive - PNNL/Documents/CODE/Coupling/files/Mosart_output/GSWP/test/datm_GSWP_2way_yr2000_test_2way.2021-03-10-11.mosart.h0.2006-11-08-00000.nc'
    #rof_grid_file = '/Users/feng779/OneDrive - PNNL/Documents/CODE/Coupling/Coupling_input_file/files/MOSART_2way_LLR3_c210324_coupling.nc'
    #ocn_grid_file = '/Users/feng779/OneDrive - PNNL/Documents/CODE/Coupling/grid_EC30to60E2r2/files/ocean.EC30to60E2r2.200908.nc'
    
    MMG = mosart_mpaso_grid(rof_grid_file, rof_mask_file, ocn_grid_file)
    MMG.ocn_grid_on = False
    MMG.topolines_on = True
    MMG.usgs_on = True
    MMG.xlim = [-80, -70]
    MMG.ylim = [32, 42]
    
    MMG.read_rof_grid()
    MMG.read_rof_mask()
    
    MMG.grid_plot()




