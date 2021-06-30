import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import MFDataset, Dataset, num2date

import sys
sys.path.append('.')
import othertime
from myearth import findNearset1D, findNearset2D

import pdb

river_names = {'01463500': 'Delaware', '01578310': 'Susquehanna', \
               '01646502': 'Potomac',  '01594440': 'Patuxeat', \
               '01491000': 'Choptank', '01668000': 'Rappahannock', \
               '01674500': 'Mattaponi','01673000': 'Pamunkey', \
               '02035000': 'James',    '02041650': 'Appomattox', \
               '01554000': 'Susquehanna'}


def readUSGS(stationid_inp):
    """
    read USGS discharge data at a station
    """
    filename_USGS = '/qfs/people/feng779/DATA/USGS/USGS_streamflow.nc'
    nc = Dataset(filename_USGS, 'r')
    for g in nc.groups:
        tmp_grp = nc.groups[g]
        tmp_Q = tmp_grp.variables['discharge']
        stationid = tmp_Q.getncattr('StationID')
        if stationid == stationid_inp:
            lon_gauge = tmp_grp.variables['longitude'][:]
            lat_gauge = tmp_grp.variables['latitude'][:]
            if stationid_inp == '01428500':
                lon_gauge = -74.936
                lat_gauge = 41.443
            elif stationid_inp == '01673000':
                lon_gauge = -77.31
                lat_gauge = 37.68
            timei_data = tmp_grp.variables['time']
            time_data = num2date(timei_data[:], timei_data.units, only_use_cftime_datetimes=False)
            discharge = tmp_grp.variables['discharge'][:,0,0,0]
    return time_data, discharge, lon_gauge, lat_gauge

def read_model_output(filename, vname):

    if len(filename) == 1:
        MFnc = Dataset(filename, 'r')
    else:
        MFnc = MFDataset(filename, 'r')

    lon = MFnc.variables['lon'][:]
    lat = MFnc.variables['lat'][:]
    vv = MFnc.variables[vname][:]
    MASK = MFnc.variables['MASK'][:]
    timei_model = MFnc.variables['time']
    time_model = num2date(timei_model[:], timei_model.units, only_use_cftime_datetimes=False)

    return vv, time_model, lon, lat

def plot_hurricane_event(stationids_inp, filenames, starttimes, endtimes, figpaths):
    """
    plot validation time series given the hurricane events
    """
   
    for k in range(len(filenames)):
        filename = filenames[k]
        starttime = starttimes[k]
        endtime = endtimes[k]
        
        plt.rcParams.update({'font.size': 14})
        if len(stationids_inp) <= 4:
            fig, axes = plt.subplots(len(stationids_inp),1, figsize=(15, 12))
        else:
            fig, axes = plt.subplots( int(np.floor(len(stationids_inp)/3.))+1,3, figsize=(15, 12))
        axes = axes.ravel()

        for i, stationid_inp in enumerate(stationids_inp):
            ## data
            time_data, discharge, lon_gauge, lat_gauge = readUSGS(stationid_inp)
            ind0 = othertime.findNearest(starttime, time_data)
            ind1 = othertime.findNearest(endtime, time_data)
            time_data, discharge = time_data[ind0:ind1+1], discharge[ind0:ind1+1]
            axes[i].plot(time_data, discharge, 'ok', label='USGS')
            ## model 
            vv, time_model, lon, lat = read_model_output(filename, 'RIVER_DISCHARGE_OVER_LAND_LIQ')
            ind2 = othertime.findNearest(starttime, time_model)
            ind3 = othertime.findNearest(endtime, time_model)
            if vv.ndim == 2:
                ind = findNearset2D(lon_gauge, lat_gauge, lon, lat)
                axes[i].plot(time_model[ind2:ind3+1], vv[ind2:ind3+1,ind], label='model')
            elif vv.ndim == 3:
                ind_lon = findNearset1D(lon_gauge, lon)[0][0]
                ind_lat = findNearset1D(lat_gauge, lat)[0][0] 
                axes[i].plot(time_model[ind2:ind3+1], vv[ind2:ind3+1, ind_lat, ind_lon], label='model')

            axes[i].set_ylabel('Streamflow [m3/s]')
            axes[i].set_title('{}'.format(river_names[stationid_inp]))
            axes[i].grid() 
        if len(stationids_inp) > 4:
            for j in range(len(stationids_inp), (int(np.floor(len(stationids_inp)/3.))+1)*3):
                axes[j].set_visible(False)
 
        plt.gcf().autofmt_xdate()
        plt.tight_layout()
        plt.savefig(figpaths[k])
        plt.close()
