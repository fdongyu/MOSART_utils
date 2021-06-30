import numpy as np
from datetime import datetime

def time_convert(intime):
    """
    function to convert the time from string to datetime
    """
    Nt = intime.shape[0]
    outtime = []
    for t in range(Nt):
        timestr = ''.join([intime[t,:][~intime[t,:].mask].data[i].decode('utf-8') for i in range(len(intime[t,:][~intime[t,:].mask].data))])
        outtime.append(datetime.strptime(timestr, '%Y-%m-%d_%H:%M:%S'))
    return outtime

def findNearset(x,y,lon,lat):
    """
    Return the J,I indices of the nearst grid cell to x,y
    """
    dist = np.sqrt( (lon - x)**2 + (lat - y)**2)

    return np.argwhere(dist==dist.min())[0][0]

def convertLatLon(latCell, lonCell):
    """
    convert MPASO coordinates to latitude and longitude
    """
    cell_lats = np.array([])
    cell_lons = np.array([])
    for lat in latCell:
        cell_lats = np.append(cell_lats, lat * (180 / np.pi)) 
    for lon in lonCell:
        cell_lons = np.append(cell_lons, lon * (180 / np.pi))  

    return cell_lats, cell_lons
