import re

station_list_search='/nfs/a90/eepdw/Data/Observations/Radiosonde_downloaded_from_NOAA_GUAN/igra-stations.txt'

def StationListParse(station_list_search):

    '''
    '''

    station_metadata=[]
    f = open(station_list_search,'r')
    for line in f:
        line = line.strip()
        line=re.sub(r'([A-Z])\s([A-Z])', r'\1_\2',line)
        line=re.sub(r'([A-Z])\s\s([A-Z])', r'\1_\2',line)
        station_metadata.append(line.split())
    f.close()
    return station_metadata

station_metadata = StationListParse(station_list_search)

def StationInfoSearch(stat):

    '''
    '''
    for line in station_metadata:
         if "%s" % stat in line: 
             st = line[2].lower().title().replace('_',' ')
             lo = float(line[3])
             la = float(line[4])
             st_height = float(line[5])
    return st,la,lo, st_height

def CalculateDistanceFromFirstStation(stat, first_station_lon, first_station_lat, station_lat, station_lon):

    fslat_rad = radians(first_station_lat)
    fslon_rad = radians(first_station_lon)
    lat_rad = radians(station_lat)
    lon_rad = radians(station_lon)

    #Haversine Formula
    
    a = sin((lat_rad-fslat_rad)/2)**2 + cos(lat_rad) * cos(fslat_rad) * sin((lon_rad-fslon_rad)/2)**2
    c = 2 * atan2(sqrt(a), sqrt(1-a))
    d = 6371 * c

    return d
