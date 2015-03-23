import datetime
from dateutil.relativedelta import relativedelta

from calendar import timegm

#date_min=datetime.datetime(1960,5,1,0,0,0)
#date_max=datetime.datetime(2014,10,1,0,0,0)

def YearlyGrid(date_min, date_max):

    delta = relativedelta(years=+1)

    d = date_min

    grid=[]
    while (d <= date_max): 
        grid.append(d.timetuple().tm_year)
        d += delta
    
    return grid\

def DiurnalHourlyGrid(date_min, date_max):

    delta = relativedelta(hours=+1)

    d = date_min

    grid=[]
    while ((d <= date_max) and (d.timetuple().tm_hour not in grid)): 
        grid.append(d.timetuple().tm_hour)
        d += delta
    
    return grid

def ClimatolWeeklyGrid(date_min, date_max):

    delta = relativedelta(weeks=+1)

    d = date_min

    grid=[]
    while ((d <= date_max) and (d.timetuple() not in grid)): 
        grid.append(d.timetuple().tm_yday)
        d += delta
    
    return grid

def ClimatolMonthlyGrid(date_min, date_max):

    delta = relativedelta(months=+1)

    d = date_min

    grid=[]
    while ((d <= date_max) and (d.timetuple() not in grid)): 
        grid.append(d.timetuple().tm_yday)
        d += delta
    
    return grid
def UnixEpochWeeklyGrid(date_min, date_max):

    delta = relativedelta(weeks=+1)

    d = date_min

    grid=[]
    while ((d <= date_max) and (d.timetuple() not in grid)): 
        grid.append(timegm(d.timetuple()))
        d += delta
    
    return grid

def Climatol3WeekGrid(date_min, date_max):

    delta = relativedelta(weeks=+3)

    d = date_min

    grid=[]
    while ((d <= date_max) and (d.timetuple() not in grid)): 
        grid.append(d.timetuple().tm_yday)
        d += delta
    
    return grid

#def BinData()
