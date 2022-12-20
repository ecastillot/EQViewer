import numpy as np
import pandas as pd
import os
import glob
import math
import geopandas as gpd
import numpy as np
import matplotlib.pyplot as plt
import string
from obspy.geodetics.base import gps2dist_azimuth

def truncate(f, n):
    return math.floor(f * 10 ** n) / 10 ** n

def inside_the_polygon(p,pol_points):
    """
    Parameters:
    -----------
    p: tuple
        Point of the event. (lon,lat)
    pol_points: list of tuples
        Each tuple indicates one polygon point (lon,lat).
    Returns: 
    --------
    True inside 
    """
    V = pol_points

    cn = 0  
    V = tuple(V[:])+(V[0],)
    for i in range(len(V)-1): 
        if ((V[i][1] <= p[1] and V[i+1][1] > p[1])   
            or (V[i][1] > p[1] and V[i+1][1] <= p[1])): 
            vt = (p[1] - V[i][1]) / float(V[i+1][1] - V[i][1])
            if p[0] < V[i][0] + vt * (V[i+1][0] - V[i][0]): 
                cn += 1  
    condition= cn % 2  
    
    if condition== 1:   
        return True
    else:
        return False

def get_xy_profile_coords_from_r_phi(r,phi,origin,d2km=114):


    if (phi < 0) or (phi>180):
        raise Exception("phi is azimuth. It can't be minor than 0 or greater than 90")

    x0,y0 = origin
    phi = phi*np.pi/180

    x1 = r*np.sin(phi)/1e3/d2km +x0
    y1 = r*np.cos(phi)/1e3/d2km +y0
    x2 = -r*np.sin(phi)/1e3/d2km +x0
    y2 = -r*np.cos(phi)/1e3/d2km +y0

    return (x1,y1),(x2,y2)

def get_several_rotated_profiles_points(r,delta_phi,origin,
                                phi_0 =0,phi_1=180 ):
    points = {}
    for phi in range(phi_0,phi_1,delta_phi):
        p1,p2= get_xy_profile_coords_from_r_phi(r,phi,origin)
        points[round(phi,2)] = (p1,p2)
    return points    

def get_several_rotated_profiles(r,delta_phi,origin,
                                phi_0 =0,phi_1=180,
                                profiles = {"name":("A","A'"),      
                                            "coords":None,
                                            "width":(-0.3,0.3),
                                            "colorline":"magenta",
                                            "color": None,
                                            "size":None,
                                            "style" :None,
                                            "pen":None,
                                            "cmap":True,
                                            "cbar_profile_args" :{"cmap":'roma',
                                                                    "color_target":"depth",
                                                                    "label":"Depth(m)",
                                                                    "overrule_bg":True,
                                                                    "reverse":False,
                                                                    "series":[0, 3e3] }
                                            }
                                
                                 ):

    points = get_several_rotated_profiles_points(r,delta_phi,origin,
                                         phi_0,phi_1)

    new_profiles = []
    abc = list(string.ascii_lowercase)
    for i,(angle,point) in enumerate(points.items()):
        profile = profiles.copy()
        profile["name"] = (abc[i],abc[i]+"'")
        profile["coords"] = point
        new_profiles.append(profile)

    return new_profiles

def make_grid_profile(min_lon,max_lon,n_lon,
                    min_lat,max_lat,n_lat,
                    r,delta_phi,
                    profiles = {"name":("A","A'"),      
                    "coords":None,
                    "width":(-0.3,0.3),
                    "colorline":"magenta",
                    "color": None,
                    "size":None,
                    "style" :None,
                    "pen":None,
                    "cmap":True,
                    "cbar_profile_args" :{"cmap":'roma',
                                            "color_target":"depth",
                                            "label":"Depth(m)",
                                            "overrule_bg":True,
                                            "reverse":False,
                                            "series":[0, 3e3] }
                    }
                                            ):
    lons = np.arange(min_lon,max_lon,n_lon)                
    lats = np.arange(min_lat,max_lat,n_lat)  

    print(lons)
    print(lats)

    all_profiles = []
    for lon in lons:
        for lat in lats:
            origin = (lon,lat)
            profs = get_several_rotated_profiles(r,delta_phi,origin,
                                                    profiles=profiles)
            all_profiles.append(profs)
    return all_profiles

def get_d_az(p1,p2):
    """
    Parameters:
    -----------
    p1: tuple
        lon, lat
    p2: tuple
        lon,lat
    
    Returns:
    --------
    Get distance and azimuth 
    """

    d2 = (p2[1]-p1[1])**2 + (p2[0]-p1[0])**2
    d = np.sqrt(d2)

    if (p2[0]-p1[0])==0:
        theta = 90
    else:
        theta = np.arctan((p2[1]-p1[1])/(p2[0]-p1[0])) * 180/np.pi

    return d,theta

def get_t_points(p1,p2,d,d2km=114):
    """
    Parameters:
    -----------
    p1 : tuple
        lon,lat 
    p2 : tuple
        lon, lat
    d : float
        distance in km of the transect
    d2km: float
        conversor degrees to kilometers
    """

    _,theta = get_d_az(p1,p2)
    alpha = 90 - theta

    x = d*np.cos(alpha * np.pi/180)
    y = d*np.sin(alpha * np.pi/180)
    x= x/d2km
    y= y/d2km

    pm = ((p2[0]+p1[0])/2 , (p2[1]+p1[1])/2)

    tp1 = (pm[0]-x ,pm[1]+y)
    tp2 = (pm[0]+x,pm[1]-y)

    return tp1, tp2

def get_centers(N,p1,p2):
    
    """
    Paramters:
    ----------
    N: int
        Number of divisions (transects)
    p1: tuple
        (lon, lat) order.
    p2: tuple
        (lon, lat) order.
    Return:
    -------
    center: np.array
        arra of points that indicate the center
    azi: float 
        azimuth between p1 and p2
    """
    d2 = (p2[1]-p1[1])**2 + (p2[0]-p1[0])**2
    d = np.sqrt(d2)
    if (p2[0]-p1[0])==0:
        theta = 90
    else:
        theta = np.arctan((p2[1]-p1[1])/(p2[0]-p1[0])) * 180/np.pi

    d = np.linspace(0,d,num=N,endpoint=True)
    x = p1[0] + d*np.cos(theta*np.pi/180)
    y = p1[1] + d*np.sin(theta*np.pi/180)
    center = np.array(list(zip(x,y)))

    azi = 90 - theta

    return center, azi

def get_line_in_map(center,distance,azimuth,d2k=114,save=None):
    """
    Parameters:
    -----------
    center: tuple
        (lon, lat) order.  
    distance: tuple
        (left_distance,rigth_distance) order in km.
    azimuth: float
        degrees
    d2k: float
        factor of conversion kilometers to degree
    save: str
        Path
    """
    cx,cy = center 
    dl,dr = distance
    azi = azimuth * np.pi/180

    xr = cx + np.sin(azi)*dr/d2k
    yr = cy + np.cos(azi)*dr/d2k
    xl = cx - np.sin(azi)*dl/d2k
    yl = cy - np.cos(azi)*dl/d2k

    alpha =  np.arctan((yr-yl)/(xr-xl) )* 180/np.pi
    # alpha = 90 -alpha
    # print(xr,yr,xl,yl,alpha)

    if save != None:
        df = pd.DataFrame({'lon': [xl, xr], 'lat': [yl, yr]})
        df.to_csv(save,sep = " ", header=False,index=False)
    
    return (xr,yr,xl,yl,alpha)

def points_parallel_to_line(line,d,upper_line=True,d2k=114):
    x0,y0 = line[0]   
    x1,y1 = line[1]   

    delta_x = (x1-x0)
    delta_y = (y1-y0)

    l = np.sqrt(delta_x**2 + delta_y**2)

    d = d/d2k #distance to degrees
    dx = delta_y*d/l
    dy = -delta_x*d/l


    if upper_line:
        x0,y0 = x0 - dx, y0-dy 
        x1,y1 = x1 - dx, y1-dy 
    else:
        x0,y0 = x0 + dx, y0+dy 
        x1,y1 = x1 + dx, y1+dy 
    return ((x0,y0),(x1,y1))

def get_cmap(n, name='hsv'):
    '''Returns a function that maps each index in 0, 1, ..., n-1 to a distinct 
    RGB color; the keyword argument name must be a standard mpl colormap name.'''
    return plt.cm.get_cmap(name, n)

def get_df(csv, between=None,
            columns=["longitude","latitude","depth"]):

    events = pd.read_csv(csv)

    events["time"] = pd.to_datetime(events["time"])
    if between != None:
        events = events[(events["time"]>between[0]) & (events["time"]<between[1]) ]

    events = events[columns]
    return events

def get_glob_gdf(folder,find,fmt="shp"):
    """
    You could use glob.glob wildcars in the 'find' parameter.
    """

    filepath = os.path.join(folder,find)
    data = []
    for path in glob.glob(f"{filepath}*.{fmt}"):
        base = os.path.basename(path)
        label = os.path.splitext(base)[0]
        gdf  = gpd.read_file(path)
        gdf["shapename"] = label

        data.append(gdf)
    data = gpd.GeoDataFrame( pd.concat( data, ignore_index=True) )
    return data

def get_gdf(gdf,attrs=None,
            proj="EPSG:4326"):
    """
    returns the geopandas dataframe with specific attributes
    """


    if proj != "EPSG:4326":
        gdf = gdf.to_crs(proj)
    else: 
        gdf = gdf.to_crs("EPSG:4326")


    if attrs != None: 
        data = []   
        for key,values in attrs.items():
            for value in values:
                df = gdf[gdf[key] == value]
                data.append(df)
    else: 
        data = [gdf]
    
    data = gpd.GeoDataFrame( pd.concat( data, ignore_index=True) )
    return data

def get_cross_profiles(coords,n_cross,theta=90,
        number_of_profile=0,
        profiles = {"name":("A","A'"),      
        "coords":None,
        "width":(-0.3,0.3),
        "colorline":"magenta",
        "color": None,
        "size":None,
        "style" :None,
        "pen":None,
        "cmap":True,
        "cbar_profile_args" :{"cmap":'roma',
                                "color_target":"depth",
                                "label":"Depth(m)",
                                "overrule_bg":True,
                                "reverse":False,
                                "series":[0, 3e3] }
        }
    ):
    """
    coords: list of tuples
        i.e. ((-73.685857,3.870372),(-73.674690,3.876694))
    """
    profiles = profiles.to_dict()
    profiles["coords"] = coords
    p1,p2 = profiles["coords"]
    l = profiles["width"]
    l = [abs(i) for i in l]
    centers,azi = get_centers(n_cross,p1,p2) ##lineas perpendiculares
    r,a,ba = gps2dist_azimuth(centers[0][1],centers[0][0],
                        centers[1][1],centers[1][0])
    w_norm = r/2/1e3
    width = (-w_norm,w_norm)
    abc = list(string.ascii_lowercase)


    cross_profiles = []
    for i,center in enumerate(centers,number_of_profile):
            xr,yr,xl,yl,alpha = get_line_in_map(center,l,azi+theta,
                                                    d2k=114,save=None)
            name = (abc[i],abc[i]+"'")
            point_1 = (xl,yl)
            point_2 = (xr,yr)

            cross_profile = profiles.copy()
            cross_profile["name"] = name
            cross_profile["coords"] = (point_1,point_2)
            cross_profile["width"] = width


            cross_profile = Profile(**cross_profile)

            cross_profiles.append(cross_profile)

    profiles= cross_profiles
    return profiles

def transform_to_fmt_catalog(csv,
    columns={"Origin time":"origin_time",
            "Latitude (deg)":"latitude",
            "Longitude (deg)":"longitude",
            "Depth (m)":"depth",
            "Mag.":"magnitude",
            "Mag. (Mw)":"magnitude",
            "Strike n1 (deg)":"strike_n1",
            "Dip n1 (deg)":"dip_n1",
            "Dip n2 (deg)":"dip_n2",
            "Rake n1 (deg)":"rake_n1",
            "Rake n2 (deg)":"rake_n2",
            "Strike n1 (deg)":"strike_n1",
            "Strike n2 (deg)":"strike_n2",
            "Dip n2 (deg)":"dip_n2",
            
            },
    origin_time_fmt = None,
    sep=","):

    df = pd.read_csv(csv,sep=sep)
    df = df.rename(columns=columns)

    if origin_time_fmt == None:
        df["origin_time"] = pd.to_datetime(df['origin_time']).dt.tz_localize(None)
    else:
        df["origin_time"] = pd.to_datetime(df['origin_time'],
                                            format=origin_time_fmt).dt.tz_localize(None)


    df = df.sort_values(["origin_time"],ignore_index=True)
    return df
if __name__ == "__main__":
    # x = get_line_in_map((-76.45,1.542),(5,5),35,d2k=114,save=None)
    # print(x)

    # wells = "/home/emmanuel/G-Ecopetrol/ecastillo/Avances/2022/Asociacion_castilla/qgis/Pozos"
    # get_glob_gdf(wells,"CLIA*",fmt="shp")


    # r = 2e3
    # phi = 45
    # origin = (-73.667240,3.820862)
    # # p1,p2 = get_xy_profile_coords_from_r_phi(r,phi,origin)
    # # print(p1,p2)

    # points = get_several_rotated_profiles(r,phi,origin)
    # print(points)

    profiles = make_grid_profile(-73.682718,-73.645,0.005,
                    3.801395,3.830791,0.005,
                    1e3,30)
    print(profiles)

    