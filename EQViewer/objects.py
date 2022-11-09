# /**
#  * @author [Emmanuel Castillo]
#  * @email [ecastillot@unal.edu.co]
#  * @create date 2022-11-03 12:07:04
#  * @modify date 2022-11-03 12:07:04
#  * @desc [description]
#  */
from operator import add
import copy
import types
import pygmt
import pandas as pd
import geopandas as gpd
import datetime as dt
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
from obspy.geodetics.base import gps2dist_azimuth
from . import utils as ut
pygmt.config(FORMAT_GEO_MAP="ddd.xx")

def get_cmap(n, name='hsv'):
    '''Returns a function that maps each index in 0, 1, ..., n-1 to a distinct 
    RGB color; the keyword argument name must be a standard mpl colormap name.'''
    return plt.cm.get_cmap(name, n)

def args_cleaner(args,rm_args=[]):
    info = args
    for key in list(info.keys()):
        if key in rm_args:
            info.pop(key,None)

    return info

def update_kwargs(args1,args2):
    """
    keys in args1 will update to arg2
    """
    for key,value in args1.items():
        if key in list(args2.keys()):
            args2[key] = value
    return args2

# class BaseText():

class BasePlot():
    def __init__(self,
            size=None,
            style="c0.2c",
            cmap=False,
            color="lightblue",
            label="data",
            transparency=0,
            pen=None,
            **plot_kwargs) -> None:
        """
        Parameters:
        -----------
        size: None or lambda function
            Equation for the size. 
            lambda x: 0.1 * np.sqrt(1.5 ** (x.magnitude*1.5))
            in this case size is assigned to the magnitude.
            
            If size is defined as lambda function,
            you must use 'style':"cc"
        style: str
            style of you data. 
            First letter assing the symbol.
            Second letter assign the measure distance
            If there is a number between the letter, means the
            size for every point.
            
            For instance, use c0.2c circle,0.2 centimeters for all data
            use cc, for circle and centimeters (this time size must be specified)
        cmap: bool
            Use Colorbar (the specifications of the colorbar always are located in '.plot' functions). 
        color: str or None
            Color from pygmt color gallery. 
            It is not considered when cmap=True
        transparency: float
            transparency of your plots
        pen : str or None
            color and size of the symbol border
        plot_kwargs: other pygmt.plot arguments
        """
        self.size = size
        self.color = color
        self.label = label
        self.style = style
        self.cmap = cmap
        self.transparency = transparency
        self.pen = pen
        self.plot_kwargs = plot_kwargs

    def get_info2pygmt(self,data=pd.DataFrame):
        rm_args = ["plot_kwargs"]
        args = self.__dict__.copy()
        
        if type(self.size) is types.LambdaType:
            if not data.empty:
                args["size"] = data.apply(args["size"],axis=1)
            else:
                raise Exception("No data to apply size")

        kwargs = update_kwargs(args,self.plot_kwargs)
        args.update(kwargs)
        args = args_cleaner(args,rm_args)
        return args

class CPT():
   def __init__(self,color_target,label,**makecpt_kwargs):
        """
        Parameters:
        -----------
        color_target: "str"
            Target name to apply the colorbar.
        label: "str"
            Label to show in the colorbar legend
        makecpt_kwargs:
            Args from pygmt.makecpt
        """
        self.color_target = color_target
        self.label = label
        self.makecpt_kwargs = makecpt_kwargs

class Catalog():
    def __init__(self,
            data,
            baseplot = BasePlot()) -> None:

        """
        Parameters:
        -----------
        data: pd.DataFrame 
            Dataframe with the next mandatory columns:
            'origin_time','latitude','longitude','depth','magnitude'
        baseplot: BasePlot
            Control plot args
        """
        self.columns = ['origin_time','latitude','longitude','depth','magnitude']
        check =  all(item in data.columns.to_list() for item in self.columns)
        if not check:
            raise Exception("There is not the mandatory columns for the data in Catalog object."\
                            +"->'origin_time','latitude','longitude','depth','magnitude'")

        data = data.drop_duplicates(subset=self.columns,ignore_index=True)
        pd.to_datetime(data.loc[:,"origin_time"]).dt.tz_localize(None)
        # self.data = data[self.columns]
        self.data = data
        if self.empty:
            raise Exception("No data in the catalog")

        self.baseplot = baseplot

    @property
    def empty(self):
        return self.data.empty

    def __len__(self):
        return len(self.data)

    def __str__(self,extended=False) -> str:
        if extended:
            timefmt = "%Y%m%dT%H:%M:%S"
            start=  self.data.origin_time.min()
            end = self.data.origin_time.max()
            region = list(map(lambda x: round(x,2),self.get_region()))
            msg = f"Catalog | {self.__len__()} events "\
                    +f"\n\tperiod: [{start.strftime(timefmt)} - {end.strftime(timefmt)}]"\
                    +f"\n\tdepth : {[round(self.data.depth.min(),2),round(self.data.depth.max(),2)]}"\
                    +f"\n\tmagnitude : {[round(self.data.magnitude.min(),2),round(self.data.magnitude.max(),2)]}"\
                    +f"\n\tregion: {region}"
        else:
            msg = f"Catalog | {self.__len__()} events "

        return msg

    def append(self, data):
        """
        append data
        """
        if isinstance(data, pd.DataFrame):

            check =  all(item in data.columns.to_list() for item in self.columns)
            if not check:
                raise Exception("There is not the mandatory columns for the data in Catalog object."\
                                +"->'origin_time','latitude','longitude','depth','magnitude'")

            data = data.drop_duplicates(subset=self.columns,ignore_index=True)
            pd.to_datetime(data.loc[:,"origin_time"]).dt.tz_localize(None)

            self.data = pd.concat([self.data,data])
        else:
            msg = 'Append only supports a single Dataframe object as an argument.'
            raise TypeError(msg)
        return self

    def remove_data(self, rowval):
        """
        remove rows to the data.

        Parameters:
        -----------
        rowval : dict
            key: 
                column name
            value: list
                values specified to remove
        """
        if not isinstance(rowval,dict):
            raise Exception("rowval must be a dictionary")
        
        mask = self.data.isin(rowval)
        mask = mask.any(axis='columns')
        self.data = self.data[~mask]
        return self
    
    def select_data(self, rowval):
        """
        select rows to the data.

        Parameters:
        -----------
        rowval : dict
            key: 
                column name
            value: list
                values specified to select
        """
        if not isinstance(rowval,dict):
            raise Exception("rowval must be a dictionary")
        mask = self.data.isin(rowval)
        mask = mask.any(axis='columns')
        self.data = self.data[mask]
        return self

    def copy(self):
        """Deep copy of the class"""
        return copy.deepcopy(self)

    def sort_values(self,**args):
        """
        Sort values. Take in mind that it could affect the order of the events plotted
        args: The parameters are the pd.DataFrame.sort_values parameters
        """
        self.data = self.data.sort_values(**args)
        return self

    def filter_datetime(self,starttime=None,endtime=None):
        """
        Filter the period of the catalog.

        Parameters:
        -----------
        starttime: datetime.datetime
            start time
        endtime: datetime.datetime
            end time
        
        """
        if starttime != None and \
            not isinstance(starttime,dt.datetime):
            raise Exception("starttime must be a datetime object")

        if endtime != None and \
            not isinstance(endtime,dt.datetime):
            raise Exception("starttime must be a datetime object")

        if isinstance(starttime,dt.datetime) \
            and isinstance(endtime,dt.datetime):
            if endtime < starttime:
                raise Exception("endtime must be greater than starttime")

        if isinstance(starttime,dt.datetime):
            self.data = self.data[self.data["origin_time"]>=starttime]
        if isinstance(endtime,dt.datetime):
            self.data = self.data[self.data["origin_time"]<=endtime]
        return self

    def filter_region(self,polygon):
        """
        Filter the region of the catalog.

        Parameters:
        -----------
        polygon: list of tuples
            Each tuple is consider a point (lon,lat).
            The first point must be equal to the last point in the polygon.
        
        """
        if polygon[0] != polygon[-1]:
            raise Exception("The first point must be equal to the last point in the polygon.")

        is_in_polygon = lambda x: ut.inside_the_polygon((x.longitude,x.latitude),polygon)
        mask = self.data[["longitude","latitude"]].apply(is_in_polygon,axis=1)
        self.data = self.data[mask]
        return self

    def get_region(self,padding=[]):
        """
        It gets the region according to the limits in the catalog
        
        Parameters:
        -----------
        padding: 4D-list or float or int
            list: Padding on each side of the region [lonw,lonw,lats,latn] in degrees.
            float or int: padding amount on each side of the region from 0 to 1,
                        where 1 is considered the distance on each side of the region.
        """
        lonw,lone = self.data.longitude.min(),self.data.longitude.max()
        lats,latn = self.data.latitude.min(),self.data.latitude.max()
        
        region = [lonw, lone, lats, latn]
        
        if isinstance(padding,list):
            if padding:
                if len(padding) != 4:
                    raise Exception("Padding parameter must be 4D")
                else:
                    region = list( map(add, region, padding) )
        elif isinstance(padding,float) or isinstance(padding,int):
            lon_distance = abs(region[1]-region[0])
            lat_distance = abs(region[3]-region[2])
            adding4lon = lon_distance*padding
            adding4lat = lat_distance*padding
            padding = [-adding4lon, adding4lon, -adding4lat, adding4lat]
            region = list( map(add, region, padding) )

        return region

    def project(self,startpoint,endpoint,
                width,verbose=True):
        """
        Project data onto a line

        Parameters:
        -----------
        startpoint: tuple
            (lon,lat)
        endpoint: tuple
            (lon,lat)
        width: tuple
            (w_left,w_right)
        
        """
        data = self.data
        data = data[["longitude","latitude","depth",
                    "origin_time","magnitude"]]
        projection = pygmt.project(
                            data=data,
                            unit=True,
                            center=startpoint,
                            endpoint=endpoint,
                            convention="pz",
                            width=width,
                            verbose=verbose
                                )
        projection = projection.rename(columns={0:"distance",
                                1:"depth",
                                2:"origin_time",
                                3:"magnitude"})
        return projection

    def plot(self,fig=None,
            cpt=None,
            show_cpt=True):
        """
        Plot the catalog.

        Parameters:
        -----------
        fig: None or pygmt.Figure
            Basemap figure
        cpt: None or CPT
            color palette table applied to the catalog
        show_cpt: bool
            Show the color palette table
        """
        data = self.data

        if fig == None:
            fig = pygmt.Figure() 
            fig.basemap(region=self.get_region(padding=0.1),
                        projection="M12c", 
                        frame=["afg","WNse"])
        
        info2pygmt = self.baseplot.get_info2pygmt(data)
        if self.baseplot.cmap:
            if cpt == None:
                zmin = data.depth.min()
                zmax = data.depth.max()
                cpt = CPT(color_target="depth",
                            label="depth",
                            cmap="rainbow",
                            series=[zmin,zmax],
                            reverse=True,
                            overrule_bg=True)

            info2pygmt["color"] = data[cpt.color_target]
            pygmt.makecpt(**cpt.makecpt_kwargs)
            
            if show_cpt:
                fig.colorbar(frame=f'af+l"{cpt.label}"',
                        position="JBC+e")
        fig.plot(
            x=data.longitude,
            y=data.latitude,
            **info2pygmt
        )

        return fig

    def matplot(self,color_target="depth",
            s=8,cpt="viridis",show_cpt=True,
            ax=None):
        """
        Quickly matplotlib figure

        Parameters:
        color_target: str
            target to apply cbar
        s: float
            marker size
        cpt: str
            Name of the colorbar
        show_cpt: bool
            Show color palette table.
        ax: axis
            existent axis
        """

        if ax == None:
            fig = plt.figure(figsize=(10, 10))
            ax = fig.add_subplot(111)
        ax.set_aspect("equal")


        if color_target == "origin_time":
            cb = ax.scatter(self.data.longitude, self.data.latitude,
                    c=mdates.date2num(self.data[color_target]), s=s, cmap=cpt)
            
            if show_cpt:
                cpt = fig.colorbar(cb)
                cpt.ax.set_ylim(cpt.ax.get_ylim()[::-1])
                cpt.set_label(f"{color_target}")
                
                loc = mdates.AutoDateLocator()
                cpt.ax.yaxis.set_major_locator(loc)
                cpt.ax.yaxis.set_major_formatter(mdates.ConciseDateFormatter(loc))
        else:
            cb = ax.scatter(self.data.longitude, self.data.latitude,
                    c=self.data[color_target], s=s, cmap=cpt)
            if show_cpt:
                cpt = fig.colorbar(cb)
                cpt.ax.set_ylim(cpt.ax.get_ylim()[::-1])
                cpt.set_label(f"{color_target}")

        ax.set_xlabel("Longitude [°]")
        ax.set_ylabel("Latitude [°]")
        return ax

class Catalogs():
    def __init__(self,catalogs=[],cpt=None,show_cpt=True):
        """
        Parameters:
        -----------
        catalogs: list
            list of Catalog objects
        cpt: None or CPT
            color palette table applied to the catalog
        """
        self.catalogs = catalogs
        self.cpt = cpt
        self.show_cpt = show_cpt

    def __iter__(self):
        return list(self.catalogs).__iter__()

    def __nonzero__(self):
        return bool(len(self.catalogs))

    def __len__(self):
        return len(self.catalogs)
    
    def __str__(self,extended=False) -> str:
        msg = f"Catalogs ({self.__len__()} catalogs)\n"
        msg += "-"*len(msg) 

        submsgs = []
        for i,catalog in enumerate(self.__iter__(),1):
            submsg = f"{i}. "+catalog.__str__(extended=extended)
            submsgs.append(submsg)
                
        if len(self.catalogs)<=20 or extended is True:
            submsgs = "\n".join(submsgs)
        else:
            three_first_submsgs = submsgs[0:3]
            last_two_subsgs = submsgs[-2:]
            len_others = len(self.catalogs) -len(three_first_submsgs) - len(last_two_subsgs)
            submsgs = "\n".join(three_first_submsgs+\
                                [f"...{len_others} other catalogs..."]+\
                                last_two_subsgs)

        return msg+ "\n" +submsgs
        
    def __setitem__(self, index, trace):
        self.catalogs.__setitem__(index, trace)

    def __getitem__(self, index):
        if isinstance(index, slice):
            return self.__class__(catalogs=self.catalogs.__getitem__(index))
        else:
            return self.catalogs.__getitem__(index)

    def __delitem__(self, index):
        return self.catalogs.__delitem__(index)

    def __getslice__(self, i, j, k=1):
        return self.__class__(catalogs=self.catalogs[max(0, i):max(0, j):k])

    def append(self, catalog):
        """
        append a catalog
        """
        if isinstance(catalog, Catalog):
            self.catalogs.append(catalog)
        else:
            msg = 'Append only supports a single Catalog object as an argument.'
            raise TypeError(msg)
        return self

    def remove_data(self,rowval):
        """
        remove rows of each data.

        Parameters:
        -----------
        rowval : dict
            key:  
                column name
            value: 
                One or more values specified to remove
        """
        catalogs = []
        for catalog in self.catalogs:
            catalogs.append(catalog.remove_data(rowval))
        self.catalogs = catalogs
        return self

    def select_data(self,rowval):
        """
        select rows of each data.

        Parameters:
        -----------
        rowval : dict
            key:  
                column name
            value: 
                One or more values specified to select
        """
        catalogs = []
        for catalog in self.catalogs:
            catalogs.append(catalog.select_data(rowval))
        self.catalogs = catalogs
        return self

    def copy(self):
        """Deep copy of the class"""
        return copy.deepcopy(self)

    def project(self,startpoint,endpoint,
                width,verbose=True):
        projections = []
        for catalog in self.catalogs:
            projection = catalog.project(startpoint,endpoint,
                                            width,verbose)
            projections.append(projection)
        return projections

    def sort_values(self,**args):
        """
        Sort values. Take in mind that it could affect the order of the events plotted
        args: The parameters are the pd.DataFrame.sort_values parameters
        """
        catalogs = []
        for catalog in self.catalogs:
            catalogs.append(catalog.sort_values(**args))
        self.catalogs = catalogs
        return self
    
    def filter_datetime(self,starttime=None,endtime=None):
        """
        Filter the period of the catalog.

        Parameters:
        -----------
        starttime: datetime.datetime
            start time
        endtime: datetime.datetime
            end time
        
        """
        catalogs = []
        for catalog in self.catalogs:
            catalogs.append(catalog.filter_datetime(starttime,endtime))
        self.catalogs = catalogs
        return self

    def filter_region(self,polygon):
        """
        Filter the region of the catalog.

        Parameters:
        -----------
        polygon: list of tuples
            Each tuple is consider a point (lon,lat).
            The first point must be equal to the last point in the polygon.
        
        """
        catalogs = []
        for catalog in self.catalogs:
            catalogs.append(catalog.filter_region(polygon))
        self.catalogs = catalogs
        return self

    def get_region(self,padding=[]):
        """
        It gets the region according to the limits of all events as a whole

        Parameters:
        -----------
        padding: 4D-list or float or int
            list: Padding on each side of the region [lonw,lonw,lats,latn] in degrees.
            float or int: padding amount on each side of the region from 0 to 1,
                        where 1 is considered the distance on each side of the region.
        """
        lons,lats = [],[]
        for catalog in self.catalogs:
            region = catalog.get_region()
            lons.append(region[0:2])
            lats.append(region[2:])
        lons = [ x for sublist in lons for x in sublist]
        lats = [ x for sublist in lats for x in sublist]
        region = [min(lons),max(lons),min(lats),max(lats)]

        if isinstance(padding,list):
            if padding:
                if len(padding) != 4:
                    raise Exception("Padding parameter must be 4D")
                else:
                    region = list( map(add, region, padding) )
        elif isinstance(padding,float) or isinstance(padding,int):
            lon_distance = abs(region[1]-region[0])
            lat_distance = abs(region[3]-region[2])
            adding4lon = lon_distance*padding
            adding4lat = lat_distance*padding
            padding = [-adding4lon, adding4lon, -adding4lat, adding4lat]
            region = list( map(add, region, padding) )

        return region

    def plot(self,fig=None):

        """
        Plot the catalog.

        Parameters:
        -----------
        fig: None or pygmt.Figure
            Basemap figure
        """

        if fig == None:
            fig = pygmt.Figure() 
            fig.basemap(region=self.get_region(padding=0.1),
                        projection="M12c", 
                        frame=["afg","WNse"])
        if self.cpt == None:
            data = []
            for catalog in self.catalogs:
                if catalog.baseplot.cmap:
                    data.append(catalog.data)
            data = pd.concat(data)
            zmin = data.depth.min()
            zmax = data.depth.max()
            self.cpt = CPT(color_target="depth",
                        label="depth",
                        cmap="rainbow",
                        series=[zmin,zmax],
                        reverse=True,
                        overrule_bg=True)

        show_catalog_cpt = []
        for catalog in self.catalogs:
            if catalog.baseplot.cmap:
                catalog.plot(fig=fig,cpt=self.cpt,show_cpt=False)
                _show_cpt = True
            else:
                catalog.plot(fig=fig,cpt=None,show_cpt=False)
                _show_cpt = False
            show_catalog_cpt.append(_show_cpt)

        if any(show_catalog_cpt):
            if self.show_cpt:
                fig.colorbar(frame=f'af+l"{self.cpt.label}"',
                        position="JBC+e")

        return fig
    
class Station():
    def __init__(self,data,
                name_in_map=True,
                color="black",
                label="stations",
                transparency = 0,
                style="i0.3c",
                pen="black",
                plot_args={},
                text_args={"font":"10p,Helvetica,black",
                            "fill":None,
                            "offset":"-0.05c/0.15c"}) -> None:

        """
        data: DataFrame
            Dataframe with the next mandatory columns:
            'station','latitude','longitude'
        name_in_map : bool
            Show the name of the station. 
        color: str or None
            Color from pygmt color gallery. 
        label: str
            Label name of the stations data
        transparency: float
            transparency of your plots
        style: str
            style of you data. 
            First letter assing the symbol.
            Second letter assign the measure distance
            If there is a number between the letter, means the
            size for every point.
            
            For instance, use c0.2c circle,0.2 centimeters for all data
        pen : str
            color and size of the symbol border
        plot_args: dict
            Arguments of the pygmt.plot
        plot_text: dict
            Arguments of the pygmt.text
        """
        self.columns = ['station','latitude','longitude']
        check =  all(item in data.columns.to_list() for item in self.columns)
        if not check:
            raise Exception("There is not the mandatory columns for the data in Station object."\
                            +"->'station','latitude','longitude'")
        # self.data = data[columns]
        self.data = data
        self.name_in_map = name_in_map
        self.color = color
        self.label = label
        self.style = style
        self.pen = pen
        self.transparency = transparency
        self.plot_args = plot_args
        self.text_args = text_args


    @property
    def empty(self):
        return self.data.empty
    
    @property
    def info2pygmt(self):
        rm_args = ["data","name_in_map","text_args",
                    "plot_args","columns"]
        args = self.__dict__.copy()
        args.update(self.plot_args)
        info_dict = args_cleaner(args,rm_args)
        return info_dict.copy()

    def get_region(self,padding=[]):
        """
        It gets the region according to the limits in the catalog
        Parameters:
        -----------
        padding: 4D-list or float or int
            list: Padding on each side of the region [lonw,lonw,lats,latn] in degrees.
            float or int: padding amount on each side of the region from 0 to 1,
                        where 1 is considered the distance on each side of the region.
        """
        lonw,lone = self.data.longitude.min(),self.data.longitude.max()
        lats,latn = self.data.latitude.min(),self.data.latitude.max()
        region = [lonw, lone, lats, latn]
        if isinstance(padding,list):
            if padding:
                if len(padding) != 4:
                    raise Exception("Padding parameter must be 4D")
                else:
                    region = list( map(add, region, padding) )
        elif isinstance(padding,float) or isinstance(padding,int):
            lon_distance = abs(region[1]-region[0])
            lat_distance = abs(region[3]-region[2])
            adding4lon = lon_distance*padding
            adding4lat = lat_distance*padding
            padding = [-adding4lon, adding4lon, -adding4lat, adding4lat]
            region = list( map(add, region, padding) )

        return region

    def __len__(self):
        return len(self.data)

    def __str__(self,extended=False) -> str:
        msg = f"Station | {self.__len__()} stations"
        if extended:
            region = list(map(lambda x: round(x,2),self.get_region()))
            msg += f"\n\tregion: {region}"
        else:
            pass
        return msg

    def remove_data(self, rowval):
        """
        remove rows to the data.

        Parameters:
        -----------
        rowval : dict
            key: 
                column name
            value: list
                values specified to remove
        """
        if not isinstance(rowval,dict):
            raise Exception("rowval must be a dictionary")
        
        mask = self.data.isin(rowval)
        mask = mask.any(axis='columns')
        self.data = self.data[~mask]
        return self

    def select_data(self, rowval):
        """
        select rows to the data.

        Parameters:
        -----------
        rowval : dict
            key: 
                column name
            value: list
                values specified to select
        """
        if not isinstance(rowval,dict):
            raise Exception("rowval must be a dictionary")
        mask = self.data.isin(rowval)
        mask = mask.any(axis='columns')
        self.data = self.data[mask]
        return self

    def append(self, data):
        """
        append data
        """
        if isinstance(data, pd.DataFrame):

            check =  all(item in data.columns.to_list() for item in self.columns)
            if not check:
                raise Exception("There is not the mandatory columns for the data in Catalog object."\
                                +"->'origin_time','latitude','longitude','depth','magnitude'")

            data = data.drop_duplicates(subset=self.columns,ignore_index=True)
            self.data = pd.concat([self.data,data])
        else:
            msg = 'Append only supports a single Dataframe object as an argument.'
            raise TypeError(msg)
        return self

    def copy(self):
        """Deep copy of the class"""
        return copy.deepcopy(self)

    def sort_values(self,**args):
        """
        Sort values. Take in mind that it could affect the order of the stations plotted
        args: The parameters are the pd.DataFrame.sort_values parameters
        """
        self.data = self.data.sort_values(**args)
        return self

    def filter_region(self,polygon):
        """
        Filter the region of the data.

        Parameters:
        -----------
        polygon: list of tuples
            Each tuple is consider a point (lon,lat).
            The first point must be equal to the last point in the polygon.
        
        """
        if polygon[0] != polygon[-1]:
            raise Exception("The first point must be equal to the last point in the polygon.")

        is_in_polygon = lambda x: ut.inside_the_polygon((x.longitude,x.latitude),polygon)
        mask = self.data[["longitude","latitude"]].apply(is_in_polygon,axis=1)
        self.data = self.data[mask]
        return self

    def plot(self,fig=None):
        """
        Plot the stations

        Parameters:
        -----------
        fig: None or pygmt.Figure
            Basemap figure
        """
        data = self.data

        if fig == None:
            fig = pygmt.Figure() 
            fig.basemap(region=self.get_region(padding=0.1),
                        projection="M12c", 
                        frame=["afg","WNse"])
        fig.plot(
                x=data["longitude"],
                y=data["latitude"],
                **self.info2pygmt
            )
        if self.name_in_map:
            fig.text(x=data["longitude"], 
                    y=data["latitude"], 
                    text=data["station"],
                    **self.text_args)
        return fig

    def matplot(self,ax=None):
        """
        Quickly matplotlib figure
        """
        if ax == None:
            fig = plt.figure(figsize=(10, 10))
            ax = fig.add_subplot(111)
        ax.set_aspect("equal")
        ax.plot(self.data.longitude, self.data.latitude, "r^", 
                ms=10, mew=1, mec="k")
        for i,row in self.data.iterrows():
            ax.text(row.longitude, row.latitude, row.station, fontsize=12)
        ax.set_xlabel("Longitude [°]")
        ax.set_ylabel("Latitude [°]")
        return ax

class Stations():
    def __init__(self,stations=[]):
        """
        Parameters:
        -----------
        stations: list
            list of Station objects
        """
        self.stations = stations

    def __iter__(self):
        return list(self.stations).__iter__()

    def __nonzero__(self):
        return bool(len(self.stations))

    def __len__(self):
        return len(self.stations)
    
    def __str__(self,extended=False) -> str:
        msg = f"Stations ({self.__len__()} stations)\n"
        msg += "-"*len(msg) 

        submsgs = []
        for i,station in enumerate(self.__iter__(),1):
            submsg = f"{i}. "+station.__str__(extended=extended)
            submsgs.append(submsg)
                
        if len(self.stations)<=20 or extended is True:
            submsgs = "\n".join(submsgs)
        else:
            three_first_submsgs = submsgs[0:3]
            last_two_subsgs = submsgs[-2:]
            len_others = len(self.stations) -len(three_first_submsgs) - len(last_two_subsgs)
            submsgs = "\n".join(three_first_submsgs+\
                                [f"...{len_others} other catalogs..."]+\
                                last_two_subsgs)

        return msg+ "\n" +submsgs

    def __getitem__(self, index):
        if isinstance(index, slice):
            return self.__class__(stations=self.stations.__getitem__(index))
        else:
            return self.stations.__getitem__(index)

    def __delitem__(self, index):
        return self.stations.__delitem__(index)

    def __getslice__(self, i, j, k=1):
        return self.__class__(stations=self.stations[max(0, i):max(0, j):k])

    def get_region(self,padding=[]):
        """
        It gets the region according to the limits of all stations as a whole

        Parameters:
        -----------
        padding: 4D-list or float or int
            list: Padding on each side of the region [lonw,lonw,lats,latn] in degrees.
            float or int: padding amount on each side of the region from 0 to 1,
                        where 1 is considered the distance on each side of the region.
        """
        lons,lats = [],[]
        for station in self.stations:
            region = station.get_region()
            lons.append(region[0:2])
            lats.append(region[2:])
        lons = [ x for sublist in lons for x in sublist]
        lats = [ x for sublist in lats for x in sublist]
        region = [min(lons),max(lons),min(lats),max(lats)]

        if isinstance(padding,list):
            if padding:
                if len(padding) != 4:
                    raise Exception("Padding parameter must be 4D")
                else:
                    region = list( map(add, region, padding) )
        elif isinstance(padding,float) or isinstance(padding,int):
            lon_distance = abs(region[1]-region[0])
            lat_distance = abs(region[3]-region[2])
            adding4lon = lon_distance*padding
            adding4lat = lat_distance*padding
            padding = [-adding4lon, adding4lon, -adding4lat, adding4lat]
            region = list( map(add, region, padding) )

        return region
    
    def append(self, stations):
        """
        append a stations
        """
        if isinstance(stations, Station):
            self.stations.append(stations)
        else:
            msg = 'Append only supports a single Stations object as an argument.'
            raise TypeError(msg)
        return self

    def remove_data(self,rowval):
        """
        remove rows of each data.

        Parameters:
        -----------
        rowval : dict
            key:  
                column name
            value: 
                One or more values specified to remove
        """
        stations = []
        for station in self.stations:
            stations.append(station.remove_data(rowval))
        self.stations = stations
        return self

    def select_data(self,rowval):
        """
        select rows of each data.

        Parameters:
        -----------
        rowval : dict
            key:  
                column name
            value: 
                One or more values specified to select
        """
        stations = []
        for station in self.stations:
            stations.append(station.select_data(rowval))
        self.stations = stations
        return self

    def copy(self):
        """Deep copy of the class"""
        return copy.deepcopy(self)    

    def sort_values(self,**args):
        """
        Sort values. Take in mind that it could affect the order of the stations plotted
        args: The parameters are the pd.DataFrame.sort_values parameters
        """
        stations = []
        for station in self.stations:
            stations.append(station.sort_values(**args))
        self.stations = stations
        return self

    def filter_region(self,polygon):
        """
        Filter the region of the stations.

        Parameters:
        -----------
        polygon: list of tuples
            Each tuple is consider a point (lon,lat).
            The first point must be equal to the last point in the polygon.
        
        """
        stations = []
        for station in self.stations:
            stations.append(station.filter_region(polygon))
        self.stations = stations
        return self

    def plot(self,fig=None):
        """
        Plot the stations.

        Parameters:
        -----------
        fig: None or pygmt.Figure
            Basemap figure
        """

        if fig == None:
            fig = pygmt.Figure() 
            fig.basemap(region=self.get_region(padding=0.1),
                        projection="M12c", 
                        frame=["afg","WNse"])

        for station in self.stations:
            station.plot(fig=fig)
        return fig

class Shape():
    def __init__(self,data,projection,**plot_kwargs):
        """
        data: GeoDataFrame
            Data of the shape file
        projection: str
            EPSG projection. Example 'EPSG:4326'
        plot_kwargs: args from Pygmt.plot()
        """
        self.projection =  projection
        self.data = data.to_crs(projection)
        self.plot_kwargs = plot_kwargs

    @property
    def empty(self):
        return self.data.empty

    @property
    def info2pygmt(self):
        args = self.__dict__.copy()
        self.plot_kwargs = update_args2(args,self.plot_kwargs)
        args.update(self.plot_kwargs)
        args = args_cleaner(args,["projection","plot_kwargs"])
        return args.copy()

    def __len__(self):
        return len(self.data)

    def __str__(self,extended=False) -> str:
        msg = f"Shape | {self.__len__()} geometries"
        if extended:
            region = list(map(lambda x: round(x,2),self.get_region()))
            msg += f"\n\tprojection: {self.data.geometry.crs}"
            msg += f"\n\tregion: {region}"
            for i,row in enumerate(self.data.geometry,1):
                msg += f"\n\tgeometry type #{i}: {row.geom_type}"
        else:
            pass
        return msg

    def get_region(self,padding=[]):
        """
        It gets the region according to the limits in the shape
        
        Parameters:
        -----------
        padding: 4D-list or float or int
            list: Padding on each side of the region [lonw,lonw,lats,latn] in degrees.
            float or int: padding amount on each side of the region from 0 to 1,
                        where 1 is considered the distance on each side of the region.
        """
        data = self.data
        df = data["geometry"].apply(lambda x: x.bounds)
        data = pd.DataFrame(df.tolist(),columns=["min_x","min_y","max_x","max_y"])                                    

        lonw,lone = data.min_x.min(),data.max_x.max()
        lats,latn = data.min_y.min(),data.max_y.max()
        
        region = [lonw, lone, lats, latn]
        
        if isinstance(padding,list):
            if padding:
                if len(padding) != 4:
                    raise Exception("Padding parameter must be 4D")
                else:
                    region = list( map(add, region, padding) )
        elif isinstance(padding,float) or isinstance(padding,int):
            lon_distance = abs(region[1]-region[0])
            lat_distance = abs(region[3]-region[2])
            adding4lon = lon_distance*padding
            adding4lat = lat_distance*padding
            padding = [-adding4lon, adding4lon, -adding4lat, adding4lat]
            region = list( map(add, region, padding) )

        return region

    def to_crs(self,projection):
        """
        projection: str
            EPSG projection. Example 'EPSG:4326'
        """
        self.data.to_crs(projection)

    def append(self, data):
        """
        append data
        """
        if isinstance(data, gpd.GeoDataFrame):
            data = data.drop_duplicates(subset=self.columns,ignore_index=True)
            self.data = pd.concat([self.data,data])
        else:
            msg = 'Append only supports a single GeoDataframe object as an argument.'
            raise TypeError(msg)
        return self

    def copy(self):
        """Deep copy of the class"""
        return copy.deepcopy(self)

    def sort_values(self,**args):
        """
        Sort values. Take in mind that it could affect the order of the geometry plotted
        args: The parameters are the pd.DataFrame.sort_values parameters
        """
        self.data = self.data.sort_values(**args)
        return self

    def remove_data(self, rowval):
        """
        remove rows to the data.

        Parameters:
        -----------
        rowval : dict
            key: 
                column name
            value: list
                values specified to remove
        """
        if not isinstance(rowval,dict):
            raise Exception("rowval must be a dictionary")
        
        mask = self.data.isin(rowval)
        mask = mask.any(axis='columns')
        self.data = self.data[~mask]
        return self
    
    def select_data(self, rowval):
        """
        select rows to the data.

        Parameters:
        -----------
        rowval : dict
            key: 
                column name
            value: list
                values specified to select
        """
        if not isinstance(rowval,dict):
            raise Exception("rowval must be a dictionary")
        mask = self.data.isin(rowval)
        mask = mask.any(axis='columns')
        self.data = self.data[mask]
        return self

    def matplot(self,**args):
        """
        args: See GeoDataFrame args.
        """
        self.data.plot(**args)

    def plot(self,fig=None):
        """
        Plot the stations

        Parameters:
        -----------
        fig: None or pygmt.Figure
            Basemap figure
        """

        if fig == None:
            fig = pygmt.Figure() 
            fig.basemap(region=self.get_region(padding=0.1),
                        projection="M12c", 
                        frame=["afg","WNse"])
        fig.plot(**self.info2pygmt)
        return fig

class Shapes():
    def __init__(self,shapes=[]):
        """
        Parameters:
        -----------
        shapes: list
            list of Shape objects
        """
        self.shapes = shapes
    
    def __iter__(self):
        return list(self.shapes).__iter__()

    def __nonzero__(self):
        return bool(len(self.shapes))

    def __len__(self):
        return len(self.shapes)
    
    def __getitem__(self, index):
        if isinstance(index, slice):
            return self.__class__(shapes=self.shapes.__getitem__(index))
        else:
            return self.shapes.__getitem__(index)

    def __delitem__(self, index):
        return self.shapes.__delitem__(index)

    def __getslice__(self, i, j, k=1):
        return self.__class__(shapes=self.shapes[max(0, i):max(0, j):k])
    
    def __str__(self,extended=False) -> str:
        msg = f"Shapes ({self.__len__()} shapes)\n"
        msg += "-"*len(msg) 

        submsgs = []
        for i,shape in enumerate(self.__iter__(),1):
            submsg = f"{i}. "+shape.__str__(extended=extended)
            submsgs.append(submsg)
                
        if len(self.shapes)<=20 or extended is True:
            submsgs = "\n".join(submsgs)
        else:
            three_first_submsgs = submsgs[0:3]
            last_two_subsgs = submsgs[-2:]
            len_others = len(self.shapes) -len(three_first_submsgs) - len(last_two_subsgs)
            submsgs = "\n".join(three_first_submsgs+\
                                [f"...{len_others} other catalogs..."]+\
                                last_two_subsgs)

        return msg+ "\n" +submsgs

    def get_region(self,padding=[]):
        """
        It gets the region according to the limits of all shapes as a whole

        Parameters:
        -----------
        padding: 4D-list or float or int
            list: Padding on each side of the region [lonw,lonw,lats,latn] in degrees.
            float or int: padding amount on each side of the region from 0 to 1,
                        where 1 is considered the distance on each side of the region.
        """
        lons,lats = [],[]
        for shape in self.shapes:
            region = shape.get_region()
            lons.append(region[0:2])
            lats.append(region[2:])
        lons = [ x for sublist in lons for x in sublist]
        lats = [ x for sublist in lats for x in sublist]
        region = [min(lons),max(lons),min(lats),max(lats)]

        if isinstance(padding,list):
            if padding:
                if len(padding) != 4:
                    raise Exception("Padding parameter must be 4D")
                else:
                    region = list( map(add, region, padding) )
        elif isinstance(padding,float) or isinstance(padding,int):
            lon_distance = abs(region[1]-region[0])
            lat_distance = abs(region[3]-region[2])
            adding4lon = lon_distance*padding
            adding4lat = lat_distance*padding
            padding = [-adding4lon, adding4lon, -adding4lat, adding4lat]
            region = list( map(add, region, padding) )

        return region

    def append(self, shapes):
        """
        append a shapes
        """
        if isinstance(shapes, Shape):
            self.shapes.append(shapes)
        else:
            msg = 'Append only supports a single Shape object as an argument.'
            raise TypeError(msg)
        return self
    
    def to_crs(self,projection):
        """
        projection: str
            EPSG projection. Example 'EPSG:4326'
        """
        shapes = []
        for shape in self.shapes:
            shapes.append(shape.to_crs(projection))
        self.shapes = shapes
        return self

    def remove_data(self,rowval):
        """
        remove rows of each data.

        Parameters:
        -----------
        rowval : dict
            key:  
                column name
            value: 
                One or more values specified to remove
        """
        shapes = []
        for shape in self.shapes:
            shapes.append(shape.remove_data(rowval))
        self.shapes = shapes
        return self

    def select_data(self,rowval):
        """
        select rows of each data.

        Parameters:
        -----------
        rowval : dict
            key:  
                column name
            value: 
                One or more values specified to select
        """
        shapes = []
        for shape in self.shapes:
            shapes.append(shape.select_data(rowval))
        self.shapes = shapes
        return self

    def copy(self):
        """Deep copy of the class"""
        return copy.deepcopy(self)  
    
    def sort_values(self,**args):
        """
        Sort values. Take in mind that it could affect the order of the shapes plotted
        args: The parameters are the pd.DataFrame.sort_values parameters
        """
        shapes = []
        for shape in self.shapes:
            shapes.append(shape.sort_values(**args))
        self.shapes = shapes
        return self

    def plot(self,fig=None):
        """
        Plot the stations.

        Parameters:
        -----------
        fig: None or pygmt.Figure
            Basemap figure
        """

        if fig == None:
            fig = pygmt.Figure() 
            fig.basemap(region=self.get_region(padding=0.1),
                        projection="M12c", 
                        frame=["afg","WNse"])

        for shape in self.shapes:
            shape.plot(fig=fig)
        return fig



class FocalMechanism():
    def __init__(self,data,
                text_in_map=True,
                color="red",
                apply_cbar=False,
                scale=1,
                ):
        """
        Parameters:
        -----------
        data: pd.DataFrame 
            Dataframe with the next mandatory columns:
            'origin_time','latitude','longitude','depth','magnitude',
            'strike','dip','rake'.
            Optionals:
                'event_name' to write text in each beachball
                'plot_latitude','plot_longitude' at which to place beachball
        color: str or None
            Color from pygmt color gallery. 
            It is not considered when cbar=True
        apply_cbar: bool
            Use Colorbar (the specifications of the colorbar are located in FocalMechanisms object).
        scale: float
            default:  1cm -> M5
            Adjusts the scaling of the radius of the beachball, 
            which is proportional to the magnitude. 
            Scale defines the size for magnitude = 5.
        """

        self.columns = ['origin_time','latitude','longitude',
                        'depth','magnitude',
                        'strike','dip','rake']
        check =  all(item in data.columns.to_list() for item in self.columns)
        if not check:
            raise Exception("There is not the mandatory columns for the data in Catalog object."\
                            +"->'origin_time','latitude','longitude','depth','magnitude'")

        data = data.drop_duplicates(subset=self.columns,ignore_index=True)
        data["origin_time"] = pd.to_datetime(data["origin_time"]).dt.tz_localize(None)

        self.data = data
        self.text_in_map = text_in_map
        self.color = color
        self.apply_cbar = apply_cbar
        self.scale = scale

        

    @property
    def empty(self):
        return self.data.empty

    def __len__(self):
        return len(self.data)

    def __str__(self) -> str:
        msg = f"FM | {self.__len__()} events "\
                +f"| start:{self.data.origin_time.min()} "\
                +f"| end:{self.data.origin_time.max()}"
        return msg

class Well():
    def __init__(self,data,name,
                color="blue",
                apply_cbar=False,
                injection = pd.DataFrame(),
                injection_cbar = None
                ) -> None:
        """
        Parameters:
        -----------
        data: pd.DataFrame 
            Dataframe with the next mandatory columns:
            'latitude','longitude','z','TVD','MD'
        name: str
            Name of the well
        color: str or None
            Color from pygmt color gallery. 
            It is not considered when cbar=True
        apply_cbar: bool
            Use Colorbar (the specifications of the colorbar are located in Wells object). 
        injection: pd.DataFrame
            Dataframe with the next mandatory columns:
            'min_depth','max_depth','depth_type','water_flow'
        injection_cbar: None or Cbar
            Colorbar for the amount of injection.
        """
        columns = ['latitude','longitude','z','TVD','MD']
        cols = list(set(columns) & set(data.columns.to_list()))
        if list(set(cols)) != list(set(columns)):
            raise Exception("There is not the mandatory columns for the data in Well object."\
                            +"->'latitude','longitude','z','TVD','MD'")
        self.data = data[columns]
        self.color = color
        self.apply_cbar = apply_cbar
        self.name = name
        self.injection = injection
        self.injection_cbar = injection_cbar

        

        @property
        def empty(self):
            return self.data.empty

        def __len__(self):
            return len(self.data)

        def __str__(self) -> str:
            msg = f"Well | {self.__len__()} points in the trajectory"
            return msg

    def plot(self,ax=None):
        if ax == None:
            fig = plt.figure()
            ax = plt.axes(projection='3d')

        ax.plot3D(self.data.longitude, self.data.latitude,
                    self.data.z, 'gray')
        ax.invert_zaxis()

class Profile():
    def __init__(self,
        name, coords, width, 
        colorline="magenta", #only for map figure
        color= "blue", # only for profile figure.
                        #color of the events in the profile. Only if cbar is False
        apply_cbar=True, # only for profile figure. # cbar controlled by cbar_profile_args
        grid=None,
        legend=False,
        ):
        """
        Parameter:
        ----------
        name: str
            Name of the profile.
        coords: 2d-tuple
            2d-Tuple of two 2d-tuples.
            ((ini_lon,ini_lat),(end_lon,end_lat))
        width: 2d-tuple
            (left_width,right_width) in km.
        color: str or None
            Color from pygmt color gallery. 
            It is not considered when cbar=True
        apply_cbar: bool
            Use Colorbar (the specifications of the colorbar
            are located in FocalMechanisms object).
        grid: 2d-tuple
            Grid of the profile plot. (x-axis in km,y-axis in km)
        legend:bool
            True to show the legends
        """
        self.name = name
        self.coords = coords
        self.width = width
        self.colorline = colorline
        self.color = color
        self.apply_cbar = apply_cbar
        self.grid=grid
        self.legend = legend

class Wells():
    def __init__(self,wells=[],cbar = None):
        """
        Parameters:
        -----------
        wells: list
            list of Well objects
        cbar: Cbar object
            Colorbar applied.
        """
        self.wells = wells
        self.cbar = cbar
        
class FocalMechanisms():
    def __init__(self,fms=[],cbar = None):
        """
        Parameters:
        -----------
        fms: list
            list of FocalMechanism objects
        cbar: Cbar object
            Colorbar applied.
        """
        self.fms = fms
        self.cbar = cbar
        
class Profiles():
    def __init__(self,profiles=[]):
        """
        Parameters:
        -----------
        profiles: list
            list of Shape objects
        """
        self.profiles = profiles
        


if __name__ == "__main__":
    cat = Catalog(data="hola")
    print(cat.to_dict())
