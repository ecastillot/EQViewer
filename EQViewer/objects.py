# /**
#  * @author [Emmanuel Castillo]
#  * @email [ecastillot@unal.edu.co]
#  * @create date 2022-11-03 12:07:04
#  * @modify date 2022-11-03 12:07:04
#  * @desc [description]
#  */
import pandas as pd
import matplotlib.dates as mdates
import matplotlib.pyplot as plt

class Viewer():
    def __init__(self,**kwargs):
        self.kwargs = kwargs
    def to_dict(self):
        return self.kwargs

class Catalog(Viewer):
    def __init__(self,
            data,
            size=None,
            style="c0.2c",
            cbar=False,
            color="lightblue",
            label="data",
            transparency=0,
            pen="black",
            ) -> None:

        """
        Parameters:
        -----------
        data: pd.DataFrame 
            Dataframe with the next mandatory columns:
            'origin_time','latitude','longitude','depth','magnitude'
        size: None or lambda function
            Equation for the size. 
            lambda x: 0.1 * np.sqrt(1.5 ** (x*1.5))
            where x always is assigned to the magnitude.
            
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
        cbar: bool
            Use Colorbar (the specifications of the colorbar are located in Catalogs object). 
        color: str or None
            Color from pygmt color gallery. 
            It is not considered when cbar=True
        transparency: float
            transparency of your plots
        pen : str
            color and size of the symbol border
        """

        
        columns = ['origin_time','latitude','longitude','depth','magnitude']
        cols = list(set(columns) & set(data.columns.to_list()))
        if list(set(cols)) != list(set(columns)):
            raise Exception("There is not the mandatory columns for the data in Catalog object."\
                            +"->'origin_time','latitude','longitude','depth','magnitude'")
        data = data.drop_duplicates(subset=columns,ignore_index=True)
        data["origin_time"] = pd.to_datetime(data["origin_time"]).dt.tz_localize(None)
        self.data = data[columns]
        self.color = color
        self.label = label
        self.size = size
        self.style = style
        self.cbar = cbar
        self.transparency = transparency
        self.pen = pen

        l = locals()
        l.pop("self")
        l.pop("__class__")
        super().__init__(**l)

    @property
    def empty(self):
        return self.data.empty

    def __len__(self):
        return len(self.data)

    def __str__(self) -> str:
        msg = f"Catalog | {self.__len__()} events "\
                +f"| start:{self.data.origin_time.min()} "\
                +f"| end:{self.data.origin_time.max()}"
        return msg

    def plot(self,fig):
        if fig == None:
            fig = plt.figure(figsize=(10, 10))

    def mplot(self,color_target="depth",
            s=8,cmap="viridis",
            ax=None):
        if ax == None:
            fig = plt.figure(figsize=(10, 10))
            ax = fig.add_subplot(111)
        ax.set_aspect("equal")


        if color_target == "origin_time":
            cb = ax.scatter(self.data.longitude, self.data.latitude,
                    c=mdates.date2num(self.data[color_target]), s=s, cmap=cmap)
            cbar = fig.colorbar(cb)
            cbar.ax.set_ylim(cbar.ax.get_ylim()[::-1])
            cbar.set_label(f"{color_target}")
            
            loc = mdates.AutoDateLocator()
            cbar.ax.yaxis.set_major_locator(loc)
            cbar.ax.yaxis.set_major_formatter(mdates.ConciseDateFormatter(loc))
        else:
            cb = ax.scatter(self.data.longitude, self.data.latitude,
                    c=self.data[color_target], s=s, cmap=cmap)
            cbar = fig.colorbar(cb)
            cbar.ax.set_ylim(cbar.ax.get_ylim()[::-1])
            cbar.set_label(f"{color_target}")

        ax.set_xlabel("Longitude [°]")
        ax.set_ylabel("Latitude [°]")
        return ax

class Station(Viewer):
    def __init__(self,data,
                name_in_map=False,
                color="black",
                label="stations",
                transparency = 0,
                style="i0.3c",
                pen="black") -> None:

        """
        data: DataFrame
            Dataframe with the next mandatory columns:
            'network','station','latitude','longitude','elevation'
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
        """
        columns = ['network','station','latitude','longitude','elevation']
        cols = list(set(columns) & set(data.columns.to_list()))
        if list(set(cols)) != list(set(columns)):
            raise Exception("There is not the mandatory columns for the data in Station object."\
                            +"->'network','station','latitude','longitude','elevation'")
        self.data = data[columns]
        self.name_in_map = name_in_map
        self.color = color
        self.label = label
        self.style = style
        self.pen = pen
        self.transparency = transparency

        l = locals()
        l.pop("self")
        l.pop("__class__")
        super().__init__(**l)

    @property
    def empty(self):
        return self.data.empty

    def __len__(self):
        return len(self.data)

    def __str__(self) -> str:
        msg = f"Station | {self.__len__()} stations"
        return msg

    def mplot(self,ax=None):
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

class Well(Viewer):
    def __init__(self,data,name,
                color="blue",
                cbar=False,
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
        cbar: bool
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
        self.cbar = cbar
        self.name = name
        self.injection = injection
        self.injection_cbar = injection_cbar

        l = locals()
        l.pop("self")
        l.pop("__class__")
        super().__init__(**l)

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
    #     ax.set_aspect("equal")
    #     cb = ax.scatter(self.data.longitude, self.data.latitude,
    #             c=self.data[color_target], s=8, cmap="viridis")
    #     cbar = fig.colorbar(cb)
    #     cbar.ax.set_ylim(cbar.ax.get_ylim()[::-1])
    #     cbar.set_label(f"{color_target}")

    #     return ax

class FocalMechanism(Viewer):
    def __init__(self,data,
                color="red",
                cbar=False,
                scale_for_m5=1,
                main_n=2,
                ):
        """
        Parameters:
        -----------
        data: pd.DataFrame 
            Dataframe with the next mandatory columns:
            'origin_time','latitude','longitude','depth','magnitude',
            'strike_n1','dip_n1','rake_n1',
           'strike_n2','dip_n2','rake_n2'.
        color: str or None
            Color from pygmt color gallery. 
            It is not considered when cbar=True
        cbar: bool
            Use Colorbar (the specifications of the colorbar are located in FocalMechanisms object).
        scale_for_m5: float
            default: M5 -> 1cm
            Adjusts the scaling of the radius of the beachball, 
            which is proportional to the magnitude. 
            Scale defines the size for magnitude = 5.
        main_n: int
            Could be 1 or 2 depending on the nodal plane.
        """

        columns = ['origin_time','latitude','longitude','depth','magnitude',
                        'strike_n1','dip_n1','rake_n1',
                    'strike_n2','dip_n2','rake_n2']
        cols = list(set(columns) & set(data.columns.to_list()))
        if list(set(cols)) != list(set(columns)):
            raise Exception("There is not the mandatory columns for the data in Catalog object."\
                            +"->'origin_time','latitude','longitude','depth','magnitude',"\
                            +"'strike_n1','dip_n1','rake_n1,'"\
                            +"'strike_n2','dip_n2','rake_n2'")
        data = data.drop_duplicates(subset=columns,ignore_index=True)
        data["origin_time"] = pd.to_datetime(data["origin_time"]).dt.tz_localize(None)
        self.data = data[columns]

        self.data = data
        self.color = color
        self.cbar = cbar
        self.scale_for_m5 = scale_for_m5
        self.main_n = main_n

        l = locals()
        l.pop("self")
        l.pop("__class__")
        super().__init__(**l)

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

class Shape(Viewer):
    def __init__(self,data,**plot_kwargs):
        """
        data: GeoDataFrame

        plot_kwargs: args from Pygmt.plot()
        """
        self.data = data
        plot_kwargs.pop("data",None)
        self.plot_kwargs = plot_kwargs
        l = plot_kwargs
        super().__init__(**l)

    @property
    def empty(self):
        return self.data.empty

    def mplot(self,**args):
        """
        args: See GeoDataFrame args.
        """
        self.data.plot(**args)

class Profile(Viewer):
    def __init__(self,
        name, coords, width, 
        colorline="magenta", #only for map figure
        color= "blue", # only for profile figure.
                        #color of the events in the profile. Only if cbar is False
        cbar=True, # only for profile figure. # cbar controlled by cbar_profile_args
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
        cbar: bool
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
        self.cbar = cbar
        self.grid=grid
        self.legend = legend

        l = locals()
        l.pop("self")
        l.pop("__class__")
        super().__init__(**l)

class Cbar(Viewer):
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
        l = locals()
        l.pop("self")
        l.pop("makecpt_kwargs")
        l.pop("__class__")
        l.update(**makecpt_kwargs)
        super().__init__(**l) 

class Catalogs(Viewer):
    def __init__(self,catalogs=[],cbar = None):
        """
        Parameters:
        -----------
        catalogs: list
            list of Catalog objects
        cbar: Cbar object
            Colorbar applied.
        """
        self.catalogs = catalogs
        self.cbar = cbar
        l = locals()
        l.pop("self")
        l.pop("__class__")
        super().__init__(**l)

    def __iter__(self):
        return list(self.catalogs).__iter__()

    def __nonzero__(self):
        return bool(len(self.catalogs))

    def __len__(self):
        return len(self.catalogs)
    
    def __str__(self) -> str:
        msg = f"Catalogs ({self.__len__()} catalogs)"
        submsgs = []
        for i,catalog in enumerate(self.__iter__(),1):
            submsg = "\t"+catalog.__str__()
            submsgs.append(submsg)
        submsgs = "\n".join(submsgs)
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
        if isinstance(catalog, Catalog):
            self.catalogs.append(catalog)
        else:
            msg = 'Append only supports a single Catalog object as an argument.'
            raise TypeError(msg)
        return self

class Stations(Viewer):
    def __init__(self,stations=[]):
        """
        Parameters:
        -----------
        stations: list
            list of Catalog objects
        """
        self.catalogs = stations
        l = locals()
        l.pop("self")
        l.pop("__class__")
        super().__init__(**l)

class Wells(Viewer):
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
        l = locals()
        l.pop("self")
        l.pop("__class__")
        super().__init__(**l)

class FocalMechanisms(Viewer):
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
        l = locals()
        l.pop("self")
        l.pop("__class__")
        super().__init__(**l)

class Shapes(Viewer):
    def __init__(self,shapes=[]):
        """
        Parameters:
        -----------
        shapes: list
            list of Shape objects
        """
        self.shapes = shapes
        l = locals()
        l.pop("self")
        l.pop("__class__")
        super().__init__(**l)

class Profiles(Viewer):
    def __init__(self,profiles=[]):
        """
        Parameters:
        -----------
        profiles: list
            list of Shape objects
        """
        self.profiles = profiles
        l = locals()
        l.pop("self")
        l.pop("__class__")
        super().__init__(**l)


if __name__ == "__main__":
    cat = Catalog(data="hola")
    print(cat.to_dict())
