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
            cmap=False,
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
        cmap: bool
            Use Colorbar (the specifications of the colorbar are located in the main function called seismic_profile.map-> cmap_args). 
        color: str or None
            Color from pygmt color gallery. 
            It is not considered when cmap=True
        transparency: float
            transparency of your plots
        pen : str
            color and size of the symbol border
        """
        self.data = data.drop_duplicates(subset=['origin_time','latitude','longitude','depth','magnitude']
                    ,ignore_index=True)
        self.color = color
        self.label = label
        self.size = size
        self.style = style
        self.cmap = cmap
        self.transparency = transparency
        self.pen = pen

        l = locals()
        l.pop("self")
        l.pop("__class__")
        super().__init__(**l)

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
            Put the name of the map in the figure
        color: str or None
            Color from pygmt color gallery. 
            It is not considered when cmap=True
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

class Well(Viewer):
    def __init__(self,data,
                color="blue,0.1p",
                cmap=False) -> None:
        self.data = data
        self.color = color
        self.cmap = cmap

        l = locals()
        l.pop("self")
        l.pop("__class__")
        super().__init__(**l)

    def add_inyection(self,data):
        self.inyection = data

class FocalMechanism(Viewer):
    def __init__(self,data,
                color="red",
                cmap=False,
                scale_for_m5=1,
                main_n=2,
                ):
        self.data = data
        self.color = color
        self.cmap = cmap
        self.scale_for_m5 = scale_for_m5
        self.main_n = main_n

        l = locals()
        l.pop("self")
        l.pop("__class__")
        super().__init__(**l)

class Shape(Viewer):
    def __init__(self,**plot_kwargs):
        """
        geodata: str
            It is the geodatafrane of the shape
        plot_kwargs: dict
            kwargs to plot 
            key: Column name
            val: Values 
        """

        # self.__dict__.update(**plot_kwargs)

        l = plot_kwargs
        super().__init__(**l)

class Profile(Viewer):
    def __init__(self,name, coords, width, 
        colorline="magenta", #only for map figure
        color= "blue", # only for profile figure.
                        #color of the events in the profile. Only if cmap is False
        # size=lambda x= 0.1 * np.sqrt(1.5 ** (x*1.5)), # only for profile figure.
        size=None, # only for profile figure.
        style ="c0.1c",# only for profile figure.
        pen="black", # only for profile figure. 
        grid=None,
        legend=False,
        cmap=True, # only for profile figure. # cmap controlled by cbar_profile_args
        cbar_profile_args = {"cmap":'viridis', 
                                "color_target":"depth",
                                "label":"Depth (m)",
                                "reverse":False,
                                "series":[0, 3e3] } # only for profile figure. # cbar profile args
        ):
        # {"name":("A","A'"),      
        # "coords":((-73.684091,3.871559),(-73.6777688,3.8750540)),
        # "width":(-0.15,0.15),
        # "colorline":"black",
        # "color": None,
        # "size":lambda x: 0.1 * np.sqrt(1.5 ** (x*1.5)),
        # "style" :"cc",
        # "grid":[200,200],
        # "legend":False,
        # "pen":None,
        # "cmap":True,
        # "cbar_profile_args" : {"cmap":'buda',
        #                         "color_target":"origin_time",
        #                         "label":"Time",
        #                         "overrule_bg":True,
        #                         "reverse":False,
        #                         "series":[dt.datetime(2019,1,1), 
        #                                 dt.datetime(2022,5,1)] }
        # }
        self.name = name
        self.coords = coords
        self.width = width
        self.colorline = colorline
        self.color = color
        self.size = size
        self.style = style
        self.pen = pen
        self.cmap = cmap
        self.grid=grid
        self.legend = legend
        self.cbar_profile_args = cbar_profile_args 

        l = locals()
        l.pop("self")
        l.pop("__class__")
        super().__init__(**l)

        


if __name__ == "__main__":
    #test
    # v = Viewer(hola="jeje")
    # print(v.to_dict())

    cat = Catalog(data="hola")
    print(cat.to_dict())