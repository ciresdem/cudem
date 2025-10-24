from cudem.fetches import fetches

class ShallowBathyEverywhere(fetches.FetchModule):
    """Shallow Bathy Everywhere

    Derived bathymetry for various locations

    < sbe >
    """
    
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.sbe_dem_url = 'https://shallowbathymetryeverywhere.com/data/dem/'
