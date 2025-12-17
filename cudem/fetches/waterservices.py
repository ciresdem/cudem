### waterservices.py
##
## Copyright (c) 2010 - 2025 Regents of the University of Colorado
##
## fetches.py is part of CUDEM
##
## Permission is hereby granted, free of charge, to any person obtaining a copy 
## of this software and associated documentation files (the "Software"), to deal 
## in the Software without restriction, including without limitation the rights 
## to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies 
## of the Software, and to permit persons to whom the Software is furnished to do so, 
## subject to the following conditions:
##
## The above copyright notice and this permission notice shall be included in all
## copies or substantial portions of the Software.
##
## THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
## INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR 
## PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE 
## FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, 
## ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
## SOFTWARE.
##
###############################################################################
### Commentary:
##
##
### Code:

from cudem.fetches import fetches

## WaterServices (USGS)
class WaterServices(fetches.FetchModule):
    """WaterServices station information from USGS

    https://waterservices.usgs.gov/

    < waterservices:printout=False >
    """
    
    def __init__(self, printout = False, **kwargs):
        super().__init__(name='waterservices', **kwargs)
        self.printout = printout
        
        ## The various waterservices URLs
        self._water_services_api_url = 'https://waterservices.usgs.gov/nwis/iv/?'
        
    def run(self):
        '''Run the WATERSERVICES fetching module'''
        
        if self.region is None:
            return([])
        
        _data = {
            'bBox': self.region.format('bbox'),
            'siteStatus': 'active',
            'format':'json',
        }
        _req = fetches.Fetch(
            self._water_services_api_url, verbose=self.verbose
        ).fetch_req(params=_data)
        if _req is not None:
            self.add_entry_to_results(
                _req.url,
                'water_services_results_{}.json'.format(self.region.format('fn')),
                'waterservices'
            )

            ## print out the water-station information
            if self.printout:
                j = _req.json()
                print(j.keys())
                out_list = j['value']['timeSeries']
                for l in out_list:
                    # geolocation
                    gl = l['sourceInfo']['geoLocation']['geogLocation']
                    y = gl['latitude']
                    x = gl['longitude']

                    # vars
                    #vs = l['variable']['variableCode']
                    v = l['variable']['variableCode'][0]['value']
                    vc = l['variable']['variableName']
                    #print(v)

                    # val
                    vl = l['values'][0]['value'][0]['value']

                    print(vc, v, x, y, vl)
            
        return(self)

### End
