### checkpoints_3dep.py
##
## Copyright (c) 2026 - 2026 Regents of the University of Colorado
##
## checkpoints_3dep.py is part of CUDEM
##
### Commentary:
##
## Fetch 3dep Checkpoint Data
##
### Code:

from cudem import utils
from cudem.fetches import fetches

## ==============================================
## Constants
## ==============================================
CHECKPOINTS_3DEP_URL = 'https://www.sciencebase.gov/catalog/file/get/67075e6bd34e969edc59c3e7?f=__disk__80%2F12%2F9e%2F80129e86d18461ed921b288f13e08c62e8590ffb'
REFERER = 'https://www.sciencebase.gov/vocab/category/item/identifier'

HEADERS = {
    'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:89.0) Gecko/20100101 Firefox/89.0',
    'referer': REFERER
}

class CheckPoints3DEP(fetches.FetchModule):
    def __init__(self, **kwargs):
        super().__init__(name='3dep_cp', **kwargs)
        self.headers = HEADERS
        
    def run(self):
        self.add_entry_to_results(
            CHECKPOINTS_3DEP_URL,
            'CheckPoints_3DEP.zip',
            'checkpoints'
        )
        
