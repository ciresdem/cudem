### factory.py - CUDEM module factory
##
## Copyright (c) 2023 - 2025 Regents of the University of Colorado
##
## factory.py is part of CUDEM
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
### Commentary:
##
## A demonstration Hello World Factory Module template.
##
### Code:

import sys
from cudem import factory

## Define your custom module
class HelloWorldModule(factory.CUDEMModule):
    """A simple module that prints a greeting."""
    
    def run(self):
        """The factory calls this method automatically after initialization.
        Parameters passed in the string become attributes of 'self'.
        """
        
        ## Set defaults if the user didn't provide them
        msg = getattr(self, 'msg', 'Hello World')
        count = getattr(self, 'count', 1)
        
        ## Ensure count is an integer (factory parses numbers as strings by default)
        count = int(count)
        
        print(f"\n--- Running Module: {self.params.get('mod')} ---")
        for i in range(count):
            print(f"{i+1}: {msg}")
        print("-----------------------------------\n")

        
## Define your Factory
class MyFactory(factory.CUDEMFactory):
    """A factory that knows about the HelloWorldModule."""
    
    ## Register the module here. 
    ## The key 'greet' is what we will use in the module string.
    _modules = {
        'greet': {
            'name': 'Greeter',
            'description': 'Prints a message multiple times',
            'call': HelloWorldModule  # The class to instantiate
        }
    }

    
## Execution examples
if __name__ == "__main__":
    
    print("Test 1: Basic Usage")
    ## String: "greet" -> Runs with defaults
    f1 = MyFactory(mod="greet")._acquire_module()
    f1.run()
    
    print("Test 2: Passing Arguments")
    ## String: "greet:msg=Custom Message:count=3"
    ## Note: 'msg' and 'count' are passed to the module instance
    f2 = MyFactory(mod="greet:msg=Custom Message:count=3")._acquire_module()
    f2.run()

    print("Test 3: Handling Lists and Types")
    ## String: "greet:msg=List Test:count=1:tags=A;B;C"
    ## The factory automatically converts "A;B;C" into a python list ['A','B','C']
    f3 = MyFactory(mod="greet:msg=List Test:count=1:tags=A;B;C")._acquire_module()

    ## We can inspect the instance created to see the list handling
    if hasattr(f3, 'tags'):
        print(f"Parsed 'tags' as list: {f3.tags}")


### End
