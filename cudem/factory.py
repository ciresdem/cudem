### factory.py - CUDEM utilities and functions
##
## Copyright (c) 2023 Regents of the University of Colorado
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
## The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
##
## THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
## INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR 
## PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE 
## FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, 
## ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
##
### Commentary:
##
## Module factory...
##
## Use factory.CUDEMFactory as a superclass for a specific module factory, and define
## a _modules dict that has the module name as keys and their values are a dict with atleast
## a 'call' key which points to a module sub-class
##
## class PerspectoFactory(factory.CUDEMFactory):
##     if has_pygmt:
##         _modules = {
##             'hillshade': {'call': Hillshade},
##             'perspective': {'call': perspective},
##             'sphere': {'call': sphere},
##             'figure1': {'call': figure1},        
##         }
##     else:
##         _modules = {
##             'hillshade': {'call': Hillshade},
##             'perspective': {'call': perspective},
##             'sphere': {'call': sphere},
##         }    
##     def __init__(self, **kwargs):
##         super().__init__(**kwargs)
##
##
## In the superclass we need at least a 'params' parameter, this dict holds the factory parameters, including
## 'mod', 'mod_args', 'kwargs'
##
## class MySuperClass:
##     def __init__(self, params={}):
##         self.params = params
##
## To add module paramers when module is run outside of the factory (only if params dict is
## wanted, though it could be easier to just use the factory if that's the case):
## make sure all super and sub classes __init__ functions only define parameters.
## Make use of a self.initialize() function or something if needed, where you first
## set the 'mod_args' parameters before defining anything futher...
##
## from cudem import factory
##
## # in __init__ of the superclass do:
## def __init__(self, p = None, params = {}, **kwargs):
##     self.p = p
##     self.params = params
##     factory._set_params(self, 'mod', 'mod_name', {})
##     self.initialize()
##
## # then, after a module has been loaded:
## def initialize(self):
##     factory._set_mod_args(self)
##     if not self.params['mod_args']:
##         for k in self.__dict__:
##             if k not in self.params['kwargs'].keys():
##                 self.params['mod_args'][k] = self.__dict__[k]
### Code:

import os
import sys
import json
from cudem import utils

def args2dict(args, dict_args={}):
    """convert list of arg strings to dict.
    
    Args:
      args (list): a list of ['key=val'] pairs
      dict_args (dict): a dict to append to

    Returns:
      dict: a dictionary of the key/values
    """
    
    for arg in args:
        #this_entry = re.findall(r'[^"\s]\S*|".+?"', arg)
        p_arg = arg.split('=')
        if len(p_arg) > 1:
            dict_args[p_arg[0]] = False if p_arg[1].lower() == 'false' else \
                True if p_arg[1].lower() == 'true' else \
                None if p_arg[1].lower() == 'none' else \
                '='.join(p_arg[1:]) if len(p_arg) > 2 else \
                p_arg[1]
        
    return(dict_args)

def dict2args(in_dict):
    out_args = ''
    for i, key in enumerate(in_dict.keys()):
        out_args += '{}={}{}'.format(key, in_dict[key], ':' if i+1 < len(in_dict.keys()) else '')
    return(out_args)

def _set_params(mc, mod=None, mod_name=None, mod_args={}):
    if 'params' not in mc.__dict__.keys():
        mc.params = {}
        
    if 'mod' not in mc.params.keys():
        mc.params['mod'] = mod

    if 'mod_name' not in mc.params.keys():
        mc.params['mod_name'] = mod_name

    if 'mod_args' not in mc.params.keys():
        mc.params['mod_args'] = {}

    if 'kwargs' not in mc.params.keys():
        mc.params['kwargs'] = mc.__dict__.copy()


def _set_mod_params(mc, mf=None, mod=None, mod_name=None):
    if mod is not None:
        mc.params['mod'] = mod
        
    if mod_name is not None:
        mc.params['mod_name'] = mod_name

    if mc.params['mod'] is None:
        mc.params['mod'] = _factory_module_check(mf, mc)

    if mc.params['mod_name'] is None:
        mc.params['mod_name'] = mc.params['mod']
        
    if not mc.params['mod_args']:
        for k in mc.__dict__:
            if k not in mc.params['kwargs'].keys():
                mc.params['mod_args'][k] = mc.__dict__[k]

def _factory_module_check(mf, mc):
    for k in mf._modules.keys():
        if isinstance(mf._modules[k]['call'](), mc):
            return(k)
                
## ==============================================
## echo cudem module options
## modules are a dictionary with the module name
## as the key and at least a 'class' key which
## points to the class/function to call for the module
## uses <class>.__doc__ as description
##
## e.g.
## _cudem_module_long_desc({'module_name': {'class': MyClass}})
## ==============================================
_cudem_module_short_desc = lambda m: ', '.join(
    ['{}'.format(key) for key in m])
_cudem_module_name_short_desc = lambda m: ',  '.join(
    ['{} ({})'.format(m[key]['name'], key) for key in m])
_cudem_module_long_desc = lambda m: '{cmd} modules:\n% {cmd} ... <mod>:key=val:key=val...\n\n  '.format(cmd=os.path.basename(sys.argv[0])) + '\n  '.join(
    ['\033[1m{:14}\033[0m{}\n'.format(str(key), m[key]['call'].__doc__) for key in m]) + '\n'

def echo_modules(module_dict, key):
    if key is None:
        sys.stderr.write(_cudem_module_long_desc(module_dict))
    else:
        if key in module_dict.keys():
            sys.stderr.write(
                _cudem_module_long_desc(
                    {k: module_dict[k] for k in (key,)}
                )
            )
        else:
            sys.stderr.write('Invalid Module Key: {}\nValid Modules: {}\n'.format(key, _cudem_module_short_desc(module_dict)))

    sys.stderr.flush()

## ==============================================
##
## Factories
##
## ==============================================    
class CUDEMModule:
    def __init__(self, params = {}, **kwargs):
        self.params = params
        for kpam, kval in kwargs.items():
            self.__setattr__(kpam, kval)

        self.run()

    def initialize(self):
        pass
            
    def __call__(self):
        print(self.__dict__)
        self.run()

    def run(self):
        raise NotImplementedError('the module {} could not be parsed by the sub factory'.format(self.params['mod']))

class CUDEMFactory:
    ## optional modules, should be dictionary of dictionaries
    ## where each key is the module name and has at least a 'call' key pointing
    ## to a function/class with a __call__ functions defined
    ## e.g from waffles { 'surface': {'name': 'surface', 'datalist-p': True, 'call':waffles.GMTSurface}, ...
    ## the function/class to call should have at least a 'params={}' paramter.
    _factory_module = {'_factory': {'call': CUDEMModule}}
    _modules = {'_factory': {'call': CUDEMModule}}
    
    def __init__(self, mod: str = None, **kwargs):
        """
        Initialize the factory default settings
        
        ----------
        Parameters

        mod - A string of a module name and optional module arguments in the format: 'mod_name:mod_arg0=mod_val0:mod_arg1=mod_val1'
        params - A parameters object. Should hold all factory global arguments
        """
        
        self.mod = mod
        self.mod_name = '_factory'
        self.mod_args = {}
        self.kwargs = kwargs
        
        if self.mod is not None:
            self._parse_mod(self.mod)
        #elif self.kwargs['fn'] is not None:
        #    self._parse_mod(self.kwargs['fn'])

    def __str__(self):
        return('<{}>'.format(self.__dict__))
    
    def __repr__(self):
        return('<{}>'.format(self.__dict__))
    
    def _parse_mod(self, mod):
        #self.add_module(self._factory_module)
        opts = mod.split(':')
        if opts[0] in self._modules.keys():
            if len(opts) > 1:
                self.mod_args = args2dict(list(opts[1:]), {})
            self.mod_name = opts[0]
        else:
            utils.echo_error_msg(
                'invalid module name `{}`'.format(opts[0])
            )
            return(None)
        return(self.mod_name, self.mod_args)

    def add_module(self, type_def={}):
        for key in type_def.keys():
            self._modules[key] = type_def[key]
    
    def _acquire_module(self):
        ## only put params in if its a defined parameter
        ##    if 'params' in self._modules[self.mod_name]['call']().__dict__.keys():
        ## KeyError: '_factory'
        #if 'params' in self._modules[self.mod_name]['call']().__dict__.keys():
        self.kwargs['params'] = self.__dict__

        m = lambda m, k: self._modules[self.mod_name]['call'](**m, **k)
        return(m(self.mod_args, self.kwargs))
        #return(self._modules[self.mod_name]['call'](self.mod_args, self.kwargs))

    def load_parameter_dict(self, param_dict: dict):
        valid_data = False

        for ky, val in param_dict.items():
            #if ky in self.__dict__:
            self.__setattr__(ky, val)
            valid_data = True
                    
        if valid_data:
            utils.echo_msg('CUDEMFactory read successfully')
        else:
            utils.echo_warning_msg('Unable to find any valid data')
            
        return(self)
        
    def write_parameter_file(self, param_file: str):
        try:
            with open(param_file, 'w') as outfile:
                json.dump(self.__dict__, outfile, indent=4)
                utils.echo_msg('New CUDEMFactory file written to {}'.format(param_file))
                
        except:
            raise ValueError('CUDEMFactory: Unable to write new parameter file to {}'.format(param_file))
        
    def open_parameter_file(self, param_file: str):
        valid_data = False
        with open(param_file, 'r') as infile:
            try:
                data = json.load(infile)
            except:
                raise ValueError('CUDEMFactory: Unable to read data from {} as json'.format(param_file))
            
            for ky, val in data.items():
                #if ky in self.__dict__:
                self.__setattr__(ky, val)
                valid_data = True
                    
            if valid_data:
                utils.echo_msg('CUDEMFactory read successfully from {}'.format(param_file))
            else:
                utils.echo_warning_msg('Unable to find any valid data in {}'.format(param_file))


if __name__ == 'main':
    f = CUDEMFactory()
    
### End
