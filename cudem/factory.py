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
###############################################################################
### Commentary:
##
## Module factory...
##
## Use factory.CUDEMFactory as a superclass for a specific module factory, and
## define a _modules dict that has the module name as keys and their values are
## a dict with atleast
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
## In the superclass we need at least a 'params' parameter, this dict holds
## the factory parameters, including
## 'mod', 'mod_args', 'kwargs'
##
## class MySuperClass:
##     def __init__(self, params={}):
##         self.params = params
##
## To add module paramers when module is run outside of the factory (only if
## params dict is wanted, though it could be easier to just use the factory
## that's the case):
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
import re
import json
from cudem import utils


## a factory module string is
## 'mod_name:mod_arg=arg_val:mod_arg1=arg_val:sub-mod_name="sub-mod_arg=sub_mod_val'
def parse_fmod(fmod):
    opts = fmod2dict(fmod, {})
    
    mod = opts['_module']
    mod_args = {i:opts[i] for i in opts if i!='_module'}

    return(opts, mod, mod_args)


def fmod2dict(fmod, dict_args: dict = {}):
    """convert factory module string to a dict
    
    Parameter:
      fmod (str): a facotory module string
      dict_args (dict): a dict to append to

    Returns:
      dict_args: a dictionary of the key/values
    """

    args_list = re.split(r':(?=(?:[^"]*"[^"]*")*[^"]*$)', fmod)
    for arg in args_list:
        p_arg = re.split(
            r'=(?=(?:[^"]*"[^"]*")*[^"]*$)',
            arg
        )
        if len(p_arg) == 1:
            if '_module' not in dict_args.keys():
                dict_args['_module'] = p_arg[0]
                
        elif len(p_arg) > 1:
            dict_args[p_arg[0]] = False \
                if p_arg[1].lower() == 'false' \
                   else True if p_arg[1].lower() == 'true' \
                        else None if p_arg[1].lower() == 'none' \
                             else '='.join(p_arg[1:]) if len(p_arg) > 2 \
                                  else p_arg[1].strip('"').split(';') if ';' in p_arg[1] \
                                       else p_arg[1].strip('"')
        
    return(dict_args)


def dict2fmod(in_dict: dict):
    """convert a dict of key:val pairs to a module factory string

    Parameter:
      in_dict (dict): the dictionary to convert

    Returns:
      out_args (str): a string representatino of in_dict, 
                      suitable for a factory cli
    """
    
    out_args = ''
    for i, key in enumerate(in_dict.keys()):
        if key == '_module':
            out_args += f'{in_dict[key]}:'
        elif isinstance(in_dict[key], list):
            out_args += f'{key}="{";".join(in_dict[key])}"'
        else:
            
            out_args += '{}={}{}'.format(
                key, in_dict[key], ':' if i+1 < len(in_dict.keys()) else ''
            )
        
    return(out_args)


def args2dict(args, dict_args: dict = {}):
    """convert list of arg strings to dict.
    
    Parameter:
      args (list): a list of ['key=val'] pairs
      dict_args (dict): a dict to append to

    Returns:
      dict_args: a dictionary of the key/values
    """

    for arg in args:
        #this_entry = re.findall(r'[^"\s]\S*|".+?"', arg)
        p_arg = arg.split('=')
        if len(p_arg) > 1:
            dict_args[p_arg[0]] = False \
                if p_arg[1].lower() == 'false' \
                   else True if p_arg[1].lower() == 'true' \
                        else None if p_arg[1].lower() == 'none' \
                             else '='.join(p_arg[1:]) if len(p_arg) > 2 \
                                  else p_arg[1]
        
    return(dict_args)


def dict2args(in_dict: dict):
    """convert a dict of key:val pairs to a module factory string

    Parameter:
      in_dict (dict): the dictionary to convert

    Returns:
      out_args (str): a string representatino of in_dict, 
                      suitable for a factory cli
    """
    
    out_args = ''
    for i, key in enumerate(in_dict.keys()):
        out_args += '{}={}{}'.format(
            key, in_dict[key], ':' if i+1 < len(in_dict.keys()) else ''
        )
        
    return(out_args)


def _set_params(
        mc: any,
        mod: any = None,
        mod_name: str = None,
        mod_args: dict = {}
):
    """set module parameters into mc module class"""
        
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

        
def _set_mod_params(
        mc: any,
        mf: any = None,
        mod: any = None,
        mod_name: str = None
):
    """set module parameters into mc module class"""
        
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

                
def _factory_module_check(mf: any, mc: any):
    """check if factory module exists in the module class

    Returns:
      k: module key
    """
    
    for k in mf._modules.keys():
        if isinstance(mf._modules[k]['call'](), mc):
            return(k)

        
###############################################################################
## echo cudem module options
## modules are a dictionary with the module name
## as the key and at least a 'class' key which
## points to the class/function to call for the module
## uses <class>.__doc__ as description
##
## e.g.
## _cudem_module_long_desc({'module_name': {'class': MyClass}})
###############################################################################
_cudem_module_short_desc = lambda m: ', '.join(
    ['{}'.format(key) for key in m])
_cudem_module_name_short_desc = lambda m: ',  '.join(
    ['{} ({})'.format(
        m[key]['name'] if 'name' in m[key].keys() else None, key
    ) for key in m]
)
_cudem_module_long_desc \
    = lambda m: '{cmd} modules:\n% {cmd} ... <mod>:key=val:key=val...\n\n  '.format(
        cmd=os.path.basename(sys.argv[0])
    ) + '\n  '.join(
        ['\033[1m{:20}\033[0m{}\n'.format(
            '{} ({})'.format(
                str(key),
                str(m[key]['name']) if 'name' in m[key].keys() else key),
            m[key]['call'].__doc__
        ) for key in m]
    ) + '\n'
_cudem_module_md \
    = lambda m: '# {cmd} modules:\n% {cmd} ... <mod>:key=val:key=val...\n\n'.format(
        cmd=os.path.basename(sys.argv[0])
    ) + '\n'.join(
        ['## {} ({})\n {}'.format(
            str(m[key]['name']), str(key), m[key]['description']
        ) for key in m]
    )
_cudem_module_md_table \
    = lambda m: '| **Name** | **Module-Key** | **Description** |\n|---|---|---|\n'.format(
        cmd=os.path.basename(sys.argv[0])
    ) + '\n'.join(
        ['| {} | {} | {} |'.format(
            str(m[key]['name']) if 'name' in m[key].keys() else key,
            str(key), m[key]['description']
        ) for key in m]
    )


def echo_modules(module_dict: dict, key: any, md: bool = False):
    """print out the existing modules from module_dict and 
    their descriptions.
    """
    
    if key is None:
        if not md:
            sys.stderr.write(_cudem_module_long_desc(module_dict))
        else:
            print(_cudem_module_md_table(module_dict))
    else:
        if key in module_dict.keys():
            sys.stderr.write(
                _cudem_module_long_desc(
                    {k: module_dict[k] for k in (key,)}
                )
            )
        else:
            for mk in module_dict.keys():
                if key == module_dict[mk]['name']:
                    key = mk
                    break

            if key in module_dict.keys():
                sys.stderr.write(
                    _cudem_module_long_desc(
                        {k: module_dict[k] for k in (key,)}
                    )
                )
            else:
                sys.stderr.write(
                    'Invalid Module Key: {}\nValid Modules: {}\n'.format(
                        key, _cudem_module_name_short_desc(module_dict)
                    )
                )

    sys.stderr.flush()

    
class CUDEMModule:
    """cudem factory module.

    Template for a compatible factory module setup, together with 
    CUDEMFactory.
    """
    
    def __init__(self, modules = {}, params = {}, **kwargs):
        self.params = params
        for kpam, kval in kwargs.items():
            self.__setattr__(kpam, kval)

        self.kwargs = kwargs
        self.modules = modules
        self.run()

    def initialize(self):
        pass
            
    def __call__(self):
        print(self.__dict__)
        self.run()

    def run(self):
        raise NotImplementedError(
            (f'the module {self.params["mod"]} could not be parsed by the sub factory, '
             f'available modules are: {_cudem_module_short_desc(self.modules)}')
        )

    
# params - A parameters object. Should hold all factory global arguments
class CUDEMFactory:
    """cudem module factory.

    Template for a compatible module factory setup, together with 
    CUDEMModule
    
    optional modules, should be dictionary of dictionaries
    where each key is the module name and has at least a 
    'name', 'description' and 'call' key pointing to a 
    function/class with a __call__ functions defined

    e.g from waffles 
    { 'surface': {'name': 'surface', 'datalist-p': True, 'call':waffles.GMTSurface}, ...

    the function/class to call should have at least a 'params={}' paramter.
    """
    
    _factory_module = {
        '_factory': {'name': 'factory',
                     'description': 'default factory setting',
                     'fmts': [],
                     'call': CUDEMModule}
    }
    _modules = {
        '_factory': {'name': 'factory',
                     'description': 'default factory setting',
                     'call': CUDEMModule}
    }

    
    def __init__(self, mod: str = None, **kwargs: any):
        """
        Initialize the factory default settings
        
        Parameters:
          mod - A string of a module name and optional module arguments in the format: 
                'mod_name:mod_arg0=mod_val0:mod_arg1=mod_val1'
          kwargs - module arguments can be passed as key-word arguments 
                   here instead of in the mod string if wanted;
                   however argumets from the mod string will over-ride arguments 
                   from the key-word arguments.
        """
        
        self.mod = mod
        self.mod_name = '_factory'
        self.mod_args = {}
        self.kwargs = kwargs
        if self.mod is not None:
            self._parse_mod(self.mod)

        self.add_module(self._factory_module)

        
    def __str__(self):
        return('<{}>'.format(self.__dict__))

    
    def __repr__(self):
        return('<{}>'.format(self.__dict__))


    def _parse_mod(self, mod: str):
        """parse the module string.

        Returns:
          (module-name, module-arguments)
        """
        
        opts = fmod2dict(mod, {})
        if opts['_module'] in self._modules.keys():
            self.mod_name = opts['_module']
            self.mod_args = {i:opts[i] for i in opts if i!='_module'}
            
        else:
            utils.echo_error_msg(
                'invalid module name `{}`'.format(opts['_module'])
            )            
            self.mod_args['modules'] = self._modules
            
        return(self.mod_name, self.mod_args)


    def _parse_mod_old(self, mod: str):
        """parse the module string.

        Returns:
          (module-name, module-arguments)
        """

        opts = mod.split(':')
        self.mod_args = {}
        if opts[0] in self._modules.keys():
            if len(opts) > 1:
                self.mod_args = args2dict(list(opts[1:]), {})
            self.mod_name = opts[0]
        else:
            utils.echo_error_msg(
                'invalid module name `{}`'.format(opts[0])
            )            
            self.mod_args['modules'] = self._modules
            
        return(self.mod_name, self.mod_args)

    
    def add_module(self, type_def: dict = {}):
        """Add a module to the factory `_modules` dict"""

        for key in type_def.keys():
            self._modules[key] = type_def[key]

            
    def _acquire_module(self):
        """Acquire the module from the factory.

        only put params in if its a defined parameter
           if 'params' in self._modules[self.mod_name]['call']().__dict__.keys():
        KeyError: '_factory'
        #if 'params' in self._modules[self.mod_name]['call']().__dict__.keys():
        """

        def _get_all_keys(d):
            for key, value in d.items():
                if key != 'params':
                    yield key
                else:
                    continue
                
                if isinstance(value, dict):
                    yield from _get_all_keys(value)
        
        self.kwargs['params'] = self.__dict__
        for k in self.kwargs.keys():
            if k in self.mod_args.keys():
                # utils.echo_warning_msg(
                #     'duplicate options! {}: {} --> {}'.format(
                #         k, self.kwargs[k], self.mod_args[k]
                #     )
                # )
                self.kwargs[k] = self.mod_args[k]
                del self.mod_args[k]
                
        if self.mod_name is not None:
            try:
                m = lambda m, k: self._modules[self.mod_name]['call'](**m, **k)
                mm = m(self.mod_args, self.kwargs)
                        
                for mod_arg in self.mod_args.keys():
                    if mod_arg not in list(_get_all_keys(mm.__dict__)):
                        utils.echo_warning_msg(
                            f'{mod_arg} is not a valid parameter...'
                        )
                
                #return(m(self.mod_args, self.kwargs))
                return(mm)
            except Exception as e:
                utils.echo_error_msg(
                    'could not acquire module, {}'.format(e)
                )
                
        #return(self._modules[self.mod_name]['call'](self.mod_args, self.kwargs))

        
    def load_parameter_dict(self, param_dict: dict):
        """load a module parameter dictionary. 

        Can be used to read a stored json representation of the module.
        """
        
        valid_data = False

        for ky, val in param_dict.items():
            #if ky in self.__dict__:
            self.__setattr__(ky, val)
            valid_data = True
                    
        if valid_data:
            utils.echo_msg(
                'CUDEMFactory read successfully'
            )
        else:
            utils.echo_warning_msg(
                'Unable to find any valid data'
            )
            
        return(self)

    
    def write_parameter_file(self, param_file: str):
        """write a module parameter file to disk."""
        
        try:
            with open(param_file, 'w') as outfile:
                json.dump(self.__dict__, outfile, indent=4)
                utils.echo_msg(
                    f'New CUDEMFactory file written to {param_file}'
                )
                
        except:
            raise ValueError(
                f'CUDEMFactory: Unable to write new parameter file to {param_file}'
            )

        
    def open_parameter_file(self, param_file: str):
        """Open and read a saved parameter file"""
        
        valid_data = False
        with open(param_file, 'r') as infile:
            try:
                data = json.load(infile)
            except:
                raise ValueError(
                    f'CUDEMFactory: Unable to read data from {param_file} as json'
                )
            
            for ky, val in data.items():
                #if ky in self.__dict__:
                self.__setattr__(ky, val)
                valid_data = True
                    
            if valid_data:
                utils.echo_msg(
                    f'CUDEMFactory read successfully from {param_file}'
                )
            else:
                utils.echo_warning_msg(
                    f'Unable to find any valid data in {param_file}'
                )

                
if __name__ == 'main':
    f = CUDEMFactory()
    
### End
