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
##
### Code:

import os
import sys
import re
import json
from typing import Any, Dict, List, Optional, Tuple, Union

try:
    from cudem import utils
except ImportError:
    pass  # Expecting calling environment to have cudem installed


## ==============================================
## Utility Functions
## ==============================================
def _parse_value_string(val_str: str) -> Any:
    """Helper to parse string values into Python types (bool, None, list)."""
    
    val_lower = val_str.lower()
    # if utils.str2bool(val_str) is not None:
    #     return utils.str2bool(val_str)
    if val_lower == 'false':
        return False
    elif val_lower == 'true':
        return True
    elif val_lower == 'none':
        return None
    elif ';' in val_str:
        return val_str.strip('"').split(';')
    else:
        return val_str.strip('"')


def parse_fmod(fmod: str) -> Tuple[Dict[str, Any], str, Dict[str, Any]]:
    """Parse a factory module string.
    
    Returns:
        Tuple containing (all_options, module_name, module_arguments)
    """
    
    opts = fmod2dict(fmod)
    mod = opts.get('_module')
    mod_args = {k: v for k, v in opts.items() if k != '_module'}
    return opts, mod, mod_args


def fmod2dict(fmod: str, dict_args: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    """Convert factory module string to a dict.

    Args:
      fmod (str): A factory module string.
      dict_args (dict, optional): A dict to append to.

    Returns:
      dict: A dictionary of the key/values.
    """
    
    if dict_args is None:
        dict_args = {}

    ## Split by colon, ignoring colons inside quotes
    args_list = re.split(r':(?=(?:[^"]*"[^"]*")*[^"]*$)', fmod)
    
    for arg in args_list:
        ## Split by equals, ignoring equals inside quotes
        p_arg = re.split(r'=(?=(?:[^"]*"[^"]*")*[^"]*$)', arg)
        
        if len(p_arg) == 1:
            if '_module' not in dict_args:
                dict_args['_module'] = p_arg[0]
        elif len(p_arg) > 1:
            key = p_arg[0]
            val_str = p_arg[1]
            
            ## If there are multiple '=' parts, rejoin the rest
            if len(p_arg) > 2:
                dict_args[key] = '='.join(p_arg[1:])
            else:
                dict_args[key] = _parse_value_string(val_str)
        
    return dict_args


def dict2fmod(in_dict: Dict[str, Any]) -> str:
    """Convert a dict of key:val pairs to a module factory string.

    Args:
      in_dict (dict): The dictionary to convert.

    Returns:
      str: A string representation of in_dict suitable for a factory CLI.
    """
    
    out_args = []
    keys = list(in_dict.keys())
    
    for key in keys:
        if key == '_module':
            out_args.append(f'{in_dict[key]}')
        elif isinstance(in_dict[key], list):
            out_args.append(f'{key}="{";".join(in_dict[key])}"')
        else:
            out_args.append(f'{key}={in_dict[key]}')
            
    return ':'.join(out_args)


def args2dict(args: List[str], dict_args: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    """Convert list of arg strings (key=val) to dict.
    
    Args:
      args (list): A list of 'key=val' strings.
      dict_args (dict, optional): A dict to append to.

    Returns:
      dict: A dictionary of the key/values.
    """
    
    if dict_args is None:
        dict_args = {}

    for arg in args:
        p_arg = arg.split('=')
        if len(p_arg) > 1:
            key = p_arg[0]
            if len(p_arg) > 2:
                val = '='.join(p_arg[1:])
            else:
                val = _parse_value_string(p_arg[1])
            dict_args[key] = val
        
    return dict_args


def dict2args(in_dict: Dict[str, Any]) -> str:
    """Convert a dict of key:val pairs to a string of quoted arguments.

    Args:
      in_dict (dict): The dictionary to convert.

    Returns:
      str: A string representation.
    """
    
    out_args = []
    for key, val in in_dict.items():
        out_args.append(f'{key}="{val}"')
    return ':'.join(out_args)


def _set_params(mc: Any, mod: Any = None, mod_name: str = None, mod_args: Optional[Dict] = None) -> None:
    """Set module parameters into mc module class."""
    
    if not hasattr(mc, 'params'):
        mc.params = {}
    
    if mod_args is None:
        mod_args = {}

    mc.params.setdefault('mod', mod)
    mc.params.setdefault('mod_name', mod_name)
    mc.params.setdefault('mod_args', {})
    
    # If kwargs not present, copy current instance dict
    if 'kwargs' not in mc.params:
        mc.params['kwargs'] = mc.__dict__.copy()


def _set_mod_params(mc: Any, mf: Any = None, mod: Any = None, mod_name: str = None) -> None:
    """Set specific module parameters based on factory lookup."""
    
    if mod is not None:
        mc.params['mod'] = mod
        
    if mod_name is not None:
        mc.params['mod_name'] = mod_name

    if mc.params.get('mod') is None:
        mc.params['mod'] = _factory_module_check(mf, mc)

    if mc.params.get('mod_name') is None:
        mc.params['mod_name'] = mc.params['mod']
        
    if not mc.params['mod_args']:
        for k, v in mc.__dict__.items():
            if k not in mc.params['kwargs']:
                mc.params['mod_args'][k] = v


def _factory_module_check(mf: Any, mc: Any) -> Optional[str]:
    """Check if factory module exists in the module class.

    Returns:
      str: Module key if found, else None.
    """
    
    for k in mf._modules.keys():
        ## Check if the module class matches the factory definition
        try:
            if isinstance(mf._modules[k]['call'](), type(mc)):
                return k
        except Exception:
            continue
    return None


## ==============================================
## Module Description Helpers
##
## modules are a dictionary with the module name
## as the key and at least a 'class' key which
## points to the class/function to call for the module
## uses <class>.__doc__ as description
## ==============================================
def get_module_short_desc(m: Dict) -> str:
    return ', '.join([f'{key}' for key in m])


def get_module_name_short_desc(m: Dict) -> str:
    return ',  '.join([
        f"{m[key].get('name', 'Unknown')} ({key})" for key in m
    ])


def get_module_long_desc(m: Dict) -> str:
    cmd = os.path.basename(sys.argv[0])
    header = f"{cmd} modules:\n% {cmd} ... <mod>:key=val:key=val...\n\n  "
    
    rows = []
    for key in m:
        name = m[key].get('name', key)
        doc = m[key]['call'].__doc__ if m[key]['call'].__doc__ else "No description available."
        rows.append(f"\033[1m{key} ({name}) \033[0m{doc}\n")
    
    return header + '\n  '.join(rows) + '\n'


def get_module_md(m: Dict) -> str:
    header = "| **Name** | **Module-Key** | **Description** |\n|---|---|---|\n"
    rows = []
    for key in m:
        name = m[key].get('name', key)
        desc = m[key].get('description', 'No description')
        rows.append(f"| {name} | {key} | {desc} |")
    return header + '\n'.join(rows)


def get_module_md_table(m: Dict) -> str:
    header = "| **Name** | **Module-Key** | **Description** |\n|---|---|---|\n"
    rows = []
    for key in m:
        name = m[key].get('name', key)
        desc = m[key].get('description', 'No description')
        rows.append(f"| {name} | {key} | {desc} |")
    return header + '\n'.join(rows)


## ==============================================
## Aliases
## ==============================================
_cudem_module_short_desc = get_module_short_desc
_cudem_module_name_short_desc = get_module_name_short_desc
_cudem_module_long_desc = get_module_long_desc
_cudem_module_md = get_module_md_table
_cudem_module_md_table = get_module_md_table


def echo_modules(module_dict: Dict, key: Any, md: bool = False) -> None:
    """Print out the existing modules from module_dict and their descriptions."""
    
    if key is None:
        if not md:
            sys.stderr.write(get_module_long_desc(module_dict))
        else:
            print(get_module_md_table(module_dict))
    else:
        selected_mod = None
        if key in module_dict:
            selected_mod = {key: module_dict[key]}
        else:
            ## Search by 'name' attribute
            for mk in module_dict:
                if key == module_dict[mk].get('name'):
                    selected_mod = {mk: module_dict[mk]}
                    break

        if selected_mod:
            sys.stderr.write(get_module_long_desc(selected_mod))
        else:
            sys.stderr.write(
                f'Invalid Module Key: {key}\n'
                f'Valid Modules: {get_module_name_short_desc(module_dict)}\n'
            )

    sys.stderr.flush()
    return 0


## ==============================================
## Default Factory Compatible Module
## ==============================================
class CUDEMModule:
    """CUDEM factory module.

    Template for a compatible factory module setup, together with 
    CUDEMFactory.
    """
    
    def __init__(self, modules: Dict = None, params: Dict = None, **kwargs):
        self.params = params if params is not None else {}
        self.modules = modules if modules is not None else {}
        self.kwargs = kwargs
        
        ## Set kwargs as attributes
        for kpam, kval in kwargs.items():
            setattr(self, kpam, kval)

        self.run()

        
    def initialize(self):
        pass

    
    def __call__(self):
        print(self.__dict__)
        self.run()

        
    def run(self):
        raise NotImplementedError(
            f'The module {self.params.get("mod")} could not be parsed by the sub factory, '
            f'available modules are: {get_module_short_desc(self.modules)}'
        )


## ==============================================
## Main Factory Class
## ==============================================
class CUDEMFactory:
    """CUDEM module factory.

    Template for a compatible module factory setup, together with 
    CUDEMModule.
    
    The _modules dictionary should contain:
    { 'key': {'name': 'surface', 'description': '...', 'call': ClassReference}, ... }
    """
    
    _factory_module = {
        '_factory': {
            'name': 'factory',
            'description': 'default factory setting',
            'fmts': [],
            'call': CUDEMModule
        }
    }
    
    ## Class-level default modules
    _modules = {
        '_factory': {
            'name': 'factory',
            'description': 'default factory setting',
            'call': CUDEMModule
        }
    }

    def __init__(
            self,
            mod: Optional[str] = None,
            mod_name: str = '_factory',
            mod_args: Optional[Dict] = None,
            **kwargs: Any
    ):
        """
        Initialize the factory default settings.
        
        Parameters:
          mod: Module string 'mod_name:arg=val'
          mod_name: Name of the module
          mod_args: Dictionary of module arguments
          kwargs: Additional arguments
        """
        
        self.mod = mod
        self.mod_name = mod_name
        self.mod_args = mod_args if mod_args is not None else {}
        self.kwargs = kwargs
        
        ## IMPORTANT: Create an instance copy of _modules so we don't modify the class variable
        self._modules = self._modules.copy()

        if self.mod is not None:
            self._parse_mod(self.mod)

        self.add_module(self._factory_module)

        
    def __str__(self):
        return f'<{self.__dict__}>'

    
    def __repr__(self):
        return f'<{self.__dict__}>'

    
    def _parse_mod(self, mod: str):
        """Parse the module string.

        Returns:
          (module-name, module-arguments)
        """
        
        opts = fmod2dict(mod, {})
        if '_module' not in opts:
            utils.echo_error_msg(f'Could not parse module string: {opts}')
            self.mod_args['modules'] = self._modules
            return None
            
        if opts['_module'] in self._modules:
            self.mod_name = opts['_module']
            self.mod_args = {i: opts[i] for i in opts if i != '_module'}
        else:
            utils.echo_error_msg(f"Invalid module name `{opts['_module']}`")            
            self.mod_args['modules'] = self._modules
            
        return self.mod_name, self.mod_args

    
    def add_module(self, type_def: Dict = None):
        """Add a module to the factory `_modules` dict."""
        
        if type_def is None:
            return
        
        for key, val in type_def.items():
            self._modules[key] = val

    def _acquire_module(self):
        """Acquire (instantiate) the module from the factory."""

        def _get_all_keys(d):
            for key, value in d.items():
                if key == 'params':
                    continue
                yield key
                if isinstance(value, dict):
                    yield from _get_all_keys(value)
        
        ## Add factory params to kwargs
        self.kwargs['params'] = self.__dict__
        
        ## Merge mod_args into kwargs, prioritizing mod_args
        keys_to_remove = []
        for k in self.kwargs:
            if k in self.mod_args:
                ## self.kwargs takes precedence or is overwritten here
                self.kwargs[k] = self.mod_args[k]
                keys_to_remove.append(k)
        
        for k in keys_to_remove:
            del self.mod_args[k]
                
        if self.mod_name is not None and self.mod_name in self._modules:
            ## Instantiate the module class
            module_class = self._modules[self.mod_name]['call']
            
            ## Call the module
            ## The module class should accept **mod_args and **kwargs for arguments
            instance = module_class(**self.mod_args, **self.kwargs)

            ## Validation check (warn if arguments were unused)
            instance_keys = list(_get_all_keys(instance.__dict__))
            for mod_arg in self.mod_args:
                if mod_arg not in instance_keys:
                    utils.echo_warning_msg(f'{mod_arg} is not a valid parameter for {self.mod_name}...')

            return instance
        else:
            return None

        
    def load_parameter_dict(self, param_dict: Dict):
        """Load a module parameter dictionary."""
        
        valid_data = False

        for ky, val in param_dict.items():
            setattr(self, ky, val)
            valid_data = True
                    
        if valid_data:
            utils.echo_msg('CUDEMFactory read successfully')
        else:
            utils.echo_warning_msg('Unable to find any valid data')
            
        return self

    
    def write_parameter_file(self, param_file: str):
        """Write a module parameter file to disk."""
        
        try:
            with open(param_file, 'w') as outfile:
                params = self.__dict__.copy()
                params.pop('_modules')
                ## Convert paths to absolute paths
                if 'mod_args' in params:
                    for key in params['mod_args']:
                        val = params['mod_args'][key]
                        if isinstance(val, str) and os.path.exists(val):
                            params['mod_args'][key] = os.path.abspath(val)

                for key in params:
                    val = params[key]
                    if isinstance(val, str) and os.path.exists(val):
                        params[key] = os.path.abspath(val)

                utils.echo_msg(params)
                json.dump(params, outfile, indent=4)                
                utils.echo_msg(f'New CUDEMFactory file written to {param_file}')
                
        except Exception as e:
            raise ValueError(
                f'CUDEMFactory: Unable to write new parameter file to {param_file}. Error: {e}'
            )

        
    def open_parameter_file(self, param_file: str):
        """Open and read a saved parameter file."""
        
        valid_data = False
        try:
            with open(param_file, 'r') as infile:
                data = json.load(infile)
        except Exception as e:
            raise ValueError(
                f'CUDEMFactory: Unable to read data from {param_file} as json. Error: {e}'
            )
        
        for ky, val in data.items():
            setattr(self, ky, val)
            valid_data = True
                
        if valid_data:
            utils.echo_msg(f'CUDEMFactory read successfully from {param_file}')
        else:
            utils.echo_warning_msg(f'Unable to find any valid data in {param_file}')


if __name__ == '__main__':
    f = CUDEMFactory()


### End
