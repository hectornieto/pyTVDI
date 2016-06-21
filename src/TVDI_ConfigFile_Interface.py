# This file is part PyTSEB, consisting of of high level pyTSEB scripting
# Copyright 2016 Hector Nieto and contributors listed in the README.md file.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import pyTVDI

class TVDI_ConfigFile_Interface():

    def __init__(self):
        
        self.input_vars=('geoflag','x_pix','y_pix','dim_cols','dim_rows',
                   'moving_window','window_size_x','window_size_y','dry_edge','ts_min','output',
                   'ndvi_step','ndvi_lower_limit','ts_min_params','ts_min_file','min_ndvi',
                   'max_ndvi', 'max_fi', 'interpolation','ndvi_file',
                   'ts_file','CLM_file','delta_file','output_dir','ndvi_mult',
                   'ts_mult','delta_mult')
        
        self.params = {}        
        
    def parseInputConfig(self,InputFile):
        ''' Parses the information contained in a configuration file into a dictionary''' 
        # look for the following variable names
        from re import match
        # create the output configuration file:
        configdata=dict()
        for var in self.input_vars:
            configdata[var]=0
        try:
            fid=open(InputFile,'r')
        except IOError:
            print('Error reading ' + InputFile + ' file')
            return configdata
        for line in fid:
            if match('\s',line):# skip empty line
                continue
            elif line[0]=='#': #skip comment line
                continue
            elif '=' in line:
                line=line.split('#')[0].rstrip(' \r\n') # Remove comments in case they exist
                field,value=line.split('=')
                for var in self.input_vars:
                    if var==field:
                        configdata[field]=value
        del fid
        return configdata
    
    def GetDataTVDI(self,configdata,isImage):
        '''Parses the parameters in a configuration file directly to TVDI variables for running tvdi'''
        self.roi_inf = {'geoflag':int(configdata['geoflag']),
                       'x_pix':float(configdata['x_pix']),'y_pix':float(configdata['y_pix']),
                       'dim_cols':int(configdata['dim_cols']),'dim_rows':int(configdata['dim_rows']), 
                       'moving_window':int(configdata['moving_window']),
                       'window_size_x':int(configdata['window_size_x']),'window_size_y':int(configdata['window_size_x']) }
        
        self.alg_inf = {'dry_edge': int(configdata['dry_edge']), 
                       'ts_min': int(configdata['ts_min']), 
                       'output': int(configdata['output']), 
                       'dry_edge_params': [float(configdata['ndvi_step']), float(configdata['ndvi_lower_limit'])], # [ndvi_step, ndvi_lower_limit]
                       'ts_min_params': float(configdata['ts_min_params']), # see constants file for explenation
                       'ts_min_file':configdata['ts_min_file'], 
                       'output_params': [float(configdata['ts_min_params']), float(configdata['max_ndvi']),
                                         float(configdata['max_fi']),
                                         int(configdata['interpolation'])] } #[min_ndvi, max_ndvi, max_fi, interpolation] 
        #default io options
        self.io_inf = {'ndvi_file': configdata['ndvi_file'],
                      'ts_file': configdata['ts_file'], 
                      'CLM_file': configdata['CLM_file'],
                      'delta_file': configdata['delta_file'],
                      'output_dir': configdata['output_dir'],
                     'ndvi_mult': float(configdata['ndvi_mult']),
                     'ts_mult': float(configdata['ts_mult']),
                     'delta_mult': float(configdata['delta_mult'])} 
    
    def RunTVDI(self, isImage):
        pyTVDI.tvdi(self.io_inf,self.roi_inf,self.alg_inf)
