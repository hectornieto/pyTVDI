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

from TVDI_ConfigFile_Interface import TVDI_ConfigFile_Interface
import pyTVDI as tvdi

class TVDI_IPython_Interface(TVDI_ConfigFile_Interface):
    def __init__(self):

        TVDI_ConfigFile_Interface.__init__(self)        
        
        '''Initialize input variables  with default  values'''
        self.geoflag=tvdi.cWHOLE_IMG
        self.x_pix=0
        self.y_pix=0
        self.dim_cols=0
        self.dim_rows=0
        self.moving_window=tvdi.cTILES
        self.window_size_x=0
        self.window_size_y=0
        self.dry_edge=tvdi.cTANG
        self.ts_min=tvdi.cMEAN
        self.output=tvdi.cTVDI
        self.ndvi_step=0.01
        self.ndvi_lower_limit=0.1
        self.ts_min_params=20
        self.ts_min_file=''
        self.min_ndvi=0.1
        self.max_ndvi=0.9
        self.max_fi=1.26
        self.interpolation=tvdi.cLINEAR
        self.ndvi_file='./Input/NDVI_example.tif'
        self.ts_file='./Input/LST_example.tif'
        self.CLM_file='./Input/Mask_example.tif'
        self.delta_file=''
        self.output_dir='./Output/'
        self.ndvi_mult=1
        self.ts_mult=1
        self.delta_mult=1
                   
    def TVDIWidget(self):
        '''Creates a jupyter notebook GUI for running TSEB for an image'''
        import ipywidgets as widgets
        from IPython.display import display
        # Load and save configuration buttons
        self.w_loadconfig=widgets.Button(description='Load Configuration File')
        self.w_saveconfig=widgets.Button(description='Save Configuration File')
        # Input and output images
        self.w_LST=widgets.Button(description='Browse LST Image', width=50)
        self.w_LSTtxt=widgets.Text(description='LST:',value=self.ts_file, width=500)
        self.w_ts_mult=widgets.FloatText(description='LST scale factor: ',value=self.ts_mult, width=50)
        self.w_VI=widgets.Button(description='Browse VI Image', width=50)
        self.w_VItxt=widgets.Text(description='VI:',value=self.ndvi_file, width=500)
        self.w_ndvi_mult=widgets.FloatText(description='NDVI scale factor: ',value=self.ndvi_mult, width=50)

        self.w_CLM=widgets.Button(description='Browse mask Image', width=50)
        self.w_CLMtxt=widgets.Text(description='Mask:',value=self.CLM_file, width=500)

        self.w_delta=widgets.Button(description='Browse Delta Image', width=50)
        self.w_deltatxt=widgets.Text(description='Delta:',value=str(self.max_fi), width=500)
        self.w_delta_mult=widgets.FloatText(description='Delta scale factor: ',value=self.delta_mult, width=50)
        self.w_delta.visible=False
        self.w_deltatxt.visible=False
        self.w_delta_mult.visible=False

        self.w_output=widgets.Button(description='Select Output Folder')
        self.w_outputtxt=widgets.Text(description='Output Folder :', value=self.output_dir, width=500)

        # Run PyTSEB button
        self.w_runmodel=widgets.Button(description='Run TVDI', background_color='green')
        # Create TSEB options widgets
        self.SelectModel()
        self.ROI()
        self.DryEdge()
        self.WetEdge()
        # Model Selection tabs
        tabs = widgets.Tab(children=[self.modelPage,self.ROIPage,self.EdgePage])
        tabs.set_title(0, 'VI-TS triangle output parameters')
        tabs.set_title(1, 'Region of Interest options')
        tabs.set_title(2, 'Edge calculation options')
        # Display widgets
        display(self.w_loadconfig)        
        display(widgets.VBox([widgets.HTML('Select Surface Temperature Image'),
            widgets.HBox([self.w_LST,self.w_LSTtxt,self.w_ts_mult]),
            widgets.HTML('Select Vegetation Index Image'),
            widgets.HBox([self.w_VI,self.w_VItxt,self.w_ndvi_mult]),
            widgets.HTML('Select Processing Mask Image'),
            widgets.HBox([self.w_CLM,self.w_CLMtxt]),
            widgets.HTML('Select Priestley-Taylor Coefficient'),
            widgets.HBox([self.w_delta,self.w_deltatxt,self.w_delta_mult])]))
        display(widgets.HBox([self.w_output,self.w_outputtxt]))
        display(tabs)
        
        display(self.w_saveconfig) 
        display(self.w_runmodel)
        widgets.HTML('Select Delta/Delta+gamma or type a constant value'),
        widgets.HBox([self.w_delta,self.w_deltatxt,self.w_delta_mult]),

        # Handle interactions
        self.w_LST.on_click(self._on_inputLST_clicked)    
        self.w_VI.on_click(self._on_inputVI_clicked)    
        self.w_CLM.on_click(self._on_inputCLM_clicked)  
        self.w_delta.on_click(self._on_inputDelta_clicked)  
        self.w_Tsmin.on_click(self._on_inputTsmin_clicked)  
        self.w_output.on_click(self._on_output_clicked)    
        self.w_loadconfig.on_click(self._on_loadconfig_clicked)
        self.w_saveconfig.on_click(self._on_saveconfig_clicked)
        self.w_runmodel.on_click(self._on_runmodel_clicked)

        self.isImage=True

    def _on_loadconfig_clicked(self,b):
        '''Reads a configuration file and parses its data into the GUI'''
        
        InputFile=self.GetInputFileName(title='Select Input Configuration File')
        if not InputFile:return
        configdata=self.parseInputConfig(InputFile)

        # Update the widget fields
        # I/O Parameters
        self.w_VItxt.value=str(configdata['ndvi_file']).strip('"')
        self.w_LSTtxt.value=str(configdata['ts_file']).strip('"')
        self.w_CLMtxt.value=str(configdata['CLM_file']).strip('"')
        self.w_deltatxt.value=str(configdata['delta_file']).strip('"')
        self.w_outputtxt.value=str(configdata['output_dir']).strip('"')
        self.w_ndvi_mult.value=configdata['ndvi_mult']
        self.w_ts_mult.value=configdata['ts_mult']
        self.w_delta_mult.value=configdata['delta_mult']

        # ROI parameters
        self.w_geoflag.value=int(configdata['geoflag'])
        self.w_x_pix.value=float(configdata['x_pix'])
        self.w_y_pix.value=float(configdata['y_pix'])
        self.w_dim_cols.value=float(configdata['dim_cols'])
        self.w_dim_rows.value=float(configdata ['dim_rows'])
        self.w_moving_window.value=int(configdata['moving_window'])
        self.w_window_size_x.value=float(configdata['window_size_x'])
        self.w_window_size_y.value=float(configdata['window_size_y'])

        # Algorithm Parameters
        self.w_dry_edge.value=int(configdata['dry_edge'])
        self.w_ts_min.value=int(configdata['ts_min'])
        self.w_model.value=int(configdata['output'])
        self.w_ndvi_step.value=configdata['ndvi_step']
        self.w_ndvi_lower_limit.value=configdata['ndvi_lower_limit']
        self.w_ts_min_params.value=configdata['ts_min_params']
        #self.w_ts_min_file.value=configdata['ts_min_file']
        self.w_min_ndvi.value=configdata['min_ndvi']
        self.w_max_fi.value=configdata['max_fi']
        self.w_interpolation.value=int(configdata['interpolation'])

       

    def _on_saveconfig_clicked(self,b):
        '''Opens a configuration file and writes the parameters in the GUI into the file'''
        OutputFile=self.GetOutputFileName(title='Select Output Configuration File')
        if not OutputFile: return
        try:
            fid=open(OutputFile,'w')
        except IOError:
            print('Could not write ' +OutputFile)
            return
        fid.write('# Input files\n')
        fid.write('ndvi_file='+str(self.w_VItxt.value)+'\n')
        fid.write('ts_file='+str(self.w_LSTtxt.value)+'\n')
        fid.write('CLM_file='+str(self.w_CLMtxt.value)+'\n')
        fid.write('delta_file='+str(self.w_deltatxt.value)+'\n')

        fid.write('\n# Output folder\n')
        fid.write('output_dir='+str(self.w_outputtxt.value)+'\n')

        fid.write('\n# Additional I/O Options\n')
        fid.write('ndvi_mult='+str(self.w_ndvi_mult.value)+'\n')
        fid.write('ts_mult='+str(self.w_ts_mult.value)+'\n')
        fid.write('delta_mult='+str(self.w_delta_mult.value)+'\n')
        
        fid.write('\n# Region Of Interest Options\n')
        fid.write('geoflag='+str(self.w_geoflag.value)+'\n')
        fid.write('x_pix='+str(self.w_x_pix.value)+'\n')
        fid.write('y_pix='+str(self.w_y_pix.value)+'\n')
        fid.write('dim_cols='+str(self.w_dim_cols.value)+'\n')
        fid.write('dim_rows='+str(self.w_dim_rows.value)+'\n')
        fid.write('moving_window='+str(self.w_moving_window.value)+'\n')
        fid.write('window_size_x='+str(self.w_window_size_x.value)+'\n')
        fid.write('window_size_y='+str(self.w_window_size_y.value)+'\n')

        fid.write('\n# Edge Retrieval Algorithm parameters\n')
        fid.write('dry_edge='+str(self.w_dry_edge.value)+'\n')
        fid.write('ts_min='+str(self.w_ts_min.value)+'\n')
        fid.write('output='+str(self.w_model.value)+'\n')
        fid.write('ndvi_step='+str(self.w_ndvi_step.value)+'\n')
        fid.write('ndvi_lower_limit='+str(self.w_ndvi_lower_limit.value)+'\n')
        fid.write('ts_min_params='+str(self.w_ts_min_params.value)+'\n')
        #fid.write('ts_min_file='+str(self.w_ts_min_file.value)+'\n')
        fid.write('min_ndvi='+str(self.w_min_ndvi.value)+'\n')
        fid.write('max_fi='+str(self.w_max_fi.value)+'\n')
        fid.write('interpolation='+str(self.w_interpolation.value)+'\n')

        fid.flush()
        fid.close()
        del fid
        print('Saved Configuration File')

           
    def SelectModel(self):
        ''' Widget to select the TVDI output model'''
        import ipywidgets as widgets
        self.w_model=widgets.ToggleButtons(description='Select Triangle output to run:',
            options={'Temperature-Vegetation Dryness Index':tvdi.cTVDI, 'Evaporative Fraction':tvdi.cEF},
            value=self.output,width=500)
        
        self.w_interpolation=widgets.ToggleButtons(description='Select Triangle Isolins interpolatio method:',
            options={'Linear':tvdi.cLINEAR, 'Quadratic':tvdi.cSQUARE},value=self.interpolation,width=50)
        
        self.w_min_ndvi=widgets.FloatText(description='Set minimum VI value to calculate the output',value=self.min_ndvi,width=50)
        self.w_max_ndvi=widgets.FloatText(description='Set maximum VI value to calculate the output',value=self.max_ndvi,width=50)
        self.w_max_fi=widgets.FloatText(description='Set Priestley-Taylor value for maximum Evapotranspiration',value=self.max_fi,width=50)
        self.w_max_fi.visible=False
        
        self.modelPage=widgets.VBox([self.w_model,
                    widgets.HBox([self.w_interpolation,self.w_min_ndvi,self.w_max_ndvi,self.w_max_fi])], background_color='#EEE')

        self.w_model.on_trait_change(self._on_model_change,'value')

    def _on_model_change(self,name, value):
        '''Behaviour when TVDI output model is changed'''
        if value==tvdi.cEF:
            self.w_delta.visible=True
            self.w_deltatxt.visible=True
            self.w_delta_mult.visible=True
            self.w_max_fi.visible=True
        else:
            self.w_delta.visible=False
            self.w_deltatxt.visible=False
            self.w_delta_mult.visible=False
            self.w_max_fi.visible=False

    def ROI(self):
        '''Widgets for site description parameters'''
        import ipywidgets as widgets
        self.w_geoflag=widgets.Dropdown(options={'Use Whole Image':tvdi.cWHOLE_IMG,
                    'Set ROI by pixel coordinates':tvdi.cPIXEL_COORD,
                    'Set ROI by map coordinates':tvdi.cGEOG_COORD},
                    description='Set type of ROI',value=tvdi.cWHOLE_IMG,width=500)
        self.w_x_pix=widgets.FloatText(value=0,description='Top left X coord.',width=50)
        self.w_y_pix=widgets.FloatText(value=0,description='Top left Y coord.',width=50)
        self.w_dim_cols=widgets.FloatText(value=0,description='ROI width',width=50)
        self.w_dim_rows=widgets.FloatText(value=0,description='ROI height',width=50)
        self.w_x_pix.visible=False
        self.w_y_pix.visible=False
        self.w_dim_cols.visible=False
        self.w_dim_rows.visible=False
        
        self.w_moving_window=widgets.Dropdown(options={'Do not split the ROI':tvdi.cNO_WINDOW,
                    'Split the ROI by tiles':tvdi.cTILES,
                    'Split the ROI by a moving window':tvdi.cMOVING_WINDOW},
                    description='Define the moving window size',value=tvdi.cNO_WINDOW,width=500)
        self.w_window_size_x=widgets.FloatText(value=0,description='Width',width=50)
        self.w_window_size_y=widgets.FloatText(value=0,description='Height',width=50)
        self.w_window_size_x.visible=False
        self.w_window_size_y.visible=False

        self.ROIPage=widgets.VBox([self.w_geoflag,
                                   widgets.HBox([self.w_x_pix,self.w_y_pix,self.w_dim_cols,self.w_dim_rows]),
                                    self.w_moving_window,
                                    widgets.HBox([self.w_window_size_x,self.w_window_size_y])], background_color='#EEE')

        self.w_geoflag.on_trait_change(self._on_geoflag_change,'value') 
        self.w_moving_window.on_trait_change(self._on_moving_window_change,'value') 
    
    def _on_geoflag_change(self,name, value):
        '''Behaviour when selecting type of ROI'''
        if value==tvdi.cWHOLE_IMG:
            self.w_x_pix.visible=False
            self.w_x_pix.value=0
            self.w_y_pix.visible=False
            self.w_y_pix.value=0
            self.w_dim_cols.visible=False
            self.w_dim_cols.value=0
            self.w_dim_rows.visible=False
            self.w_dim_rows.value=0
        else:
            self.w_x_pix.visible=True
            self.w_y_pix.visible=True
            self.w_dim_cols.visible=True
            self.w_dim_rows.visible=True

    def _on_moving_window_change(self,name, value):
        '''Behaviour when selecting whether to split ROI in tiles/moving windows'''
        if value==tvdi.cNO_WINDOW:
            self.w_window_size_x.visible=False
            self.w_window_size_x.value=0
            self.w_window_size_y.visible=False
            self.w_window_size_y.value=0
        else:
            self.w_window_size_x.visible=True
            self.w_window_size_y.visible=True
   
    def DryEdge(self):
        '''Widgets for site description parameters'''
        import ipywidgets as widgets
        self.w_dry_edge=widgets.Dropdown(options={'Simple':tvdi.cSIMPLE,
                    'Tang':tvdi.cTANG},
                    description='Set method for dry edge calculation',value=tvdi.cSIMPLE,width=50)
        self.w_ndvi_step=widgets.FloatText(description='Set size of VI bins:',value=self.ndvi_step,width=50)
        self.w_ndvi_lower_limit=widgets.FloatText(description='Set minimum value VI to use in Edge calculation:',
                                                  value=self.ndvi_lower_limit,width=50)
       
    def WetEdge(self):
        '''Widgets for site spectral properties'''
        import ipywidgets as widgets
        self.w_ts_min=widgets.Dropdown(options={'Mean':tvdi.cMEAN,
                    'Median':tvdi.cMEDIAN,'Variable Max. VI':tvdi.cVAR_MAX_NDVI,
                    'Constant Min. Temp.':tvdi.cCONST_TS_MIN},
                    description='Set method for wet edge calculation',value=tvdi.cMEAN,width=50)
        self.w_ts_min_params=widgets.FloatText(value=self.ts_min_params,
                           description='Number of VI bins to be used in calculation',width=50)

        self.w_Tsmin=widgets.Button(description='Browse Minimum Teperature Image')
        self.w_Tsmin.visible=False
        self.w_Tsmintxt=widgets.Text(description='T min:',value=self.ts_min_file, width=500)
        self.w_Tsmintxt.visible=False
        
        self.w_ts_min.on_trait_change(self._on_ts_min_change,'value') 

        self.EdgePage=widgets.VBox([widgets.HBox([self.w_dry_edge,self.w_ndvi_step,self.w_ndvi_lower_limit]),
                                    widgets.HBox([self.w_ts_min,self.w_ts_min_params]),
                                    widgets.HBox([ self.w_Tsmin,self.w_Tsmintxt])], background_color='#EEE')


    def _on_ts_min_change(self,name, value):
        '''Behaviour when selecting wet edge method'''
        if value==tvdi.cMEAN or tvdi.cMEDIAN:
            self.w_ts_min_params.description='Number of VI bins to be used in calculation'
            self.w_Tsmin.visible=False
            self.w_Tsmintxt.visible=False
        elif value==tvdi.cVAR_MAX_NDVI:
            self.w_ts_min_params.description='Maximum value of VI to extrapolate the dry edge'
            self.w_Tsmin.visible=False
            self.w_Tsmintxt.visible=False
        elif value==tvdi.cCONST_TS_MIN:
            self.w_ts_min_params.description='Ts value of wet edge'
            self.w_ts_min_params.visible=True
            self.w_Tsmin.visible=True
            self.w_Tsmintxt.visible=True
         
    def _on_inputLST_clicked(self,b):
        '''Behaviour when clicking the LST input file button'''
        self.w_LSTtxt.value=str(self.GetInputFileName(title='Select Surface Temperature Image'))

    def _on_inputVI_clicked(self,b):
        '''Behaviour when clicking the LST VZA file button'''
        self.w_VItxt.value=self.GetInputFileName(title='Select Vegetation Index Image')

    def _on_inputCLM_clicked(self,b):
        '''Behaviour when clicking the LAI input file button'''
        self.w_CLMtxt.value=self.GetInputFileName(title='Select Mask Image')

    def _on_inputDelta_clicked(self,b):
        '''Behaviour when clicking the Fc input file button'''
        self.w_deltatxt.value=self.GetInputFileName(title='Select Delta Image')

    def _on_inputTsmin_clicked(self,b):
        '''Behaviour when clicking the Fg input file button'''
        self.w_Tsmintxt.value=self.GetInputFileName(title='Select Minimum Temperature Image')

    def GetInputFileName(self, title='Select Input File'):
        root, askopenfilename, _ = self._setup_tkinter()
        InputFile = askopenfilename(parent = root, title=title) # show an "Open" dialog box and return the path to the selected file
        root.destroy() # Destroy the GUI
        return InputFile
    
    def _on_output_clicked(self,b):
        '''Behaviour when clicking the output file button'''
        self.w_outputtxt.value=self.GetOutputFolder()

    def GetOutputFolder(self, title='Select Output Folder'):
        root, _, askdirectory = self._setup_tkinter()
        OutputFolder = askdirectory(title=title) # show an "Open" dialog box and return the path to the selected file
        root.destroy()  # Destroy the GUI
        return OutputFolder

    def _setup_tkinter(self):
        '''Creates a Tkinter input file dialog'''
        import sys
        # Import Tkinter GUI widgets
        if sys.version_info.major==2:
            from tkFileDialog import askopenfilename, askdirectory
            import Tkinter as tk
        else:
            from tkinter.filedialog import askopenfilename, askdirectory
            import tkinter as tk
        
        # Code below is to make sure the file dialog appears above the 
        # terminal/browser
        # Based on http://stackoverflow.com/questions/3375227/how-to-give-tkinter-file-dialog-focus

        # Make a top-level instance and hide since it is ugly and big.
        root=tk.Tk()
        root.withdraw()

        # Make it almost invisible - no decorations, 0 size, top left corner.
        root.overrideredirect(True)
        root.geometry('0x0+0+0')
        
        # Show window again and lift it to top so it can get focus,
        # otherwise dialogs will end up behind the terminal.
        root.deiconify()
        root.lift()
        root.focus_force()
        
        return root, askopenfilename, askdirectory
    
    def _on_runmodel_clicked(self, b):
        # Change the colour of the button to know it is running
        self.w_runmodel.background_color='yellow'
        # Get the data from the widgets
        self.GetDataTVDIWidgets()
        #run TSEB
        tvdi.tvdi(self.io_inf,self.roi_inf,self.alg_inf)
        # Change the colour of the button to know it is running
        self.w_runmodel.background_color='green'

    def GetDataTVDIWidgets(self):
        '''Parses the parameters in the GUI to TVDI variables for running tvdi'''
        self.roi_inf = {'geoflag':self.w_geoflag.value,
                       'x_pix':self.w_x_pix.value,'y_pix':self.w_y_pix.value,
                       'dim_cols':self.w_dim_cols.value,'dim_rows':self.w_dim_rows.value, 
                       'moving_window':self.w_moving_window.value,
                       'window_size_x':self.w_window_size_x.value,'window_size_y':self.w_window_size_y.value }
        
        self.alg_inf = {'dry_edge': self.w_dry_edge.value, 
                       'ts_min': self.w_ts_min.value, 
                       'output': self.w_model.value, 
                       'dry_edge_params': [self.w_ndvi_step.value, self.w_ndvi_lower_limit.value], # [ndvi_step, ndvi_lower_limit]
                       'ts_min_params': self.w_ts_min_params.value, # see constants file for explenation
                       'ts_min_file':"", 
                       'output_params': [self.w_min_ndvi.value, self.w_max_ndvi.value, 
                                         self.w_max_fi.value,  
                                         self.w_interpolation.value] } #[min_ndvi, max_ndvi, max_fi,  interpolation] 
        #default io options
        self.io_inf = {'ndvi_file': self.w_VItxt.value,
                      'ts_file': self.w_LSTtxt.value, 
                      'CLM_file': self.w_CLMtxt.value,
                      'delta_file': self.w_deltatxt.value,
                      'output_dir': self.w_outputtxt.value,
                     'ndvi_mult': self.w_ndvi_mult.value,
                     'ts_mult': self.w_ts_mult.value,
                     'delta_mult': self.w_delta_mult.value} 
