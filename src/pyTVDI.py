'''
Main TVDI program
Adapted for Python by Hector Nieto 
''' 
#==============================================================================
# CONSTANTS
#==============================================================================
# methods for calcualting ts min
cMEAN = 0
cMEDIAN = 1
cVAR_MAX_NDVI = 2
cCONST_TS_MIN = 3
cDAILY_TS_MIN = 4

# methods for calculating dry edge
cSIMPLE = 0
cTANG = 1

# output variable
cTVDI = 0
cEF = 1

# interpolation
cSQUARE = 0
cLINEAR = 1

# geoflag values
cGEOG_COORD = 0
cPIXEL_COORD = 1
cWHOLE_IMG = 2

# moving window or tiles
cNO_WINDOW = 0
cTILES = 1
cMOVING_WINDOW = 2

#==============================================================================
# Setting parameters for cloud mask data
#==============================================================================
#Value of cloud mask in cloud-free areas over land
CLM_OK = 1

roi_inf = {'geoflag':cWHOLE_IMG,
               'x_pix':0,'y_pix':0,
               'dim_cols':0,'dim_rows':0, 
               'moving_window':cNO_WINDOW,
               'window_size_x':0,'window_size_y':0 }

alg_inf = {'dry_edge': cTANG, 
               'ts_min': cMEAN, 
               'output': cTVDI, 
               'dry_edge_params': [0.005, 0.6], # [ndvi_step, ndvi_lower_limit]
               'ts_min_params': 20, # see constants file for explenation
               'ts_min_file':"", 
               'output_params': [0.1, 0.9, 1.26, cLINEAR] } #[min_ndvi, max_ndvi, max_fi, use_1/delat_for_max_fi, interpolation] 
#default io options
io_inf = {'ndvi_file': '/media/hector/TOSHIBA EXT/Projects/UAV_ET/Database/Cardena/Hyper/20150429/L3/20150429_NDVI_0p7.tif',
              'ts_file': '/media/hector/TOSHIBA EXT/Projects/UAV_ET/Database/Cardena/TIR/20150429/L3/20150429_thermal_subset.tif', 
              'CLM_file': '',
              'delta_file': '',
              'output_dir': '/home/hector/',
             'ndvi_mult': 1,
             'ts_mult': 100,
             'delta_mult': 1} 
         
def tvdi(io_inf,roi_inf,alg_inf):
    '''
    Calculates TVDI/Evaporative Fraction from input files, and saves the results in a GeoTIFF file
     
    Parameters
    ----------
    io_inf : dict
        Input/Output parameters.
        ndvi_file : str
            Input Vegetation Index file
        ts_file : str
            Input Surface Temperature file
        CLM_file : str, optional
            Input Cloud Mask file
        delta_file : str, optional
            Input Delta/Delta+psicr, for calculation of Evaporative Fraction
        output_dir : str
            Output folder
        ndvi_mult : float
            Vegetation Index scale factor
        ts_mult : float
            Surface Temperature scale factor
        delta_mult : float
            Delta/Delta+psicr scale factor
    roi_inf : dict
        specifies the Region of Interest. 
        - geoflag = cGEOG_COORD if you use geographic coordinates.
        - geoflag = cPIXEL_COORD if you use ns, nl.
        - geoflag = cWHOLE_IMG if you want to use the whole picture.
        x_pix, y_pix : float
            the coordinate of the top left corner of ROI, units depend on geoflag.
        dim_rows, dim_cols : int
            the size, in pixels, of the ROI.
        moving_window : bool
            set to 1 if moving window is to be used, 0 if tiles are to be used.
        window_size_x, window_size_y : int
            size of the window/tiles, in pixels, if 0 then the whole ROI is one window.
    alg_inf : dict
        specifies the algorithm to use to calcualte the triangle
        Set dry_edge = cSIMPLE or cTANG depending on which dry edge algorithm is to be used
        Set ts_min = cMEAN, cMEDIAN, cVAR_MAX_NDVI or cCONST_TS_MIN depending on which ts min algorithm is to be used
        dry_edge_params - array containing paramenters required for the selected dry edge algorithm:
            - [ndvi step, ndvi lower limit]
        ts_min_params - array containing parameters required for the selected ts min algorithm:
            [[number of max NDVI values to use, N/A]              - cMEAN
            [number of max NDVI values to use, N/A]              - cMEDIAN
            [given value of NDVI, scaling factor for max NDVI]   - cVAR_MAX_NDVI
            [given value of TS min, N/A]]                        - cCONST_TS_MIN 
            [daily Ts min value read from file, N/A]]             - cDAILY_TS_MIN                      
        output_params - array containing paremeters required for calculating EF:
            [min NDVI, max NDVI, max FI, interpolation type]
   
    Returns
    -------
    None
    
    References
    ----------
    
    '''
    
    import gdal
    #import constants
    from os.path import isfile,basename,splitext,isdir
    from numpy import zeros,logical_or
    from os import mkdir
    

    # open files
    if not isfile(io_inf['ndvi_file']):
        print(io_inf['ndvi_file'] +' not found')
        return        
    if not isfile(io_inf['ts_file']):
        print(io_inf['ts_file'] +' not found')
        return                   
    if not isfile(io_inf['CLM_file']):
        clm_fid=-1
    if alg_inf['output'] == cEF and isfile(io_inf['delta_file']):
        delta_fid=gdal.Open(io_inf['delta_file'],gdal.GA_ReadOnly)
        delta=delta_fid.GetRasterBand(1).ReadAsArray()
        del delta_fid
    else:
        delta=0
    if not isdir(io_inf['output_dir']):
        mkdir(io_inf['output_dir'])

    ndvi_fid=gdal.Open(io_inf['ndvi_file'],gdal.GA_ReadOnly)
    ts_fid=gdal.Open(io_inf['ts_file'],gdal.GA_ReadOnly)
    clm_fid=gdal.Open(io_inf['CLM_file'],gdal.GA_ReadOnly)
    # Get Geotransform
    geo=ts_fid.GetGeoTransform()
    prj=ts_fid.GetProjection()    
    #Read and process data
    #Read the ROI data and get rid of bad pixels
    if roi_inf['geoflag']==cWHOLE_IMG:
        ndvi = ndvi_fid.GetRasterBand(1).ReadAsArray()
        ts = ts_fid.GetRasterBand(1).ReadAsArray()
        if clm_fid:
            cld_mask = clm_fid.GetRasterBand(1).ReadAsArray()
        else:
            cld_mask=zeros(ndvi.shape)+CLM_OK
            
        dims=[0,0,ndvi.shape[0]-1,ndvi.shape[1]-1]
        
    else:
        dims=[roi_inf['y_pix'],roi_inf['x_pix'],
              roi_inf['y_pix']+roi_inf['dim_rows'],roi_inf['x_pix']+roi_inf['dim_cols']]
        ndvi = ndvi_fid.GetRasterBand(1).ReadAsArray()[dims[0]:dims[2],dims[1]:dims[3]]
        ts = ts_fid.GetRasterBand(1).ReadAsArray()[dims[0]:dims[2],dims[1]:dims[3]]
        if clm_fid:
            cld_mask = clm_fid.GetRasterBand(1).ReadAsArray()[dims[0]:dims[2],dims[1]:dims[3]]
        else:
            cld_mask=zeros(ndvi.shape)+CLM_OK
        if alg_inf['output'] == cEF and isfile(io_inf['delta_file']):
            delta=delta[dims[0]:dims[2],dims[1]:dims[3]]

    # Closing NDVI and LST files
    del ndvi_fid,ts_fid
    if clm_fid: del clm_fid

    ndvi = ndvi/float(io_inf['ndvi_mult'])
    ts = ts/float(io_inf['ts_mult'])
    
    #; bad pixel if ndvi being lower then the lower limit set in algorithm options or higher then 1
    # or temperature is negative
    ndvi[ndvi < alg_inf['dry_edge_params'][1]] = 0
    ts[ndvi < alg_inf['dry_edge_params'][1]] = 0
    ndvi[ts<0]=0
    ts[ts<0]=0
        
    ndvi[cld_mask != CLM_OK] = 0.0
    ts[cld_mask != CLM_OK] = 0.0
   
    # Create a text file where the edge line equations will be saved
    filename=splitext(basename(io_inf['ndvi_file']))[0]
    tvdi_out_file = io_inf['output_dir']+'/'+filename + "_TVDI" 
    lun=open(tvdi_out_file+'_line_equations.txt', 'w')
    string='image \t dry_edge_intercept \t dry_edge_slope \t wet_edge_intercept\t wet_edge_slope \t r^2 \t dry_edge_points \t total_points \t points_per_bin \t ndvi_range_ratio\n'
    lun.write(string)
    lun.flush()
    lun.close()
    # create output file name
    

    # Calculate TVDI using the method specified
    tvdi=triangle_window(ndvi, ts, delta,dims, roi_inf, alg_inf, tvdi_out_file)

    # get rid of bad or cloudy pixels in the output
    tvdi[logical_or(ndvi < alg_inf['dry_edge_params'][1],ts < 0)] = -1
    tvdi[cld_mask != CLM_OK] = -1
    
    #write GDAL tvdi file
    #set output name
    out_name = tvdi_out_file+  '.tif'
    #cut off the black "border" of width of half the window size from the final image
    if roi_inf['moving_window'] == cMOVING_WINDOW:
        tvdi = tvdi[roi_inf['window_size_y']/2:dims[2]-dims[0]-roi_inf['window_size_y']/2, roi_inf['window_size_x']/2:dims[3]-dims[1]-roi_inf['window_size_x']/2]

    # Create the output GeoTIFF file
    driver = gdal.GetDriverByName("GTiff")
    src = driver.Create(out_name, tvdi.shape[1], tvdi.shape[0], 1, gdal.GDT_Float32)
    band=src.GetRasterBand(1)
    band.WriteArray(tvdi)
    band.SetNoDataValue(-1)
    src.SetGeoTransform(geo)
    src.SetProjection(prj)    

    del band
    del src
    print('Done : TVDI saved in '+out_name)    

def triangle_window(ndvi, ts, delta, dims, roi_inf, alg_inf, tvdi_out_file): 
    ''' Calculates the TVDI/EF from numpy arrays and returns the TVDI as another array.
    It also writes an ASCII file with the TVDI statistics
    
    Parameters
    ----------
    ndvi : array_like
        Vegetation Index array.
    ts : array_like
        Surface Temperature array.
    delta : array_like
        Delta/Deta+psicr array
    dims : list
        Dimensions of the Region of Interest: [row min, column min, row max, column max]
    roi_inf : dict
        specifies the Region of Interest. 
        - geoflag = cGEOG_COORD if you use geographic coordinates.
        - geoflag = cPIXEL_COORD if you use ns, nl.
        - geoflag = cWHOLE_IMG if you want to use the whole picture.
        x_pix, y_pix : float
            the coordinate of the top left corner of ROI, units depend on geoflag.
        dim_rows, dim_cols : int
            the size, in pixels, of the ROI.
        moving_window : bool
            set to 1 if moving window is to be used, 0 if tiles are to be used.
        window_size_x, window_size_y : int
            size of the window/tiles, in pixels, if 0 then the whole ROI is one window.
    alg_inf : dict
        specifies the algorithm to use to calcualte the triangle
        Set dry_edge = cSIMPLE or cTANG depending on which dry edge algorithm is to be used
        Set ts_min = cMEAN, cMEDIAN, cVAR_MAX_NDVI or cCONST_TS_MIN depending on which ts min algorithm is to be used
        dry_edge_params - array containing paramenters required for the selected dry edge algorithm:
            - [ndvi step, ndvi lower limit]
        ts_min_params - array containing parameters required for the selected ts min algorithm:
            [[number of max NDVI values to use, N/A]              - cMEAN
            [number of max NDVI values to use, N/A]              - cMEDIAN
            [given value of NDVI, scaling factor for max NDVI]   - cVAR_MAX_NDVI
            [given value of TS min, N/A]]                        - cCONST_TS_MIN 
            [daily Ts min value read from file, N/A]]             - cDAILY_TS_MIN                      
        output_params - array containing paremeters required for calculating EF:
            [min NDVI, max NDVI, max FI, use 1/delta for max FI, interpolation type]
    tvdi_out_file : str
        Name of the the output ASCII file

    Returns
    -------
    output : array_like
        TVDI  or Evaporative Fraction array   
    
    '''

    from numpy import zeros,amin
    from re import search
    from os.path import dirname
    # Creat the output array
    output=zeros(ndvi.shape)
  
    ts_min = 0

    #Read algorithm options for alg_inf struct
    dry_edge_method = alg_inf['dry_edge']
    ts_min_method = alg_inf['ts_min']
    ts_min_params =alg_inf['ts_min_params']
    ndvi_step = alg_inf['dry_edge_params'][0]
    ndvi_lower_limit = alg_inf['dry_edge_params'][1]
    output_var = alg_inf['output']
    
    #EF options
    ndvi_min = alg_inf['output_params'][0]
    ndvi_max = alg_inf['output_params'][1]
    fi_max = alg_inf['output_params'][2]
    
    #use square or linear interpolation
    if alg_inf['output_params'][3] == cSQUARE: 
        power = 2
    else: power = 1
  
    #Initialize counting variables for loops in sub-images
  
    row = dims[0]
    col = dims[1]
    
    # Use the whole image if window_size_y and window_size_y are set to 0
    if roi_inf['window_size_y']==0:
        roi_inf['window_size_y']=ndvi.shape[0]
    if roi_inf['window_size_x']==0:
        roi_inf['window_size_x']=ndvi.shape[1]
    #Vertical loop in sub-images
    while row < dims[2]:
        
        #Horizontal loop in sub-images
        while col < dims[3]:
            
            window_size_y = amin([roi_inf['window_size_y'], dims[2]-row])
            window_size_x = amin([roi_inf['window_size_x'], dims[3]-col])
            
            #Select the window data
            window_col = col - dims[1]
            window_row = row - dims[0]
            print('col: '+str(window_col)+'   row: ' + str(window_row))
            ndvi_window = ndvi[window_row:window_row+window_size_y, window_col:window_col+window_size_x]
            ts_window = ts[window_row:window_row+window_size_y, window_col:window_col+window_size_x]
            output_window = zeros(ndvi_window.shape)
            if delta:
                delta_window = delta[window_row:window_row+window_size_y, window_col:window_col+window_size_x]
       

#==============================================================================
#             Calculate the dry edge and ts_min and TVDI
#==============================================================================
            #set plot output file name
            out_name = tvdi_out_file + '_cols_'+str(col)+'_'+str(col+window_size_x-1)+'_rows_'+str(row)+'_'+str(row+window_size_y-1)+'.png'
        
            if dry_edge_method == cSIMPLE:
                [lin_fit, ts_min, fit_stats]=calc_triangle_simple(ndvi_window, 
                    ts_window, ndvi_step, ndvi_lower_limit, ts_min_method, ts_min_params, plot_out_file=out_name)
            if dry_edge_method == cTANG:
                [lin_fit, ts_min, fit_stats]=calc_triangle_tang(ndvi_window, 
                    ts_window, ndvi_step, ndvi_lower_limit, ts_min_method, ts_min_params, plot_out_file=out_name)

            #if dry edge couldn't be found just output error value
            if lin_fit[0] == 0 and lin_fit[1] == 0:
                print('Could not calcualte dry edge for ' + out_name)
            
                output_window = zeros(ndvi_window.shape)
                output_window = -1

            #otherwise calcualte the required output
            else:
  
                #write edge line equations to a text file
  
                lun=open(tvdi_out_file+'_line_equations.txt', 'a')

                string=out_name+ '\t' +str(lin_fit[1])+ '\t' +str(lin_fit[0]) + '\t' +str(ts_min)+ '\t' +'0'+ '\t' +str(fit_stats['r'])+ '\t' +str(fit_stats['dry_edge_points'])+ '\t'+ str(fit_stats['total_points'])+ '\t'+ str(fit_stats['points_per_bin'])+ '\t'+ str(fit_stats['ndvi_range_ratio']) +'\n'
                lun.write(string)
                lun.flush()
                lun.close()
    
                # if the output variable is to be TVDI then calcualte it
                if output_var == cTVDI:
                    tvdi_window = (ts_window -ts_min)/(lin_fit[1]+lin_fit[0]*ndvi_window - ts_min)
                    tvdi_window[tvdi_window > 1.0] = 1.0
                    output_window = tvdi_window
                
                # else calculate FI or EF
                else:
                    fi_min = fi_max * ((ndvi_window - ndvi_min)/(ndvi_max - ndvi_min))**power
                    ts_max = lin_fit[1]+lin_fit[0]*ndvi_window
                    fi_window = (ts_max-ts_window)/(ts_max - ts_min)*(fi_max - fi_min) + fi_min
                    
                    #if delta file is present output EF, ohterwise output fi
                    if delta: 
                        output_window = fi_window  * delta_window
                    else:
                        output_window = fi_window
            #incorporate the tvdi_window into the appropriate place in the large tvdi matrix or just use the central pixel if moving window is used
            if roi_inf['moving_window'] == cMOVING_WINDOW:
                output[row+window_size_y/2+1-dims[0], 
                       col+window_size_x/2+1-dims[1]] = output_window[window_size_y/2+1, 
                                                                window_size_x/2+1]
                col = col+1;
            else:
                output[row-dims[0]:row-dims[0]+window_size_y, 
                       col-dims[1]:col-dims[1]+window_size_x] = output_window
                col = col+roi_inf['window_size_x']
    
        col = dims[0]
        # increase row by 1 if moving window otherwise by window size
        if roi_inf['moving_window'] == cMOVING_WINDOW:
            row = row+1;
        else:
            row = row+roi_inf['window_size_y']
    
    return output


def calc_triangle_simple(ndvi, ts, ndvi_step, ndvi_lower_limit, ts_min_method, 
                         ts_min_params, plot_out_file=None):
    ''' Estimates the dry and wet edges of the LST-VI triangle using the 
    simple method for the dry edge.
    
    Parameters
    ----------
    ndvi, ts, ndvi_step, ndvi_lower_limit, ts_min_method, 
                         ts_min_params, plot_out_dir, plot_out_file
    '''    
        
    from numpy import amax, amin, floor, zeros, isfinite, mean, where,arange,size, logical_and
    from scipy.stats import linregress as linfit

    # Dry edge value plot
    max_ndvi = amax(ndvi)
    range_ndvi = max_ndvi - ndvi_lower_limit

    lin_fit=[0,0]
    ts_min=0
    fit_stats=dict()
    
    if range_ndvi < 2*ndvi_step:
        print('Range in VI sample is too small, skipping Ts/VI space')

        return lin_fit,ts_min,fit_stats

    steps = floor(range_ndvi/ndvi_step)
    ndvi_val_old = ndvi_lower_limit
    ts_max_arr = zeros(steps)
    ts_min_arr = zeros(steps)
    ndvi_arr = zeros(steps)
  
    #fill into one dimensional array
    for j in arange(steps):
        ndvi_val_new = ndvi_step + ndvi_val_old
        ndvi_ind=where(logical_and(ndvi >= ndvi_val_old,ndvi < ndvi_val_new))
       
        
        #need at least two points in each bin
        if size(ndvi_ind) >=2:
            ts_max_arr[j] = amax(ts[ndvi_ind])
            ts_min_arr[j] = amin(ts[ndvi_ind])
        else:
            ts_max_arr[j] = 0
            ts_min_arr[j] = 0

        ndvi_arr[j] = ndvi_val_new
        ndvi_val_old = ndvi_val_new

#==============================================================================
#     Making arrays where points where ts_max is identically zero are removed
#==============================================================================
    zeroindex = where(ts_max_arr != 0)
    if size(zeroindex)==0: return lin_fit,ts_min,fit_stats
    ts_nozero = ts_max_arr[zeroindex]
    ndvi_nozero = ndvi_arr[zeroindex]
  
 
#==============================================================================
# Cutting away points to the left of the maximum LST value and points lying lower
# than the mean minimum temperature
#==============================================================================
    ts_first = ts_nozero
    ndvi_first = ndvi_nozero
    if size(ts_first) < 2: return lin_fit,ts_min,fit_stats
  
    #determine temporary Ts_min
    ts_min = mean(ts_min_arr[ts_min_arr > 0])
  
    #Determine size of reduced arrays
    array_size=size(ts_first)
  
    #Finding maximum LST cut-off
    max_plot_ind=where(ts_first == amax(ts_first))
    max_plot_ind = amin(max_plot_ind)

    # Filling arrays after maximum LST value cut-off
    ts_second = ts_first[max_plot_ind:array_size]
    ndvi_second = ndvi_first[max_plot_ind:array_size]
  
    #Cutting off points lower than the mean minimum temperature
    ts_min_index=where(ts_second > ts_min)
    if size(ts_min_index) > 0:
        ts_third=ts_second[ts_min_index]
        ndvi_third = ndvi_second[ts_min_index]
    else:        
        ndvi_third = ndvi_second
        ts_third=ts_second
  

#==============================================================================
#     Making linear  fit
#==============================================================================
    lin_fit= linfit(ndvi_third,ts_third)                                      
    if not isfinite(lin_fit[0]):
        lin_fit = [0, 0]
        return lin_fit,ts_min,fit_stats
 
    fit_stats = calc_fit_stats(lin_fit, ts_third, ndvi_third, ts, ndvi, steps, range_ndvi)                                 
    
    # calculate ts_min
    ts_min = calc_ts_min(ts_min_method, ts_min_arr, ts_min_params, lin_fit, max_ndvi)
 
    if plot_out_file:   
        # call plot routine . Create Scatter plot
        yrange = [amin(ts[ts > 0]), amax(ts)]
        plot_pro(plot_out_file, ndvi, ts, ndvi_third, ts_third, ts_min, lin_fit, yrange)
  
    return lin_fit, ts_min, fit_stats

def calc_triangle_tang(ndvi, ts, ndvi_step, ndvi_lower_limit, ts_min_method, 
                       ts_min_params, plot_out_file=None):

#==============================================================================
# for algorithm details see Tang et al., Remote Sensing of Environment 114 (2010)
# "An application of the Ts-VI triangle method with enhaced edges determination for evapotranspiration estimation ..."
#==============================================================================

    #import constants
    from numpy import amax,amin, floor,zeros, nan, isfinite, mean ,std, round, sqrt, where, arange,size, logical_and
    from scipy.stats import linregress as linfit
    
    lin_fit=[0,0]
    ts_min=0
    fit_stats=dict()
    
    max_ndvi = amax(ndvi)
    range_ndvi = max_ndvi - ndvi_lower_limit
    if range_ndvi < 2*ndvi_step:
        print('Range in VI sample is too small, skipping file')
        return lin_fit,ts_min,fit_stats

    intervals = floor(range_ndvi/ndvi_step)
    subintervals = 5
    ts_max_arr = zeros(intervals)
    ts_min_arr = zeros(intervals)
    
    ndvi_arr = zeros(intervals)
    ndvi_val_old = ndvi_lower_limit
  
  
    #for steps descriptions see the referenced paper, page 543
    #step (i)
    for j in arange(intervals):
        subint_max_ts = zeros(subintervals)
        for i in arange(subintervals):
            ndvi_val_new = float(ndvi_step/subintervals) + ndvi_val_old
            
            #step (ii)
            #need at least two points in each bin
            ndvi_ind=where(logical_and(ndvi >= ndvi_val_old,ndvi < ndvi_val_new))
            if size(ndvi_ind)>2:
                subint_max_ts[i]=amax(ts[ndvi_ind])
            else:
                subint_max_ts[i] = nan

            ndvi_val_old = ndvi_val_new
    
        #step (iii)
        not_NaN=where(isfinite(subint_max_ts))
        if size(not_NaN)==0:
            mean_max_ts = nan
            std_max_ts = nan
        else:
            subint_max_ts = subint_max_ts[not_NaN]
            mean_max_ts = mean(subint_max_ts)
            std_max_ts = std(subint_max_ts)

        while True:

            if not isfinite(std_max_ts): break
        
            #step (iv)
            num_el = size(subint_max_ts)
            subint_max_ts = subint_max_ts[subint_max_ts >= mean_max_ts - std_max_ts]
            if size(subint_max_ts) == num_el: break
            
            #step (v)
            mean_max_ts = mean(subint_max_ts)
            std_max_ts = std(subint_max_ts)
            
            #step (vi)
            if size(subint_max_ts) <= round(subintervals/2.0) or std_max_ts <= 4.0: break
    
        #step (vii)
        if isfinite(mean_max_ts): ts_max_arr[j] = mean_max_ts #otherwise ts_max_arr[j] stays as 0
        ndvi_arr[j] = ndvi_lower_limit + j*ndvi_step
        ndvi_val_old = ndvi_lower_limit + (j+1.0)*ndvi_step
        
        #finding ts_min_arr is not in the Tang algorithm
        #only do this if ts_min is set using the mean or median methods
        if ts_min_method == cMEAN or ts_min_method == cMEDIAN:
            ndvi_ind=where((ndvi >= ndvi_val_old-ndvi_step) & (ndvi < ndvi_val_old))
            if size(ndvi_ind)>0:
                ts_min_arr[j] = amin(ts[ndvi_ind])

  
    #the following is not part of the original Tang algorithm
    #remove all the invalid values before fitting the line
    valid_pix = where(ts_max_arr > 0)
    # need at least two pixels to fit a straight line
    if size(valid_pix) <= 1: return lin_fit,ts_min,fit_stats
    ts_max_arr = ts_max_arr[valid_pix]
    ndvi_arr = ndvi_arr[valid_pix]
    
    #Finding maximum LST cut-off
    max_plot_ind=where(ts_max_arr == amax(ts_max_arr))
    max_plot_ind = amin(max_plot_ind)                                        
  
    #Filling arrays after maximum LST value cut-off                     
    ts_max_arr = ts_max_arr[max_plot_ind:]                     
    ndvi_arr = ndvi_arr[max_plot_ind:] 
    
    #need at least to pixels to fit a straight line
    valid_pix = where(ts_max_arr > 0)
    if size(valid_pix) <= 1 : return lin_fit,ts_min,fit_stats
    
    while True:

        #step (viii)
        lin_fit= linfit(ndvi_arr,ts_max_arr)
        fitted_ts = lin_fit[1] + lin_fit[0]*ndvi_arr
        rmse = sqrt(sum((ts_max_arr - fitted_ts)**2)/size(ndvi_arr))
    
        #step (ix)
        good_points = where(abs(fitted_ts-ts_max_arr) <= 2*rmse)
        if size(good_points) == size(ndvi_arr) : break
        ndvi_arr = ndvi_arr[good_points]
        ts_max_arr = ts_max_arr[good_points]
        if size(ndvi_arr) < 5: break
    
    #step (x)
    lin_fit= linfit(ndvi_arr,ts_max_arr)

    fit_stats = calc_fit_stats(lin_fit, ts_max_arr, ndvi_arr, ts, ndvi, intervals, range_ndvi)   

    #calculate ts_min
    ts_min = calc_ts_min(ts_min_method, ts_min_arr, ts_min_params, lin_fit, max_ndvi)
  
    #set output name
    if plot_out_file: 
        #call plot routine . Create Scatter plot
        ts_nonzero = where(ts > 0)
        yrange = [amin(ts[ts_nonzero]), amax(ts)]
        plot_pro(plot_out_file, ndvi, ts, ndvi_arr, ts_max_arr, ts_min, lin_fit, yrange)

    return lin_fit, ts_min, fit_stats

def calc_ts_min(ts_min_method, ts_min_arr, ts_min_params, lin_fit, max_ndvi):

    #import constants
    from numpy import where, mean, median, amin, size
  
    a = lin_fit[1]
    b = lin_fit[0]

    if ts_min_method == cVAR_MAX_NDVI:
        return a + b*amin([ts_min_params, max_ndvi])
    if ts_min_method == cCONST_TS_MIN :
        return ts_min_params
    if ts_min_method == cDAILY_TS_MIN : 
        return ts_min_params
        
    if ts_min_method == cMEAN or ts_min_method == cMEDIAN:
        #take the mean or median of the min Ts values of specified number of largest NDVIs
        point_num = int(ts_min_params)
        temp = where(ts_min_arr > 0)[0]
        temp_size=size(temp)
        if temp_size<=0: return 0
        if temp_size > point_num : temp = temp[temp_size-point_num:temp_size]
        if ts_min_method == cMEAN : return mean(ts_min_arr[temp])
        if ts_min_method == cMEDIAN : return median(ts_min_arr[temp])


def calc_fit_stats(lin_fit, ts_select, ndvi_select, ts_all, ndvi_all, bins, range_ndvi):
    
    from numpy import where, size, amax, amin

    fit_stats=dict()
    fit_stats['r'] =float(lin_fit[2])
    fit_stats['dry_edge_points'] = float(size(ts_select))
    fit_stats['total_points'] = float(size(where((ndvi_all != 0) & (ts_all != 0))))
    fit_stats['points_per_bin'] = float(fit_stats['total_points']/bins)
    fit_stats['ndvi_range_ratio'] = float((amax(ndvi_select)-amin(ndvi_select))/range_ndvi)
  
    print('Linear fit parameters, y=A+Bx')                           
    print('A= '+str(lin_fit[1])+ ' B= '+str(lin_fit[0])+ ' r= '+str(fit_stats['r'])
        + ' dep= '+str(fit_stats['dry_edge_points'])+ ' tp= ' +str(fit_stats['total_points'])
        + ' ppb= '+str(fit_stats['points_per_bin'])+ ' nrr= '+str(fit_stats['ndvi_range_ratio']))
  
    return fit_stats

def plot_pro(out_name, x1array, y1array, x2array, y2array, ts_min, linefit, yrange):

    import matplotlib.pylab as plt
    from numpy import amin, arange,zeros, size
    
    # Turn interactive plotting off
    plt.ioff()
    fig=plt.figure()
    
    #set up the plot and do a scatterplot of x1array, y1array
    # do scatterplotts of the other arrays if they are present
    plt.plot(x1array,y1array, 'o',markersize=1, alpha=0.5, markeredgecolor='black', color='white')
    if size(x2array) >0 and size(y2array) > 0: plt.plot(x2array, y2array,'^', markersize=6, color='red')
    plt.xlabel('NDVI')
    plt.ylabel('TS')
    plt.title('Triangle')
    plt.ylim(yrange)
    
    # plot Dry edge
    x_dummy=arange(20)/19.0
    # scale the line to go through the whole range of x
    x_min=amin(x1array)
    x_dummy = x_dummy*(1. - x_min) + x_min
    y = linefit[1] + x_dummy*linefit[0]
    plt.plot(x_dummy,y,'k-', lw=2)
    if linefit[0]>0:
        string='y='+str(linefit[1])+ '+' +str(linefit[0])+'*x'
    else:
        string='y='+str(linefit[1])+ str(linefit[0])+'*x'
    label_pos_x=0.2
    label_pos_y=yrange[0]+0.9*(yrange[1]-yrange[0])
    
    
    # plot ts_min
    x_dummy=[x_min, 1]
    if ts_min > 0:
        plt.plot(x_dummy, zeros(2)+ts_min,'k-', lw=2)
        string=string+'\n TS min = '+str(ts_min)
    
    plt.text(label_pos_x,label_pos_y,string, color='red')
    # Save the file
    plt.savefig(out_name)
    plt.close(fig)

