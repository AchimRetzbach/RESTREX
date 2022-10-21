# -*- coding: utf-8 -*-
"""
Created on Wed Dec  8 11:27:23 2021

@author: Achim Retzbach

@version: 1.3.8
"""


import numpy as np
from os import listdir
from os import remove
from shutil import copyfile
import matplotlib.pyplot as plt
from tabulate import tabulate

plt.rcParams['axes.titley'] = 1.03
plt.rcParams['axes.formatter.limits'] = (-3,4)
plt.rcParams['figure.figsize'] = plt.rcParamsDefault['figure.figsize']
#plt.rcParams["figure.figsize"] = (12,8)

#Set variables / fields:
output_location='Desktop/RESTREX_Starterkit/Output'
current_location='Desktop/RESTREX_Starterkit/RESTREX_1.3.8.py'
kr83_kr84_norm=0.201750
kr83_kr86_norm=0.664740
rb85_rb87_norm=2.592310
sr88_sr86_norm=8.375209 #krabbenhoeft
sr88_sr86_srm987=8.37861 #error 325, nist 2007
sr87_sr86_srm987=0.71034 #error 26, nist 2007
sr84_sr86_srm987=0.05655 #error 14, nist 2007
m84_sr= 83.9134191 #error 13
m86_sr= 85.9092606 #error 12
m87_sr= 86.9088775 #error 12
m88_sr= 87.9056125 #error 12
m85_rb= 84.9117897379 #error 54
m87_rb= 86.9691805310 #error 60
m83_kr= 82.91412716 #error 32
m84_kr= 83.9114977282 #error 44
m86_kr= 85.9106106269 #error 41


def ripPath(filepath):
    '''
    Takes a string that may be a path of a file or simply the destination of a
    directory, separated by /. Scans for the last / of the string and returns
    the rest after it. Useful for just getting a directory/file name if its
    complete path is given.
    Returns the argument itsef if no / is found.

    Parameters
    ----------
    filepath : string
        Path-like string from which to rip the end off.

    Returns
    -------
    str
        Everything that comes after the last / of filepath.

    '''
    for i in range(0,len(filepath)-1):
        if filepath[len(filepath)-1-i]=='/':
            return filepath[len(filepath)-i:len(filepath)]
    return filepath


def std(arr):
    '''
    Returns the empiric variance of an input array as float.
    Note the difference to numpy.std, which gives the population variance!!!

    Parameters
    ----------
    arr : int/float-array
        The array to calculate the empiric variance of.

    Returns
    -------
    s : float
        The input array's empiric variance.'

    '''
    s=0
    for k in range(0,len(arr)):
        s=s+(arr[k]-np.mean(arr))**2
        if k==len(arr)-1:
            s=np.sqrt(s/(len(arr)-1))
    return s
        

def outlierMask(input_array,sigma_clipping=2):
    '''
    Function to return boolean mask array indicating if the corresponding value
    is an outlier.

    Parameters
    ----------
    input_array : string/float/int-array
        The array to check for outliers, 1-dimensional! If the elements are strings,
        a float-conversion is attempted!
    sigma: int/float, optional
        The distance from array mean in standard deviations after which a data
        point is considered an outlier. The default is 2.

    Returns
    -------
    mask : boolean-array
        Array of the same length as input_array indicating if the corresponding
        element in input_array is an outlier (True) or not (False).

    '''
    
    mean=np.mean(input_array)
    stdev=std(input_array)
    mask=[]
    for r in input_array:
        if np.abs(r-mean)>=sigma_clipping*stdev:
            mask.append(True)
        else:
            mask.append(False)
    return mask


def outlierCheck(input_array,return_converted=False,debug=False):
    '''
    Function to check if Neptune outliers match with own calculation and, return
    if specified.

    Parameters
    ----------
    input_array : string/float/int-array
        The array to check for outliers, 1-dimensional! If the elements are strings,
        a float-conversion is attempted!
    return_converted: boolean
        Whether to return the input_array clipped of 'X'-prefixes and converted
        to floats as second output. The default is False.
    debug : boolean, optional
        Whether to show additional debug information. Mostly only relevant if
        check_outliers is also True. The default is False.

    Returns
    -------
    error: boolean
        Error code, True if Neptune outliers didn't match own calculation somewhere.
    converted: float-array of same shape as input_array, optional
        If return_converted is True, the input_array clipped of 'X'-prefixes and
        converted to floats is returned as second output.

    '''
    if debug==True:
        print('--Outlier calculation started')
    
    #Convert array to float and calculate neptune outliers, if specified
    outliers_neptune=[]
    converted=[]
    for r in input_array:
        if r[0]=='X':
            r=r[1:]
            outliers_neptune.append(float(r))
        converted.append(float(r))
            
    #Calculate Python outliers
    outliers_python=[]
    mean=np.mean(converted)
    stdev=std(converted)
    for r in converted:
        if np.abs(r-mean)>=2*stdev:
            outliers_python.append(r)
    
    #Check if outliers match
    if debug==True:
        print('---- Neptune outliers: ',outliers_neptune)
        print('---- mean: ',mean,' std: ',std)
        print('---- Python outliers: ',outliers_python)
    error=False
    for o in outliers_neptune:
        if (o in outliers_python)==False:
            print('=================================')
            print('CAUTION: OUTLIERS DO NOT MATCH!!!')
            print('Not found: ', o, ' (N) in... ')
            print(outliers_python)
            print('=================================')
            error=True    
    for o in outliers_python:
        if (o in outliers_neptune)==False:
            print('=================================')
            print('CAUTION: OUTLIERS DO NOT MATCH!!!')
            print('Not found: ', o, ' (P) in... ')
            print(outliers_neptune)
            print('=================================')
            error=True
    if debug=='True':
        print('---- Error code: ',error)
            
    #Ausgaben
    if return_converted==False:
        return error
    if return_converted==True:
        return error,converted
    
    
def linearIsotopeCorrection(isotope_of_interest,correction_isotope,correction_ratio,debug=False):
    '''Function to apply a linear isotope interference correction to the isotope_of_interest array.
    Needs an array of the isotope used for correction of equal length as well as the isotope ratio
    as inputs.
    Outputs the corrected isotope data array.'''
    output=[]
    for i in range(0,len(isotope_of_interest)):
        output.append(isotope_of_interest[i]-correction_isotope[i]/correction_ratio)
        #output.append(isotope_of_interest[i]-np.mean(correction_isotope)/correction_ratio)
        if debug==True and i==0:
            print('Linear isotope correction (first element): ')
            print(isotope_of_interest[0], ' - ', correction_isotope[0]/correction_ratio, ' = ', output[0])
    return output


def thermoMassFractionation(measured_ratios,correction_ratios,normalizing_factor,measured_masses,normalizing_masses):
    '''
    Function to apply exponential law mass fractionation correction to isotope ratios.
    Should do the same as massFractionationCorrection(), but the formula is taken from Thermo AN 30158.

    Parameters
    ----------
    measured_ratios : float-array
        Ratios to correct.
    correction_ratios : float-array
        Ratios of known true value to correct with.
    normalizing_factor : TYPE
        Known true value of correction_ratios
    measured_masses : two-element array of ints/floats
        Masses of the isotopes of measured_ratios, numerator first, denominator second.
    normalizing_masses : two-element array of ints/floats
        Masses of the isotopes of correction_ratios, numerator first, denominator second.

    Returns
    -------
    corrected_ratios : float-array
        Fractionation-corrected alteration of measured_ratios

    '''
    corrected_ratios=[]
    for i in range(0,len(measured_ratios)):
        a=np.log(measured_masses[0]/measured_masses[1])/np.log(normalizing_masses[0]/normalizing_masses[1])
        corrected_ratios.append(measured_ratios[i]*np.power(normalizing_factor/correction_ratios[i],a))
    return corrected_ratios  


def clippedMeanError(input_array,mask=[],sigma_clipping=2):
    '''
    Function to calculate the mean and error if an array with outlier mask is given.

    Parameters
    ----------
    input_array : float/int-array, 1-dimensional!
        The array to clip and perform calculations with.
    mask : boolean array, same shape as arr, optional
        The mask indicating if the corresponding value in arr is an outlier and should
        be neglected (True) or if it should be respected (False). If not specified,
        the array is being sigma-clipped anew using outlierMask().
        The default is [].
    sigma_clipping : int/float, optional
        In case no mask is specified, this value specifies which sigma_clipping
        is being given over to outlierMask().

    Returns
    -------
    mean : float
        The mean of the masked array
    error : float
        The error of the masked array.

    '''
    clipped_array=[]
    if mask!=[]:
        if len(input_array)!=len(mask):
            print('WARNING: sigma-clipped calculation of mean & error detected different array lengths for input & mask!')
        for i in range(0,len(input_array)):
            if mask[i]==False:
                clipped_array.append(input_array[i])
    if mask==[]:
        mask=outlierMask(input_array,sigma_clipping=sigma_clipping)
    for i in range(0,len(mask)):
        if mask[i]==False:
            clipped_array.append(input_array[i])
    mean=np.mean(clipped_array)
    error=std(clipped_array)/np.sqrt(len(clipped_array))
    return mean,error


def ratioCheck(first_array,second_array,digit,errcode=False,identifier='',max_messages=10):
    '''
    Function to check if elements of first_array match elements of second_array if rounded to
    the digit given in digit.

    Parameters
    ----------
    first_array : int/float-array
        Array to compare with second_array.
    second_array : int/float-array
        Array to compare with first_array.
    digit : int
        The digit to which the input arrays have to be rounded
    errcode : boolean, optional
        Whether to return an error code (True: mismatches found, False: no mismatches)
    identifier : string, optional
        Identifier to show in debug messages when starting the method. The default is ''.
    max_messages : int, optional
        Maximum of error messages showing individual mismatches.
        After that, mismatches will only be counted. The default is 10.

    Returns
    -------
    (error code) : boolean, optional
        If specified in errcode, this is the corresponding output.
        

    '''
    print('--Ratio check started ('+str(identifier)+')')
    e=0
    for i in range(0,len(first_array)):
        if e<max_messages:
            if round(first_array[i],digit)!=round(second_array[i],digit) and e<10:
                e=e+1
                print('----ERROR in line ', i)
                print('----',round(first_array[i],digit), ' vs. ', round(second_array[i],digit))
        if e==max_messages:
            e=e+1
            print('----POSSIBLY MORE ERRORS')   
        if e>=max_messages:
            e=e+1
        if i==len(first_array)-1 and (round(first_array[i],digit)==round(second_array[i],digit) or round(first_array[i],digit)!=round(second_array[i],digit)):
            print('----Ratio check completed, total errors: ',e)
    print()
    if e!=0 and errcode==True: return True
    elif e==0 and errcode==True: return False
    
    
def allanVariance(ratio_array):
    '''
    Function to calculate the Allan variance of an input ratio array.
    Formula taken from https://en.wikipedia.org/wiki/Allan_variance#Definitions .
    (Last inspect: 04.02.2022, 15:09)

    Parameters
    ----------
    ratio_array : int/float-array
        The array to calculate the Allan variance of.

    Returns
    -------
    output : float-array
        Array of allan variances if summed up to the corresponding element of ratio_array. 

    '''
    linsum=0
    means=[]
    for i in range(0,len(ratio_array)):
       linsum=linsum+ratio_array[i] 
       means.append(linsum/(i+1))
    output=[]
    quadsum=0
    for N in range(1,len(means)):
        quadsum=quadsum+(means[N]-means[N-1])**2
        output.append(quadsum/(2*N))
    return output


def evaluateMeasurement(path,file,output='none',plots=[],tables=[],sigma_clipping=2,wash_mode=False, blank_correction=[], flip_negatives=False, max_eic_steps=20, ratio_check=[],marker_size=5,debug=False):
    '''
    Function mainly used to calculate Sr-isotope ratios from raw isotope data,
    correct them und return them as output. Can also plot various data if
    specified accordingly.

    Parameters
    ----------
    path : string
        Path of the used file's directory. Relative path from Python's home
        directory also works (typically user directory).
    file : string
        Name of the file in the directory specified in path. Separated from
        path to allow for easier plot naming.
    output : string, optional
        String-literal that determines which data are returned by the method.
        Possible options are:
            'arrays_standard' : Returns fully-corrected (84/86,87/86) / EIC-
                corrected (88/86) or LIC-corrected (if wash mode enabled)
                arrays of the whole measurement as well as the individual
                ratios given by the Neptune over the whole measurement.
            'means_errors_standard' : Returns tuples (two-element arrays) of
                the sigma-clipped mean and error for own-calculated ratios
                (84/86,87/86: fully corrected, 88/86: EIC-only) as well as the
                corresponding tuples given by the Neptune.
            'none' : Doesn't return anything useful. ( ...abusive behavior :( )
        The default is 'none'.
    plots : string or string-array, optional
        Which plots to show. Possible options are:
            'sr84_sr86'
            'sr87_sr86'
            'sr88_sr86'
            'sr84_sr86_diff'
            'sr87_sr86_diff'
            'sr88_sr86_diff'
            'sr88_sr86_diff_kr83'
            'sr88_sr86_diff_rb85'
            'sr87_sr86_diff_rb85'
            'kr83'
            'sr84'
            'rb85'
            'sr86'
            'sr87'
            'sr88'
            'sr84_sr86_allan'
            'sr87_sr86_allan'
            'sr88_sr86_allan'
        The default is [].
    tables: string or string-array, optional
        Specifies which additional tables are saved in (separate) text files.
        The identifiers are the same as for measurement_plots. Is given over to
        evaluateSequence(). If the input is [], none are saved, if it is
        'plots', the measurement_plots input is taken in. The default is
        [].
    sigma_clipping : integer, optional
        Determines after how many multiples of sigma distance from the mean
        outliers are discarded. The default is 2.
    wash_mode : boolean, optional
        If set to True, certain operations which cannot be performed on
        negative numbers (e.g. EIC, fractionation correction) are not
        calculated nor plotted. Also the output changes! The default is False.
    blank_correction : float-tuple-array, optional
        Array of float-tuples. Each entry represents the mean blank signal that
        is going to be subtracted from the corresponding isotope and its error.
        Has to be of equal length as the number of isotopes read out / number
        of cups. In this version of the code, this corresponds to the following
        entries: 82Kr, 83Kr, 84Sr, 85Rb, 86Sr, 87Sr, 88Sr
        If set to [], no blank correction is performed. The default is [].
    flip_negatives : boolean, optional
        If set to True, negative numbers in the isotope raw data are made
        positive (absolute value operation). The default is False.
    max_eic_steps : integer, optional
        The maximum number of steps that has to be carried out when performing
        iterative element interference correction (EIC), and given that the
        correction was not already stopped before because the results did no
        longer change. When reaching this limit, the current correction will be
        taken as final one.
        The default is 20.
    ratio_check: array of three int-values
        Whether to check if the isotope ratios of 84/86, 87/86 and 88/86
        match the Neptune's values up to a specific decimal. The decimal is
        given by the corresponding entry in the array (i.e. first entry for
        84/86 and so on).
        If wash mode is set to True, ratios calculated using the linear isotope
        correction data are compared. If set to False, for 84/86 and 87/86
        ratios calculated using EIC & mass fractionation correction and for
        88/86 using EIC only are compared with the Neptune's output.
        The default is [].
    marker_size : float, optional
        Sets the size of the markers in the plots specified in plots. Is given
        over to matplotlib.pyplot. The default is 5.
    debug : boolean, optional
        Whether to show additional debug information in the console. The
        default is False.

    Returns
    -------
    float-array or float-tuple
        Sigma-clipped mean and error or whole measurement series of 84/86 data,
        corrected with own routine (usually EIC + fractionation correction, LIC 
        for wash mode).
    float-array or float-tuple
        Sigma-clipped mean and error or whole measurement series of 87/86 data,
        corrected with own routine (usually EIC + fractionation correction, LIC 
        for wash mode).
    float-array or float-tuple
        Sigma-clipped mean and error or whole measurement series of 88/86 data,
        corrected with own routine (usually EIC, LIC for wash mode).
    float-array or float-tuple
        Sigma-clipped mean and error or whole measurement series of 84/86 data,
        as given out by the Neptune.
    float-array or float-tuple
        Sigma-clipped mean and error or whole measurement series of 87/86 data,
        as given out by the Neptune.
    float-array or float-tuple
        Sigma-clipped mean and error or whole measurement series of 88/86 data,
        as given out by the Neptune.
    array of 7 float-typles, optional ('means_errors_isotopemeans' setting)
        Mean and its error for each (uncorrected) isotope signal

    '''
    print('==================================================')
    print('Refining data for file: '+path+'/'+file)
    print()
    #Read-in of data
    error=True
    prefixlines=1
    suffixlines=0
    while error==True:    
        try:
            daten = np.loadtxt(open(path+'/'+file).readlines()[:-suffixlines],dtype='str',skiprows=prefixlines)
            if daten[len(daten)-1][0] == '***':
                raise ValueError('Found possible configuration for prefix/suffix lines, but last line is still a suffix (i.e. beginning with ***)')
            print('--Prefix/Suffix lines ',prefixlines,' / ',suffixlines)
            error=False
        except:
            if prefixlines==0:
                prefixlines=suffixlines+1
                suffixlines=0
            else:
                prefixlines=prefixlines-1
                suffixlines=suffixlines+1
    cycle=[]
    time=[]
    kr82=[]
    kr83=[]
    rb85=[]
    sr84=[]
    sr86=[]
    sr87=[]
    sr88=[]
    sr84_sr86_nep=[]
    sr87_sr86_nep=[]
    sr88_sr86_nep=[]
    for i,line in enumerate(daten):
        cycle.append(line[0])
        time.append(line[1])
        kr82.append(line[2])
        kr83.append(line[3])
        sr84.append(line[4])
        rb85.append(line[5])
        sr86.append(line[6])
        sr87.append(line[7])
        sr88.append(line[8])
        sr84_sr86_nep.append(line[9])
        sr87_sr86_nep.append(line[10])
        sr88_sr86_nep.append(line[11])
        
    #Float conversion and outlier mask/check
    kr82=np.array(kr82).astype(float)
    kr83=np.array(kr83).astype(float)
    sr84=np.array(sr84).astype(float)
    rb85=np.array(rb85).astype(float)
    sr86=np.array(sr86).astype(float)
    sr87=np.array(sr87).astype(float)
    sr88=np.array(sr88).astype(float)   
    errcode,sr84_sr86_nep=outlierCheck(sr84_sr86_nep,return_converted=True,debug=debug)
    errcode,sr87_sr86_nep=outlierCheck(sr87_sr86_nep,return_converted=True,debug=debug)
    errcode,sr88_sr86_nep=outlierCheck(sr88_sr86_nep,return_converted=True,debug=debug)  
    print()
    
    #Flip negatives (if on):
    if flip_negatives==True:
        kr82=np.absolute(kr82)
        kr83=np.absolute(kr83)
        sr84=np.absolute(sr84)
        rb85=np.absolute(rb85)
        sr86=np.absolute(sr86)
        sr87=np.absolute(sr87)
        sr88=np.absolute(sr88)
    
    #Raw data:
    sr84_sr86_raw=[]
    sr87_sr86_raw=[]
    sr88_sr86_raw=[]
    for i in range(0,len(sr87)):
        sr84_sr86_raw.append(sr84[i]/sr86[i])
        sr87_sr86_raw.append(sr87[i]/sr86[i])
        sr88_sr86_raw.append(sr88[i]/sr86[i])
        
    #Calculate mean signals
    if 'isotopemeans' in output:
        isotopemeans=[
            [*clippedMeanError(kr82,sigma_clipping=sigma_clipping)],
            [*clippedMeanError(kr83,sigma_clipping=sigma_clipping)],
            [*clippedMeanError(sr84,sigma_clipping=sigma_clipping)],
            [*clippedMeanError(rb85,sigma_clipping=sigma_clipping)],
            [*clippedMeanError(sr86,sigma_clipping=sigma_clipping)],
            [*clippedMeanError(sr87,sigma_clipping=sigma_clipping)],
            [*clippedMeanError(sr88,sigma_clipping=sigma_clipping)],
        ]
        
    #Blank correction
    if blank_correction!=[]:
        print('--Performing blank correction, Sr88: ',blank_correction[6])
        for i in range(0,len(sr88)):
            kr82[i]=kr82[i]-blank_correction[0][0]
            kr83[i]=kr83[i]-blank_correction[1][0]
            sr84[i]=sr84[i]-blank_correction[2][0]
            rb85[i]=rb85[i]-blank_correction[3][0]
            sr86[i]=sr86[i]-blank_correction[4][0]
            sr87[i]=sr87[i]-blank_correction[5][0]
            sr88[i]=sr88[i]-blank_correction[6][0]
        print()
    
    #Linear isotope correction (LIC / IEC)
    sr84_lic=linearIsotopeCorrection(sr84,kr83,kr83_kr84_norm)
    sr86_lic=linearIsotopeCorrection(sr86,kr83,kr83_kr86_norm)
    sr87_lic=linearIsotopeCorrection(sr87,rb85,rb85_rb87_norm)
    if debug==True:
        print()
        sr84_lic=linearIsotopeCorrection(sr84,kr83,kr83_kr84_norm,debug=True)
        sr86_lic=linearIsotopeCorrection(sr86,kr83,kr83_kr86_norm,debug=True)
        sr87_lic=linearIsotopeCorrection(sr87,rb85,rb85_rb87_norm,debug=True)
        print()
        print('Signals without corrections (84,86,87):')
        print(sr84[0],sr86[0],sr87[0])
        print('Signals with linear isotope correction (84,86,87):')
        print(sr84_lic[0],sr86_lic[0],sr87_lic[0])
        print()
        print('Ratios without corrections:')
        print(sr84_sr86_raw[0],sr87_sr86_raw[0],sr88_sr86_raw[0])
    
    #Calculation of isotope ratios    
    sr84_sr86_lic=[]
    sr87_sr86_lic=[]
    sr88_sr86_lic=[]
    for i in range(0,len(sr87)):
        sr84_sr86_lic.append(sr84_lic[i]/sr86_lic[i])
        sr87_sr86_lic.append(sr87_lic[i]/sr86_lic[i])
        sr88_sr86_lic.append(sr88[i]/sr86_lic[i])
    if debug==True:
        print('Linear isotope correction only:')
        print(sr84_sr86_lic[0],sr87_sr86_lic[0],sr88_sr86_lic[0])
    
    #EIC using mass bias
    if wash_mode==False:
        steps=max_eic_steps
        j=0
        a=np.log(m83_kr/m86_kr)/np.log(m88_sr/m86_sr)
        sr88_sr86_eic=np.array(sr88_sr86_lic)
        sr86_eic=np.zeros(len(sr88_sr86_eic))
        previous_ratio=np.zeros(len(sr88_sr86_eic))
        while j<steps:
            print('--Perfoming 88/86 IEC, step: ',j+1)
            for i in range(0,len(sr88_sr86_eic)):
                previous_ratio[i]=sr88_sr86_eic[i]
                kr83_kr86_correction=kr83_kr86_norm*np.power(sr88_sr86_eic[i]/sr88_sr86_norm,a)
                sr86_eic[i]=sr86[i]-kr83[i]/kr83_kr86_correction
                sr88_sr86_eic[i]=sr88[i]/sr86_eic[i]
            if np.array_equal(previous_ratio, sr88_sr86_eic)==True:
                j=steps
            j=j+1
        b=np.log(m83_kr/m84_kr)/np.log(m88_sr/m86_sr)
        c=np.log(m85_rb/m87_rb)/np.log(m88_sr/m86_sr)
        sr87_eic=[]
        sr84_eic=[]
        sr87_sr86_eic=[]
        sr84_sr86_eic=[]
        for i in range(0,len(sr88_sr86_eic)):
            kr83_kr84_correction=kr83_kr84_norm*np.power(sr88_sr86_eic[i]/sr88_sr86_norm,b)
            rb85_rb87_correction=rb85_rb87_norm*np.power(sr88_sr86_eic[i]/sr88_sr86_norm,c)
            current_sr87=sr87[i]-rb85[i]/rb85_rb87_correction
            current_sr84=sr84[i]-kr83[i]/kr83_kr84_correction
            sr87_eic.append(current_sr87)
            sr84_eic.append(current_sr84)
            sr87_sr86_eic.append(current_sr87/sr86_eic[i])
            sr84_sr86_eic.append(current_sr84/sr86_eic[i])
        print()
        
    #Internal mass fractionation correction
    if wash_mode==False:
        sr87_sr86_ter=thermoMassFractionation(sr87_sr86_lic, sr88_sr86_lic, sr88_sr86_norm, [m87_sr,m86_sr], [m88_sr,m86_sr])
        sr84_sr86_ter=thermoMassFractionation(sr84_sr86_lic, sr88_sr86_lic, sr88_sr86_norm, [m84_sr,m86_sr], [m88_sr,m86_sr])
        if debug==True:    
            print('Linear Correction -> Thermo mass fractionation')
            print(sr84_sr86_lic[0],sr87_sr86_ter[0],sr88_sr86_lic[0])
            print()
        sr87_sr86_ter2=thermoMassFractionation(sr87_sr86_eic, sr88_sr86_eic, sr88_sr86_norm, [m87_sr,m86_sr], [m88_sr,m86_sr])
        sr84_sr86_ter2=thermoMassFractionation(sr84_sr86_eic, sr88_sr86_eic, sr88_sr86_norm, [m84_sr,m86_sr], [m88_sr,m86_sr])
        if debug==True:    
            print('Iterative EIC -> Thermo mass fractionation')
            print(sr84_sr86_lic[0],sr87_sr86_ter[0],sr88_sr86_lic[0])
            print()
            
    #Calculate sigma-clapped means and errors:
    sr84_sr86_raw_mean,sr84_sr86_raw_error=clippedMeanError(sr84_sr86_raw,sigma_clipping=sigma_clipping)
    sr87_sr86_raw_mean,sr87_sr86_raw_error=clippedMeanError(sr87_sr86_raw,sigma_clipping=sigma_clipping)
    sr88_sr86_raw_mean,sr88_sr86_raw_error=clippedMeanError(sr88_sr86_raw,sigma_clipping=sigma_clipping)
    sr84_sr86_lic_mean,sr84_sr86_lic_error=clippedMeanError(sr84_sr86_lic,sigma_clipping=sigma_clipping)
    sr87_sr86_lic_mean,sr87_sr86_lic_error=clippedMeanError(sr87_sr86_lic,sigma_clipping=sigma_clipping)
    sr88_sr86_lic_mean,sr88_sr86_lic_error=clippedMeanError(sr88_sr86_lic,sigma_clipping=sigma_clipping)
    sr84_sr86_nep_mean,sr84_sr86_nep_error=clippedMeanError(sr84_sr86_nep,sigma_clipping=sigma_clipping)
    sr87_sr86_nep_mean,sr87_sr86_nep_error=clippedMeanError(sr87_sr86_nep,sigma_clipping=sigma_clipping)
    sr88_sr86_nep_mean,sr88_sr86_nep_error=clippedMeanError(sr88_sr86_nep,sigma_clipping=sigma_clipping)
    if wash_mode==False:
        sr84_sr86_ter_mean,sr84_sr86_ter_error=clippedMeanError(sr84_sr86_ter,sigma_clipping=sigma_clipping)
        sr87_sr86_ter_mean,sr87_sr86_ter_error=clippedMeanError(sr87_sr86_ter,sigma_clipping=sigma_clipping)
        sr84_sr86_eic_mean,sr84_sr86_eic_error=clippedMeanError(sr84_sr86_eic,sigma_clipping=sigma_clipping)
        sr87_sr86_eic_mean,sr87_sr86_eic_error=clippedMeanError(sr87_sr86_eic,sigma_clipping=sigma_clipping)
        sr88_sr86_eic_mean,sr88_sr86_eic_error=clippedMeanError(sr88_sr86_eic,sigma_clipping=sigma_clipping)
        sr84_sr86_ter2_mean,sr84_sr86_ter2_error=clippedMeanError(sr84_sr86_ter2,sigma_clipping=sigma_clipping)
        sr87_sr86_ter2_mean,sr87_sr86_ter2_error=clippedMeanError(sr87_sr86_ter2,sigma_clipping=sigma_clipping)
         
    #Check with machine, round to specified decimal:
    if ratio_check!=[] and wash_mode==False:
        ratioCheck(sr84_sr86_ter2, sr84_sr86_nep, digit=ratio_check[0], errcode= False, identifier='84/86', max_messages=10)
        ratioCheck(sr87_sr86_ter2, sr87_sr86_nep, digit=ratio_check[1], errcode= False, identifier='87/86', max_messages=10)
        ratioCheck(sr88_sr86_eic, sr88_sr86_nep, digit=ratio_check[2], errcode= False, identifier='88/86', max_messages=10)
    if ratio_check!=[] and wash_mode==True:
        ratioCheck(sr84_sr86_lic, sr84_sr86_nep, digit=ratio_check[0], errcode= False, identifier='84/86', max_messages=10)
        ratioCheck(sr87_sr86_lic, sr87_sr86_nep, digit=ratio_check[1], errcode= False, identifier='87/86', max_messages=10)
        ratioCheck(sr88_sr86_lic, sr88_sr86_nep, digit=ratio_check[2], errcode= False, identifier='88/86', max_messages=10)
        
    #Calculate Allan variance, if needed
    if (('sr84_sr86_allan' in plots) or ('sr87_sr86_allan' in plots) or ('sr88_sr86_allan' in plots)) and wash_mode==False:
        sr84_sr86_nep_allan=allanVariance(sr84_sr86_nep)
        sr87_sr86_nep_allan=allanVariance(sr87_sr86_nep)
        sr88_sr86_nep_allan=allanVariance(sr88_sr86_nep)
        sr84_sr86_ter2_allan=allanVariance(sr84_sr86_ter2)
        sr87_sr86_ter2_allan=allanVariance(sr87_sr86_ter2)
        sr88_sr86_eic_allan=allanVariance(sr88_sr86_eic)
        sr84_sr86_raw_allan=allanVariance(sr84_sr86_raw)
        sr87_sr86_raw_allan=allanVariance(sr87_sr86_raw)
        sr88_sr86_raw_allan=allanVariance(sr88_sr86_raw)
    
    #Fix plots and tables input
    if tables=='plots':
        tables=plots
    
    #Kr82 plot/table
    if 'kr82' in plots:
        plt.figure()
        plt.plot(range(0,len(kr82)),kr82,'x',markersize=marker_size)
        plt.title('Kr82 signal in '+ripPath(path)+'/'+file[:len(file)-4])
        plt.xlabel('Cycle')
        plt.ylabel('Kr82 [V]')
        plt.tight_layout()
        plt.savefig(output_location+'/'+ripPath(path)+'_'+file[:len(file)-4]+'_kr82.jpg',dpi=200,bbox_inches='tight')        
    if 'kr82' in tables:
        table=np.matrix.transpose(np.array([cycle, kr82]))
        headerlist=['Cycle','Kr82']
        headerstring=''
        for h in headerlist:
            headerstring=headerstring+h+'\u0009'
        np.savetxt(output_location+'/neptune_auswertung_'+file[:len(file)-4]+'_kr82.csv',table,delimiter='\u0009',fmt="%s",header=headerstring)

    #Kr83 plot/table
    if 'kr83' in plots:
        plt.figure()
        plt.plot(range(0,len(kr83)),kr83,'x',markersize=marker_size)
        plt.title('Kr83 signal in '+ripPath(path)+'/'+file[:len(file)-4])
        plt.xlabel('Cycle')
        plt.ylabel('Kr83 [V]')
        plt.tight_layout()
        plt.savefig(output_location+'/'+ripPath(path)+'_'+file[:len(file)-4]+'_kr83.jpg',dpi=200,bbox_inches='tight')        
    if 'kr83' in tables:
        table=np.matrix.transpose(np.array([cycle, kr83]))
        headerlist=['Cycle','Kr83']
        headerstring=''
        for h in headerlist:
            headerstring=headerstring+h+'\u0009'
        np.savetxt(output_location+'/neptune_auswertung_'+file[:len(file)-4]+'_kr83.csv',table,delimiter='\u0009',fmt="%s",header=headerstring)

    #Sr84 plot/table
    if 'sr84' in plots:
        plt.figure()
        plt.plot(range(0,len(sr84)),sr84,'x',markersize=marker_size)
        plt.title('Sr84 signal in '+ripPath(path)+'/'+file[:len(file)-4])
        plt.xlabel('Cycle')
        plt.ylabel('Sr84 [V]')
        plt.tight_layout()
        plt.savefig(output_location+'/'+ripPath(path)+'_'+file[:len(file)-4]+'_sr84.jpg',dpi=300,bbox_inches='tight')        
    if 'sr84' in tables:
        table=np.matrix.transpose(np.array([cycle, sr84]))
        headerlist=['Cycle','Sr84 (uncorrected)']
        headerstring=''
        for h in headerlist:
            headerstring=headerstring+h+'\u0009'
        np.savetxt(output_location+'/neptune_auswertung_'+file[:len(file)-4]+'_sr84.csv',table,delimiter='\u0009',fmt="%s",header=headerstring)

    #Rb85 plot/table
    if 'rb85' in plots:
        plt.figure()
        plt.plot(range(0,len(rb85)),rb85,'x',markersize=marker_size)
        plt.title('Rb85 signal in '+ripPath(path)+'/'+file[:len(file)-4])
        plt.xlabel('Cycle')
        plt.ylabel('Rb85 [V]')
        plt.tight_layout()
        plt.savefig(output_location+'/'+ripPath(path)+'_'+file[:len(file)-4]+'_rb85.jpg',dpi=200,bbox_inches='tight')        
    if 'rb85' in tables:
        table=np.matrix.transpose(np.array([cycle, rb85]))
        headerlist=['Cycle','Rb85']
        headerstring=''
        for h in headerlist:
            headerstring=headerstring+h+'\u0009'
        np.savetxt(output_location+'/neptune_auswertung_'+file[:len(file)-4]+'_rb85.csv',table,delimiter='\u0009',fmt="%s",header=headerstring)

    #Sr86 plot/table
    if 'sr86' in plots:
        plt.figure()
        plt.plot(range(0,len(sr86)),sr86,'x',markersize=marker_size)
        plt.title('Sr86 signal in '+ripPath(path)+'/'+file[:len(file)-4])
        plt.xlabel('Cycle')
        plt.ylabel('Sr86 [V]')
        plt.tight_layout()
        plt.savefig(output_location+'/'+ripPath(path)+'_'+file[:len(file)-4]+'_sr86.jpg',dpi=200,bbox_inches='tight')        
    if 'sr86' in tables:
        table=np.matrix.transpose(np.array([cycle, sr86]))
        headerlist=['Cycle','Sr86 (uncorrected)']
        headerstring=''
        for h in headerlist:
            headerstring=headerstring+h+'\u0009'
        np.savetxt(output_location+'/neptune_auswertung_'+file[:len(file)-4]+'_sr86.csv',table,delimiter='\u0009',fmt="%s",header=headerstring)

    #Sr87 plot/table
    if 'sr87' in plots:
        plt.figure()
        plt.plot(range(0,len(sr87)),sr87,'x',markersize=marker_size)
        plt.title('Sr87 signal in '+ripPath(path)+'/'+file[:len(file)-4])
        plt.xlabel('Cycle')
        plt.ylabel('Sr87 [V]')
        plt.tight_layout()
        plt.savefig(output_location+'/'+ripPath(path)+'_'+file[:len(file)-4]+'_sr87.jpg',dpi=200,bbox_inches='tight')        
    if 'sr87' in tables:
        table=np.matrix.transpose(np.array([cycle, sr87]))
        headerlist=['Cycle','Sr87 (uncorrected)']
        headerstring=''
        for h in headerlist:
            headerstring=headerstring+h+'\u0009'
        np.savetxt(output_location+'/neptune_auswertung_'+file[:len(file)-4]+'_sr87.csv',table,delimiter='\u0009',fmt="%s",header=headerstring)

    #Sr88 plot/table
    if 'sr88' in plots:
        plt.figure()
        plt.plot(range(0,len(sr88)),sr88,'x',markersize=marker_size)
        plt.title('Sr88 signal in '+ripPath(path)+'/'+file[:len(file)-4])
        plt.xlabel('Cycle')
        plt.ylabel('Sr88 [V]')
        plt.tight_layout()
        plt.savefig(output_location+'/'+ripPath(path)+'_'+file[:len(file)-4]+'_sr88.jpg',dpi=200,bbox_inches='tight')        
    if 'sr88' in tables:
        table=np.matrix.transpose(np.array([cycle, sr88]))
        headerlist=['Cycle','Sr88']
        headerstring=''
        for h in headerlist:
            headerstring=headerstring+h+'\u0009'
        np.savetxt(output_location+'/neptune_auswertung_'+file[:len(file)-4]+'_sr88.csv',table,delimiter='\u0009',fmt="%s",header=headerstring)

    #Sr84/Sr86 plot/table
    if 'sr84_sr86' in plots:
        plt.figure()
        plt.plot(range(0,len(sr84_sr86_raw)),sr84_sr86_raw,'x',markersize=marker_size,label='Raw')
        plt.plot(range(0,len(sr84_sr86_nep)),sr84_sr86_nep,'x',markersize=marker_size,label='Neptune')
        if wash_mode==True:
            plt.plot(range(0,len(sr84_sr86_lic)),sr84_sr86_lic,'x',markersize=marker_size,label='LIC only') 
        if wash_mode==False:
            #plt.plot(range(0,len(sr84_sr86_eic)),sr84_sr86_eic,'x',markersize=marker_size,label='EIC (iteratively)',color='black')
            plt.plot(range(0,len(sr84_sr86_ter2)),sr84_sr86_ter2,'x',markersize=marker_size,label='EIC (iteratively) + frac',color='purple')        
            plt.plot(range(0,len(sr84_sr86_ter)),sr84_sr86_ter,'x',markersize=marker_size,label='LIC + frac') 
        plt.title('Sr84/Sr86 in '+ripPath(path)+'/'+file[:len(file)-4])
        plt.xlabel('Cycle')
        plt.ylabel('Sr84/Sr86')
        plt.tight_layout()
        plt.legend()
        plt.savefig(output_location+'/'+ripPath(path)+'_'+file[:len(file)-4]+'_sr84_sr86.jpg',dpi=200,bbox_inches='tight')
    if 'sr84_sr86' in tables:
        if wash_mode==True:
            table=np.matrix.transpose(np.array([cycle, sr84_sr86_lic, sr84_sr86_nep, sr84_sr86_raw]))
            headerlist=['Cycle','Sr84/Sr86 (LIC)','Sr84/Sr86 (Neptune)','Sr84/Sr86 (raw)']
        if wash_mode==False:
            table=np.matrix.transpose(np.array([cycle, sr84_sr86_ter2, sr84_sr86_nep, sr84_sr86_raw]))
            headerlist=['Cycle','Sr84/Sr86 (EIC&TER)','Sr84/Sr86 (Neptune)','Sr84/Sr86 (raw)']
        headerstring=''
        for h in headerlist:
            headerstring=headerstring+h+'\u0009'
        np.savetxt(output_location+'/neptune_auswertung_'+file[:len(file)-4]+'_sr84_sr86.csv',table,delimiter='\u0009',fmt="%s",header=headerstring)

    #Sr87/Sr86 plot/table     
    if 'sr87_sr86' in plots:
        plt.figure()
        #plt.plot(range(0,len(sr87_sr86_raw)),sr87_sr86_raw,'x',markersize=marker_size,label='Raw')
        plt.plot(range(0,len(sr87_sr86_nep)),sr87_sr86_nep,'x',markersize=marker_size,label='Neptune')
        if wash_mode==True:
            plt.plot(range(0,len(sr87_sr86_lic)),sr87_sr86_lic,'x',markersize=marker_size,label='LIC only') 
        if wash_mode==False:
            #plt.plot(range(0,len(sr87_sr86_eic)),sr87_sr86_eic,'x',markersize=marker_size,label='EIC (iteratively)',color='black')
            plt.plot(range(0,len(sr87_sr86_ter2)),sr87_sr86_ter2,'x',markersize=marker_size,label='EIC (iteratively) + frac',color='purple')        
            plt.plot(range(0,len(sr87_sr86_ter)),sr87_sr86_ter,'x',markersize=marker_size,label='LIC + frac')  
        if sr87_sr86_nep[0]>=0.65:
            plt.axhline(y=0.71034, color='grey', linestyle=':',label='SRM987 (NIST)')
            #plt.axhline(y=0.710224, color='r', linestyle='-',label='Aridus result')
            #plt.axhline(y=0.710237, color='b', linestyle='-',label='HC result')
        plt.title('Sr87/Sr86 in '+ripPath(path)+'/'+file[:len(file)-4])
        plt.xlabel('Cycle')
        plt.ylabel('Sr87/Sr86')
        plt.tight_layout()
        plt.legend()
        plt.savefig(output_location+'/'+ripPath(path)+'_'+file[:len(file)-4]+'_sr87_sr86.jpg',dpi=200,bbox_inches='tight')
    if 'sr87_sr86' in tables:
        if wash_mode==True:
            table=np.matrix.transpose(np.array([cycle, sr87_sr86_lic, sr87_sr86_nep, sr87_sr86_raw]))
            headerlist=['Cycle','Sr87/Sr86 (LIC)','Sr87/Sr86 (Neptune)','Sr87/Sr86 (raw)']
        if wash_mode==False:
            table=np.matrix.transpose(np.array([cycle, sr87_sr86_ter2, sr87_sr86_nep, sr87_sr86_raw]))
            headerlist=['Cycle','Sr87/Sr86 (EIC&TER)','Sr87/Sr86 (Neptune)','Sr87/Sr86 (raw)']
        headerstring=''
        for h in headerlist:
            headerstring=headerstring+h+'\u0009'
        np.savetxt(output_location+'/neptune_auswertung_'+file[:len(file)-4]+'_sr87_sr86.csv',table,delimiter='\u0009',fmt="%s",header=headerstring)
        
    #Sr88/Sr86 plot/table
    if 'sr88_sr86' in plots:
        plt.figure()
        plt.plot(range(0,len(sr88_sr86_nep)),sr88_sr86_nep,'+',markersize=marker_size,label='Neptune')
        plt.plot(range(0,len(sr88_sr86_raw)),sr88_sr86_raw,'x',markersize=marker_size,label='Raw')
        if wash_mode==True:
            plt.plot(range(0,len(sr88_sr86_lic)),sr88_sr86_lic,'x',markersize=marker_size,label='LIC only',color='r')
        if wash_mode==False:
            plt.plot(range(0,len(sr88_sr86_eic)),sr88_sr86_eic,'x',markersize=marker_size,label='EIC (iteratively)',color='black')
        plt.title('Sr88/Sr86 in '+ripPath(path)+'/'+file[:len(file)-4])
        plt.xlabel('Cycle')
        plt.ylabel('Sr88/Sr86')
        plt.tight_layout()
        plt.legend()
        #print(np.mean(sr88_sr86_nep),np.mean(sr88_sr86_lic))
        plt.savefig(output_location+'/'+ripPath(path)+'_'+file[:len(file)-4]+'_sr88_sr86.jpg',dpi=200,bbox_inches='tight')       
    if 'sr88_sr86' in tables:
        if wash_mode==True:
            table=np.matrix.transpose(np.array([cycle, sr88_sr86_lic, sr88_sr86_nep, sr88_sr86_raw]))
            headerlist=['Cycle','Sr88/Sr86 (LIC)','Sr88/Sr86 (Neptune)','Sr88/Sr86 (raw)']
        if wash_mode==False:
            table=np.matrix.transpose(np.array([cycle, sr88_sr86_eic, sr88_sr86_nep, sr88_sr86_raw]))
            headerlist=['Cycle','Sr88/Sr86 (EIC)','Sr88/Sr86 (Neptune)','Sr88/Sr86 (raw)']
        headerstring=''
        for h in headerlist:
            headerstring=headerstring+h+'\u0009'
        np.savetxt(output_location+'/neptune_auswertung_'+file[:len(file)-4]+'_sr88_sr86.csv',table,delimiter='\u0009',fmt="%s",header=headerstring)
        
    #Sr84/Sr86 differences plot/table
    if 'sr84_sr86_diff' in plots and wash_mode==False:
        plt.figure()
        diff=[]
        for i in range(0,len(sr84_sr86_nep)):
            diff.append(sr84_sr86_nep[i]-sr84_sr86_ter2[i])
        plt.plot(range(0,len(diff)),diff,'.',markersize=marker_size)
        plt.title('Sr84/Sr86 difference NEP-EIC&TER in '+ripPath(path)+'/'+file[:len(file)-4])
        plt.xlabel('Cycle')
        plt.ylabel('Difference NEP-EIC&TER')
        plt.tight_layout()
        #print(np.mean(sr84_sr86_nep),np.mean(sr84_sr86_ter2))
        plt.savefig(output_location+'/'+ripPath(path)+'_'+file[:len(file)-4]+'_sr84_sr86_diff.jpg',dpi=200,bbox_inches='tight')    
    if 'sr84_sr86_diff' in tables and wash_mode==False:
        diff=[]
        for i in range(0,len(sr84_sr86_nep)):
            diff.append(sr84_sr86_nep[i]-sr84_sr86_ter2[i])
        table=np.matrix.transpose(np.array([cycle, diff]))
        headerlist=['Cycle','Sr84/Sr86 difference NEP-EIC&TER']
        headerstring=''
        for h in headerlist:
            headerstring=headerstring+h+'\u0009'
        np.savetxt(output_location+'/neptune_auswertung_'+file[:len(file)-4]+'_sr84_sr86_diff.csv',table,delimiter='\u0009',fmt="%s",header=headerstring)
        
    #Sr87/Sr86 differences plot/table
    if 'sr87_sr86_diff' in plots and wash_mode==False:
        plt.figure()
        diff=[]
        for i in range(0,len(sr87_sr86_nep)):
            diff.append(sr87_sr86_nep[i]-sr87_sr86_ter2[i])
        plt.plot(range(0,len(diff)),diff,'.',markersize=marker_size)
        plt.title('Sr87/Sr86 difference NEP-EIC&TER in '+ripPath(path)+'/'+file[:len(file)-4])
        plt.xlabel('Cycle')
        plt.ylabel('Difference NEP-EIC&TER')
        plt.tight_layout()
        #print(np.mean(sr87_sr86_nep),np.mean(sr87_sr86_ter2))
        plt.savefig(output_location+'/'+ripPath(path)+'_'+file[:len(file)-4]+'_sr87_sr86_diff.jpg',dpi=200,bbox_inches='tight')    
    if 'sr87_sr86_diff' in tables and wash_mode==False:
        diff=[]
        for i in range(0,len(sr87_sr86_nep)):
            diff.append(sr87_sr86_nep[i]-sr87_sr86_ter2[i])
        table=np.matrix.transpose(np.array([cycle, diff]))
        headerlist=['Cycle','Sr87/Sr86 difference NEP-EIC&TER']
        headerstring=''
        for h in headerlist:
            headerstring=headerstring+h+'\u0009'
        np.savetxt(output_location+'/neptune_auswertung_'+file[:len(file)-4]+'_sr87_sr86_diff.csv',table,delimiter='\u0009',fmt="%s",header=headerstring)
        
    #Sr88/Sr86 differences plot/table
    if 'sr88_sr86_diff' in plots and wash_mode==False:
        plt.figure()
        diff=[]
        for i in range(0,len(sr88_sr86_nep)):
            diff.append(sr88_sr86_nep[i]-sr88_sr86_eic[i])
        plt.plot(range(0,len(diff)),diff,'.',markersize=marker_size)
        plt.title('Sr88/Sr86 difference NEP-EIC in '+ripPath(path)+'/'+file[:len(file)-4])
        plt.xlabel('Cycle')
        plt.ylabel('Difference NEP-EIC')
        plt.tight_layout()
        #print(np.mean(sr88_sr86_nep),np.mean(sr88_sr86_eic))
        plt.savefig(output_location+'/'+ripPath(path)+'_'+file[:len(file)-4]+'_sr88_sr86_diff.jpg',dpi=200,bbox_inches='tight')    
    if 'sr88_sr86_diff' in tables and wash_mode==False:
        diff=[]
        for i in range(0,len(sr88_sr86_nep)):
            diff.append(sr88_sr86_nep[i]-sr88_sr86_eic[i])
        table=np.matrix.transpose(np.array([cycle, diff]))
        headerlist=['Cycle','Sr88/Sr86 difference NEP-EIC']
        headerstring=''
        for h in headerlist:
            headerstring=headerstring+h+'\u0009'
        np.savetxt(output_location+'/neptune_auswertung_'+file[:len(file)-4]+'_sr88_sr86_diff.csv',table,delimiter='\u0009',fmt="%s",header=headerstring)
        
    #Sr88/Sr86 differences over Kr83 plot/table
    if 'sr88_sr86_diff_kr83' in plots and wash_mode==False:
        plt.figure()
        diff=[]
        for i in range(0,len(sr88_sr86_nep)):
            diff.append(sr88_sr86_nep[i]-sr88_sr86_eic[i])
        plt.plot(kr83,diff,'.',markersize=marker_size)
        plt.title('Sr88/Sr86 difference NEP-EIC over Kr83 in '+ripPath(path)+'/'+file[:len(file)-4])
        plt.xlabel('Kr83')
        plt.ylabel('Difference NEP-LIC')
        plt.tight_layout()
        plt.savefig(output_location+'/'+ripPath(path)+'_'+file[:len(file)-4]+'_sr88_sr86_diff_kr83.jpg',dpi=200,bbox_inches='tight')        
    if 'sr88_sr86_diff_kr83' in tables and wash_mode==False:
        diff=[]
        for i in range(0,len(sr88_sr86_nep)):
            diff.append(sr88_sr86_nep[i]-sr88_sr86_eic[i])
        table=np.matrix.transpose(np.array([cycle, diff, kr83]))
        headerlist=['Cycle','Sr88/Sr86 difference NEP-EIC over Kr83','Kr83']
        headerstring=''
        for h in headerlist:
            headerstring=headerstring+h+'\u0009'
        np.savetxt(output_location+'/neptune_auswertung_'+file[:len(file)-4]+'_sr88_sr86_diff_kr83.csv',table,delimiter='\u0009',fmt="%s",header=headerstring)
         
    #Sr87/Sr86 differences over Kr83 plot/table
    if 'sr87_sr86_diff_kr83' in plots and wash_mode==False:
        plt.figure()
        diff=[]
        for i in range(0,len(sr87_sr86_nep)):
            diff.append(sr87_sr86_nep[i]-sr87_sr86_ter2[i])
        plt.plot(kr83,diff,'.',markersize=marker_size)
        plt.title('Sr87/Sr86 difference NEP-EIC&TER over Kr83 in '+ripPath(path)+'/'+file[:len(file)-4])
        plt.xlabel('Kr83')
        plt.ylabel('Difference NEP-LIC')
        plt.tight_layout()
        plt.savefig(output_location+'/'+ripPath(path)+'_'+file[:len(file)-4]+'_sr87_sr86_diff_kr83.jpg',dpi=200,bbox_inches='tight')        
    if 'sr87_sr86_diff_kr83' in tables and wash_mode==False:
        diff=[]
        for i in range(0,len(sr87_sr86_nep)):
            diff.append(sr87_sr86_nep[i]-sr87_sr86_ter2[i])
        table=np.matrix.transpose(np.array([cycle, diff, kr83]))
        headerlist=['Cycle','Sr87/Sr86 difference NEP-EIC&TER over Kr83','Kr83']
        headerstring=''
        for h in headerlist:
            headerstring=headerstring+h+'\u0009'
        np.savetxt(output_location+'/neptune_auswertung_'+file[:len(file)-4]+'_sr87_sr86_diff_kr83.csv',table,delimiter='\u0009',fmt="%s",header=headerstring)
        
    #Sr88/Sr86 differences over Rb85 plot/table
    if 'sr88_sr86_diff_rb85' in plots and wash_mode==False:
        plt.figure()
        diff=[]
        for i in range(0,len(sr88_sr86_nep)):
            diff.append(sr88_sr86_nep[i]-sr88_sr86_eic[i])
        plt.plot(rb85,diff,'.',markersize=marker_size)
        plt.title('Sr88/Sr86 difference NEP-EIC over Rb85 in '+ripPath(path)+'/'+file[:len(file)-4])
        plt.xlabel('Rb85')
        plt.ylabel('Difference NEP-LIC')
        plt.tight_layout()
        plt.savefig(output_location+'/'+ripPath(path)+'_'+file[:len(file)-4]+'_sr88_sr86_diff_rb85.jpg',dpi=200,bbox_inches='tight')        
    if 'sr88_sr86_diff_rb85' in tables and wash_mode==False:
        diff=[]
        for i in range(0,len(sr88_sr86_nep)):
            diff.append(sr88_sr86_nep[i]-sr88_sr86_eic[i])
        table=np.matrix.transpose(np.array([cycle,diff, rb85]))
        headerlist=['Cycle','Sr88/Sr86 difference NEP-EIC over Rb85','Rb85']
        headerstring=''
        for h in headerlist:
            headerstring=headerstring+h+'\u0009'
        np.savetxt(output_location+'/neptune_auswertung_'+file[:len(file)-4]+'_sr88_sr86_diff_rb85.csv',table,delimiter='\u0009',fmt="%s",header=headerstring)
           
    #Sr87/Sr86 differences over Rb85 plot/table
    if 'sr87_sr86_diff_rb85' in plots and wash_mode==False:
        plt.figure()
        diff=[]
        for i in range(0,len(sr87_sr86_nep)):
            diff.append(sr87_sr86_nep[i]-sr87_sr86_ter2[i])
        plt.plot(rb85,diff,'.',markersize=marker_size)
        plt.title('Sr87/Sr86 difference NEP-EIC&TER over Rb85 in '+ripPath(path)+'/'+file[:len(file)-4])
        plt.xlabel('Rb85')
        plt.ylabel('Difference NEP-LIC')
        plt.tight_layout()
        plt.savefig(output_location+'/'+ripPath(path)+'_'+file[:len(file)-4]+'_sr87_sr86_diff_rb85.jpg',dpi=200,bbox_inches='tight')        
    if 'sr87_sr86_diff_rb85' in tables and wash_mode==False:
        diff=[]
        for i in range(0,len(sr87_sr86_nep)):
            diff.append(sr87_sr86_nep[i]-sr87_sr86_ter2[i])
        table=np.matrix.transpose(np.array([cycle, diff, rb85]))
        headerlist=['Cycle','Sr87/Sr86 difference NEP-EIC&TER over Rb85','Rb85']
        headerstring=''
        for h in headerlist:
            headerstring=headerstring+h+'\u0009'
        np.savetxt(output_location+'/neptune_auswertung_'+file[:len(file)-4]+'_sr87_sr86_diff_rb85.csv',table,delimiter='\u0009',fmt="%s",header=headerstring)
        
    #Sr84/Sr88 Allan variance plot/table
    if 'sr84_sr86_allan' in plots and wash_mode==False:
        plt.figure()
        plt.plot(range(1,len(sr84_sr86_nep_allan)+1),sr84_sr86_nep_allan,label='Neptune')
        plt.plot(range(1,len(sr84_sr86_ter2_allan)+1),sr84_sr86_ter2_allan,label='Own calculation',linestyle=':')
        plt.plot(range(1,len(sr84_sr86_raw_allan)+1),sr84_sr86_raw_allan,label='Uncorrected')
        plt.title('Allan variance of Sr84/Sr86 in '+ripPath(path)+'/'+file[:len(file)-4])
        plt.xlabel('Cycle')
        plt.xscale('log')
        plt.yscale('log')
        plt.ylabel('$(\sigma_x)^2$')
        plt.tight_layout()
        plt.legend()
        plt.savefig(output_location+'/'+ripPath(path)+'_'+file[:len(file)-4]+'_sr84_sr86_allan.jpg',dpi=200,bbox_inches='tight')   
    if 'sr84_sr86_allan' in tables and wash_mode==False:
        table=np.matrix.transpose(np.array([cycle[0:-1],sr84_sr86_ter2_allan,sr84_sr86_nep_allan,sr84_sr86_raw_allan]))
        headerlist=['Cycle','Allan variance of Sr84/Sr86 (own)','Allan variance of Sr84/Sr86 (Neptune)','Allan variance of Sr84/Sr86 (uncorrected)']
        headerstring=''
        for h in headerlist:
            headerstring=headerstring+h+'\u0009'
        np.savetxt(output_location+'/neptune_auswertung_'+file[:len(file)-4]+'_sr84_sr86_allan.csv',table,delimiter='\u0009',fmt="%s",header=headerstring)
            
    #Sr87/Sr88 Allan variance plot/table
    if 'sr87_sr86_allan' in plots and wash_mode==False:
        plt.figure()
        plt.plot(range(1,len(sr87_sr86_nep_allan)+1),sr87_sr86_nep_allan,label='Neptune')
        plt.plot(range(1,len(sr87_sr86_ter2_allan)+1),sr87_sr86_ter2_allan,label='Own calculation',linestyle=':')
        plt.plot(range(1,len(sr87_sr86_raw_allan)+1),sr87_sr86_raw_allan,label='Uncorrected')
        plt.title('Allan variance of Sr87/Sr86 in '+ripPath(path)+'/'+file[:len(file)-4])
        plt.xlabel('Cycle')
        plt.xscale('log')
        plt.yscale('log')
        plt.ylabel('$(\sigma_x)^2$')
        plt.tight_layout()
        plt.legend()
        plt.savefig(output_location+'/'+ripPath(path)+'_'+file[:len(file)-4]+'_sr87_sr86_allan.jpg',dpi=200,bbox_inches='tight')
    if 'sr87_sr86_allan' in tables and wash_mode==False:
        table=np.matrix.transpose(np.array([cycle[0:-1],sr87_sr86_ter2_allan,sr87_sr86_nep_allan,sr87_sr86_raw_allan]))
        headerlist=['Cycle','Allan variance of Sr87/Sr86 (own)','Allan variance of Sr87/Sr86 (Neptune)','Allan variance of Sr87/Sr86 (uncorrected)']
        headerstring=''
        for h in headerlist:
            headerstring=headerstring+h+'\u0009'
        np.savetxt(output_location+'/neptune_auswertung_'+file[:len(file)-4]+'_sr87_sr86_allan.csv',table,delimiter='\u0009',fmt="%s",header=headerstring)
        
    #Sr88/Sr86 Allan variance plot/table
    if 'sr88_sr86_allan' in plots and wash_mode==False:
        plt.figure()
        plt.plot(range(1,len(sr88_sr86_nep_allan)+1),sr88_sr86_nep_allan,label='Neptune')
        plt.plot(range(1,len(sr88_sr86_eic_allan)+1),sr88_sr86_eic_allan,label='Own calculation',linestyle=':')
        plt.plot(range(1,len(sr88_sr86_raw_allan)+1),sr88_sr86_raw_allan,label='Uncorrected')
        plt.title('Allan variance of Sr88/Sr86 in '+ripPath(path)+'/'+file[:len(file)-4])
        plt.xlabel('Cycle')
        plt.xscale('log')
        plt.yscale('log')
        plt.ylabel('$(\sigma_x)^2$')
        plt.tight_layout()
        plt.legend()
        plt.savefig(output_location+'/'+ripPath(path)+'_'+file[:len(file)-4]+'_sr88_sr86_allan.jpg',dpi=200,bbox_inches='tight')    
    if 'sr88_sr86_allan' in tables and wash_mode==False:
        table=np.matrix.transpose(np.array([cycle[0:-1],sr88_sr86_eic_allan,sr88_sr86_nep_allan,sr88_sr86_raw_allan]))
        headerlist=['Cycle','Allan variance of Sr88/Sr86 (own)','Allan variance of Sr88/Sr86 (Neptune)','Allan variance of Sr88/Sr86 (uncorrected)']
        headerstring=''
        for h in headerlist:
            headerstring=headerstring+h+'\u0009'
        np.savetxt(output_location+'/neptune_auswertung_'+file[:len(file)-4]+'_sr88_sr86_allan.csv',table,delimiter='\u0009',fmt="%s",header=headerstring)
        
    #Outputs
    if output=='arrays_standard':
        if wash_mode==False:
            #return sr84_sr86_ter2,sr87_sr86_ter2,sr88_sr86_eic,sr84_sr86_nep,sr87_sr86_nep,sr88_sr86_nep,sr84_sr86_nep_mask,sr87_sr86_nep_mask,sr88_sr86_nep_mask
            return sr84_sr86_ter2,sr87_sr86_ter2,sr88_sr86_eic,sr84_sr86_nep,sr87_sr86_nep,sr88_sr86_nep
        if wash_mode==True:
            #return sr84_sr86_lic,sr87_sr86_lic,sr88_sr86_lic,sr84_sr86_nep,sr87_sr86_nep,sr88_sr86_nep,sr84_sr86_nep_mask,sr87_sr86_nep_mask,sr88_sr86_nep_mask
            return sr84_sr86_lic,sr87_sr86_lic,sr88_sr86_lic,sr84_sr86_nep,sr87_sr86_nep,sr88_sr86_nep
    if output=='means_errors_standard':
        if wash_mode==False:
            print('--Outputs:')
            print('----Sr84/Sr86:',sr84_sr86_ter2_mean,'+-',sr84_sr86_ter2_error,'(Own calculation)')
            print('----Sr84/Sr86:',sr84_sr86_nep_mean,'+-',sr84_sr86_nep_error,'(Neptune)')
            print('----Sr87/Sr86:',sr87_sr86_ter2_mean,'+-',sr87_sr86_ter2_error,'(Own calculation)')
            print('----Sr87/Sr86:',sr87_sr86_nep_mean,'+-',sr87_sr86_nep_error,'(Neptune)')
            print('----Sr88/Sr86:',sr88_sr86_eic_mean,'+-',sr88_sr86_eic_error,'(Own calculation)')
            print('----Sr88/Sr86:',sr88_sr86_nep_mean,'+-',sr88_sr86_nep_error,'(Neptune)')
            print()
            return [sr84_sr86_ter2_mean,sr84_sr86_ter2_error],[sr87_sr86_ter2_mean,sr87_sr86_ter2_error],[sr88_sr86_eic_mean,sr88_sr86_eic_error],[sr84_sr86_nep_mean,sr84_sr86_nep_error],[sr87_sr86_nep_mean,sr87_sr86_nep_error],[sr88_sr86_nep_mean,sr88_sr86_nep_error]
        if wash_mode==True:
            print('--Outputs:')
            print('----Sr84/Sr86:',sr84_sr86_lic_mean,'+-',sr84_sr86_lic_error,'(Own calculation)')
            print('----Sr84/Sr86:',sr84_sr86_nep_mean,'+-',sr84_sr86_nep_error,'(Neptune)')
            print('----Sr87/Sr86:',sr87_sr86_lic_mean,'+-',sr87_sr86_lic_error,'(Own calculation)')
            print('----Sr87/Sr86:',sr87_sr86_nep_mean,'+-',sr87_sr86_nep_error,'(Neptune)')
            print('----Sr88/Sr86:',sr88_sr86_lic_mean,'+-',sr88_sr86_lic_error,'(Own calculation)')
            print('----Sr88/Sr86:',sr88_sr86_nep_mean,'+-',sr88_sr86_nep_error,'(Neptune)')
            print()
            return [sr84_sr86_lic_mean,sr84_sr86_lic_error],[sr87_sr86_lic_mean,sr87_sr86_lic_error],[sr88_sr86_lic_mean,sr88_sr86_lic_error],[sr84_sr86_nep_mean,sr84_sr86_nep_error],[sr87_sr86_nep_mean,sr87_sr86_nep_error],[sr88_sr86_nep_mean,sr88_sr86_nep_error]
    if output=='means_errors_isotopemeans':
        if wash_mode==False:
            print('--Outputs:')
            print('----Sr84/Sr86:',sr84_sr86_ter2_mean,'+-',sr84_sr86_ter2_error,'(Own calculation)')
            print('----Sr84/Sr86:',sr84_sr86_nep_mean,'+-',sr84_sr86_nep_error,'(Neptune)')
            print('----Sr87/Sr86:',sr87_sr86_ter2_mean,'+-',sr87_sr86_ter2_error,'(Own calculation)')
            print('----Sr87/Sr86:',sr87_sr86_nep_mean,'+-',sr87_sr86_nep_error,'(Neptune)')
            print('----Sr88/Sr86:',sr88_sr86_eic_mean,'+-',sr88_sr86_eic_error,'(Own calculation)')
            print('----Sr88/Sr86:',sr88_sr86_nep_mean,'+-',sr88_sr86_nep_error,'(Neptune)')
            print('----Mean Sr88: ', isotopemeans[6][0],'+-',isotopemeans[6][1])
            print('----Mean Kr83: ', isotopemeans[1][0],'+-',isotopemeans[1][1])
            print('----Mean Rb85: ', isotopemeans[3][0],'+-',isotopemeans[3][1])
            print()
            return [sr84_sr86_ter2_mean,sr84_sr86_ter2_error],[sr87_sr86_ter2_mean,sr87_sr86_ter2_error],[sr88_sr86_eic_mean,sr88_sr86_eic_error],[sr84_sr86_nep_mean,sr84_sr86_nep_error],[sr87_sr86_nep_mean,sr87_sr86_nep_error],[sr88_sr86_nep_mean,sr88_sr86_nep_error],isotopemeans
        if wash_mode==True:
            print('--Outputs:')
            print('----Sr84/Sr86:',sr84_sr86_lic_mean,'+-',sr84_sr86_lic_error,'(Own calculation)')
            print('----Sr84/Sr86:',sr84_sr86_nep_mean,'+-',sr84_sr86_nep_error,'(Neptune)')
            print('----Sr87/Sr86:',sr87_sr86_lic_mean,'+-',sr87_sr86_lic_error,'(Own calculation)')
            print('----Sr87/Sr86:',sr87_sr86_nep_mean,'+-',sr87_sr86_nep_error,'(Neptune)')
            print('----Sr88/Sr86:',sr88_sr86_lic_mean,'+-',sr88_sr86_lic_error,'(Own calculation)')
            print('----Sr88/Sr86:',sr88_sr86_nep_mean,'+-',sr88_sr86_nep_error,'(Neptune)')
            print('----Mean Sr88: ', isotopemeans[6][0],'+-',isotopemeans[6][1])
            print('----Mean Kr83: ', isotopemeans[1][0],'+-',isotopemeans[1][1])
            print('----Mean Rb85: ', isotopemeans[3][0],'+-',isotopemeans[3][1])
            print()
            return [sr84_sr86_lic_mean,sr84_sr86_lic_error],[sr87_sr86_lic_mean,sr87_sr86_lic_error],[sr88_sr86_lic_mean,sr88_sr86_lic_error],[sr84_sr86_nep_mean,sr84_sr86_nep_error],[sr87_sr86_nep_mean,sr87_sr86_nep_error],[sr88_sr86_nep_mean,sr88_sr86_nep_error],isotopemeans


def standardBracket(alldata,measurement_types,uncorrected_index,correction_index,correction_factor=sr88_sr86_norm,cycle_loop=False):
    '''
    Function to perform standard bracketing on a dataset of standard before and after sample measurements.

    Parameters
    ----------
    alldata : array of arrays of float-tuples (=two-element arrays)
        The dataset to perform standard bracketing on. Has to be structured as follows:
        -> array of measurements
            -> contains arrays of different corrected isotope ratios
                -> contain the ratio and its error as tuple/two-element array
    measurement_types : string-array
        Array of strings, determines the type of measurement at the corresponding position in data.
        Possible inputs are: 'wash', 'water', 'hfwash', 'standard', 'sample', 'drift', 'blank'
    uncorrected_index : int
        Index of the ratio within a measurement (i.e. second index of data) of the ratio to correct.
    correction_index : int
        Index of the ratio within a measurement (i.e. second index of data) of the ratio to correct with.
    correction_factor : float, optional
        The known natural ratio corresponding to the ratio of correction_index.
        The default is sr88_sr86_norm.
    cycle_loop : boolean, optional
        This is a special option in case the structure of alldata looks differently: Instead of performing
        standard bracketing on already calculated means of a measurement, one can loop over all individual
        cycles and bracket each cycle. Only usable for measurement that switch between sample and standard
        in each cycle. The default is False.

    Returns
    -------
    corrected_ratios : float-array
        Array of the corrected 'sample'-type measurements' ratios
    d_corrected_ratios : float-array
        Array of the corrected 'sample'-type measurements ratio's errors

    '''
    sample_indices=[]
    corrected_ratios=[]
    d_corrected_ratios=[]
    #Find samples
    for i in range(0,len(alldata)):
        if measurement_types[i]=='sample':
            sample_indices.append(i)
    for p in sample_indices:
        #Find standards
        standard_index_1=-1
        standard_index_2=-1
        j=p-1
        while j>=0:
            if measurement_types[j]=='standard':
                standard_index_1=j
                break
            elif j==0:
                print('WARNING: No standard before sample measurement ',p,' found!')
            j=j-1
        j=p+1
        while j<len(measurement_types):
            if measurement_types[j]=='standard':
                standard_index_2=j
                break
            elif j==len(measurement_types)-1:
                print('WARNING: No standard after sample measurement ',p,' found!')
            j=j+1
        #Calculate bracketing
        if cycle_loop==True:
            cycles=[]
            d_cycles=[]
            for i in range(0,len(alldata[p][uncorrected_index][0])):
                cycles.append(alldata[p][uncorrected_index][0][i]*correction_factor/np.sqrt(alldata[standard_index_1][correction_index]*alldata[standard_index_2][correction_index]))
                d_cycles.append(np.sqrt((correction_factor*alldata[p][uncorrected_index][1]/np.sqrt(alldata[standard_index_1][correction_index][0]*alldata[standard_index_2][correction_index][0]))**2+(correction_factor*alldata[p][uncorrected_index][0]*alldata[standard_index_1][correction_index][1]/(2*np.sqrt(alldata[standard_index_2][correction_index][0]*alldata[standard_index_1][correction_index][0]**3)))**2+(correction_factor*alldata[p][uncorrected_index][0]*alldata[standard_index_2][correction_index][1]/(2*np.sqrt(alldata[standard_index_1][correction_index][0]*alldata[standard_index_2][correction_index][0]**3)))**2))
            corrected_ratios.append(cycles)
            d_corrected_ratios.append(d_cycles)
        elif cycle_loop==False:
            corrected_ratios.append(alldata[p][uncorrected_index][0]*correction_factor/np.sqrt(alldata[standard_index_1][correction_index][0]*alldata[standard_index_2][correction_index][0]))
            error=np.sqrt((correction_factor*alldata[p][uncorrected_index][1]/np.sqrt(alldata[standard_index_1][correction_index][0]*alldata[standard_index_2][correction_index][0]))**2+(correction_factor*alldata[p][uncorrected_index][0]*alldata[standard_index_1][correction_index][1]/(2*np.sqrt(alldata[standard_index_2][correction_index][0]*alldata[standard_index_1][correction_index][0]**3)))**2+(correction_factor*alldata[p][uncorrected_index][0]*alldata[standard_index_2][correction_index][1]/(2*np.sqrt(alldata[standard_index_1][correction_index][0]*alldata[standard_index_2][correction_index][0]**3)))**2)
            d_corrected_ratios.append(error)
    return corrected_ratios,d_corrected_ratios
            

def evaluateSequence(path,measurement_types=[],blank_correction=True,sequence_plots=[],sequence_tables=[],plot_selection='standards_samples',show_neptune=False,measurement_plots=[],measurement_tables=[],sigma_clipping=2,boundaries=[23,11]):
    '''
    Function used to evaluate an entire sequence of Sr isotope data, stored
    separately in its own folder. Uses evaluateMeasurement() to gain corrected
    isotope ratios and standardBracket() to perform standard bracketing on
    them. Is able to plot specified data (and give plot specifications over to
    evaluateMeasurement()).

    Parameters
    ----------
    path : string
        The path of the folder in which the raw isotope data are stored. They
        need to be .exp-files as exported by the Neptune Evaluate software.
        It is required to not keep any other files in this folder!
    measurement_types : string-array, optional
        An array with the same length as the number of measurements in path.
        Specifies the type of the corresponding measurement (if sorted top-down
        by name). Possible elements are:
            'standard'
            'sample'
            'blank'
            'wash'
            'hfwash'
            'water'
        If not specified, the method will attempt to generate such an array
        from the suffix of the file name, separated by an _ from the rest.
        The default is [].
    blank_correction : boolean
        Whether to perform a blank correction or not. For that, the mean
        signals of the previous blank measurements are subtracted before each
        sample or standard type measurement. The error of that is assumed to
        be negligible. The default is True.
    sequence_plots : string or string-array, optional
        Specifies which data are to be plotted, over the whole sequence as
        x-axis. Depending on the objective, either a string or an array of 
        strings from the following inputs are needed:
            'sr84_sr86'
            'sr87_sr86'
            'sr87_sr86_abserror'
            'sr87_sr86_relerror'
            'sr87_sr86_difference'
            'sr88_sr86'
        The default is [].
    sequence_tables : string-array
        Specifies which additional tables are saved in (separate) text files.
        Possible identifiers are:
            'sr84_sr86'
            'sr87_sr86'
            'sr88_sr86'
        Is given over to evaluateSequence(). If the input is [], none are
        saved. The default is [].
    plot_selection : string, optional
        Which measurements to show in the plots specified in sequence_plots.
        This is especially useful if one wants to exclude washes from being
        plotted. Note that this also excludes them from being calculated.
        Possible options are:
        'all'
        'standards_samples'
        'samples'
        The default is 'standards_samples'.
    show_neptune : boolean, optional
        Whether to show the Neptune's ratios as well in the sequence plots.
        Note that that may not have any blank correction (at least none made
        in this code). The default is False.
    measurement_plots : string or string-array, optional
        Specifies which data are to be plotted in one plot for each measurement
        with the x-axis representing the individual cycles. Essentially this is
        simply an additional plot input (as in sequence_plots) that is given
        over to evaluateMeasurement(). The default is [].
    measurement_tables: string or string-array, optional
        Specifies which additional tables are saved in (separate) text files.
        - Errors are always saved as 1 sigma! -
        The identifiers are the same as for measurement_plots. Is given over to
        evaluateSequence(). If the input is [], none are saved, if it is
        'plots', the measurement_plots input is taken in. The default is
        [].
    sigma_clipping : integer, optional
        The distance from array mean in standard deviations after which a data
        point is considered an outlier and therefore not considered. Is given
        over to evaluateMeasurement() and its dependencies. The default is 2.

    Returns
    -------
    None. Writes results into neptune_auswertung.csv in output directory and
    makes specified plots and tables, hopefully :)

    '''
    #Remove zombie output
    for file in listdir(output_location):
        remove(output_location+'/'+file)
    
    #Select correct files
    files=[]
    for file in listdir(path):
        if file[-4:] == '.exp':
            files.append(file)
    
    #Correct missing measurement types
    if measurement_types==[]:
        for file in files:
            if 'wash' in file or 'Wash' in file:
                measurement_types.append('wash')
            if 'hfw' in file or 'hfwash' in file or 'HF-Wash' in file:
                measurement_types.append('hfwash')
            if 'blank' in file or 'BLK' in file:
                measurement_types.append('blank')
            if 'STD' in file or 'std' in file or 'standard' in file:
                measurement_types.append('standard')
            if 'SMP' in file or 'smp' in file or 'sample' in file:
                measurement_types.append('sample')
    
    #Get data & perform blank correction
    alldata=[]
    for f in files:
        alldata.append([])
    nextblanks=[[0,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0]]
    washmodetypes=['hfwash','wash','water','blank']
    for m in range(0,len(files)):
        if (measurement_types[m] in washmodetypes) == True:
            out=[*evaluateMeasurement(path,files[m],output='means_errors_isotopemeans',plots=measurement_plots,tables=measurement_tables,sigma_clipping=2,wash_mode=True,blank_correction=[],flip_negatives=False,ratio_check=[8,6,9],debug=False)]
            alldata[m]=out
            nextblanks=out[len(out)-1]
        if (measurement_types[m] in washmodetypes) == False and blank_correction == True:
            out=[*evaluateMeasurement(path,files[m],output='means_errors_isotopemeans',plots=measurement_plots,tables=measurement_tables,sigma_clipping=2,wash_mode=False,blank_correction=nextblanks,flip_negatives=False,ratio_check=[],debug=False)]
            alldata[m]=out
        if (measurement_types[m] in washmodetypes) == False and blank_correction == False:
            out=[*evaluateMeasurement(path,files[m],output='means_errors_isotopemeans',plots=measurement_plots,tables=measurement_tables,sigma_clipping=2,wash_mode=False,blank_correction=[],flip_negatives=False,ratio_check=[],debug=False)]
            alldata[m]=out
    
    #Standard bracketing
    true_sr84_sr86,d_true_sr84_sr86=standardBracket(alldata=alldata, measurement_types=measurement_types, uncorrected_index=0, correction_index=0, correction_factor=sr84_sr86_srm987, cycle_loop=False)
    true_sr87_sr86,d_true_sr87_sr86=standardBracket(alldata=alldata, measurement_types=measurement_types, uncorrected_index=1, correction_index=1, correction_factor=sr87_sr86_srm987, cycle_loop=False)
    true_sr88_sr86,d_true_sr88_sr86=standardBracket(alldata=alldata, measurement_types=measurement_types, uncorrected_index=2, correction_index=2, correction_factor=sr88_sr86_srm987, cycle_loop=False)      
    
    #Measurement selection for plots
    sample_indices=[]
    standards_samples_indices=[]
    for i in range(0,len(measurement_types)):
        if measurement_types[i]=='standard' or measurement_types[i]=='sample':
            standards_samples_indices.append(i)
    for i in range(0,len(measurement_types)):
        if measurement_types[i]=='sample':
            sample_indices.append(i)
    if plot_selection=='all':
        plot_indices=range(0,len(measurement_types))
    if plot_selection=='standards_samples':
        plot_indices=standards_samples_indices
    if plot_selection=='samples':
        plot_indices=sample_indices
      
    #Plot Sr84/Sr88 means:
    if ('sr84_sr86' in sequence_plots) or ('sr84_sr86' in sequence_tables):
        y1=[]
        y2=[]
        y3=[]
        d_y1=[]
        d_y2=[]
        d_y3=[]
        j=0
        for i in plot_indices:
            y1.append(alldata[i][0][0])
            y2.append(alldata[i][3][0])
            d_y1.append(alldata[i][0][1])
            d_y2.append(alldata[i][3][1])
            if i in sample_indices:
                y3.append(true_sr84_sr86[j])
                d_y3.append(d_true_sr84_sr86[j])
                j=j+1
    if 'sr84_sr86' in sequence_plots:
        plt.figure()
        plt.errorbar(plot_indices,y2,2*np.array(d_y2),fmt='.',label='Neptune')
        plt.errorbar(plot_indices,y1,2*np.array(d_y1),fmt='.',label='Own Calculation')
        if sample_indices!=[]:
            plt.errorbar(sample_indices,y3,2*np.array(d_y3),fmt='.',label='Standard bracketing')
        plt.axhline(y=sr84_sr86_srm987, color='grey', linestyle=':',label='SRM987 (NIST)')
        plt.title('Sr84/Sr86 in '+ripPath(path))
        plt.xlabel('No. of measurement')
        plt.ylabel('Sr84/Sr86')
        plt.tight_layout()
        plt.legend()
        plt.savefig(output_location+'/ALL_'+ripPath(path)+'_sr84_sr86.jpg',dpi=200)
    if 'sr84_sr86' in sequence_tables:
        table=[]
        j=0
        for i in plot_indices:
            if i in sample_indices:
                table.append([i,alldata[i][0][0],2*alldata[i][0][1],alldata[i][3][0],2*alldata[i][3][1],y3[j],2*d_y3[j]])
                j=j+1
            else:
                table.append([i,alldata[i][0][0],2*alldata[i][0][1],alldata[i][3][0],2*alldata[i][3][1],alldata[i][0][0],2*alldata[i][0][1]])
        headerlist=['Cycle','Own ratio','(2s error)','Neptune ratio','(2s error)','Bracketed own','(2s error)']
        headerstring=''
        for h in headerlist:
            headerstring=headerstring+h+'\u0009'
        np.savetxt(output_location+'/ALL_'+ripPath(path)+'_sr84_sr86.csv',table,delimiter='\u0009',fmt="%s",header=headerstring)
        
    #Plot Sr84/Sr86 errors (relative):    
    if 'sr84_sr86_relerror' in sequence_plots:
        y1=[]
        y2=[]
        y3=[]
        d_y1=[]
        d_y2=[]
        d_y3=[]
        j=0
        for i in plot_indices:
            y1.append(alldata[i][0][0])
            y2.append(alldata[i][3][0])
            d_y1.append(alldata[i][0][1])
            d_y2.append(alldata[i][3][1])
            if i in sample_indices:
                y3.append(true_sr84_sr86[j])
                d_y3.append(d_true_sr84_sr86[j])
                j=j+1
        plt.figure()
        plt.errorbar(plot_indices,2*np.array(d_y2)/np.array(y2),fmt='.',label='Neptune')
        plt.errorbar(plot_indices,2*np.array(d_y1)/np.array(y1),fmt='.',label='Own Calculation')
        if sample_indices!=[]:
            plt.errorbar(sample_indices,np.array(d_y3)/np.array(y3),fmt='.',label='Standard bracketing')
        plt.title('Sr84/Sr86 relative error ($2\sigma$) in '+ripPath(path))
        plt.xlabel('No. of measurement')
        plt.ylabel('$\Delta$Sr84/Sr86 (rel.) ($2\sigma$)')
        plt.tight_layout()
        plt.legend()
        plt.savefig(output_location+'/ALL_'+ripPath(path)+'_sr84_sr86_relerror.jpg',dpi=200)

    #Plot Sr87/Sr86 means:    
    if ('sr87_sr86' in sequence_plots) or ('sr87_sr86' in sequence_tables):
        y1=[]
        y2=[]
        y3=[]
        d_y1=[]
        d_y2=[]
        d_y3=[]
        j=0
        for i in plot_indices:
            y1.append(alldata[i][1][0])
            y2.append(alldata[i][4][0])
            d_y1.append(alldata[i][1][1])
            d_y2.append(alldata[i][4][1])
            if i in sample_indices:
                y3.append(true_sr87_sr86[j])
                d_y3.append(d_true_sr87_sr86[j])
                j=j+1
    if 'sr87_sr86' in sequence_plots:
        plt.figure()
        plt.errorbar(plot_indices,y2,2*np.array(d_y2),fmt='.',label='Neptune')
        plt.errorbar(plot_indices,y1,2*np.array(d_y1),fmt='.',label='Own Calculation')
        if sample_indices!=[]:
            plt.errorbar(sample_indices,y3,2*np.array(d_y3),fmt='.',label='Standard bracketing')
        plt.axhline(y=0.71034, color='grey', linestyle=':',label='SRM987 (NIST)')
        plt.axhline(y=0.709171, color='steelblue', linestyle=':',label='IAPSO seawater (Voigt 2015)')
        plt.title('Sr87/Sr86 in '+ripPath(path))
        plt.xlabel('No. of measurement')
        plt.ylabel('Sr87/Sr86')
        plt.tight_layout()
        plt.legend()
        plt.savefig(output_location+'/ALL_'+ripPath(path)+'_sr87_sr86.jpg',dpi=200)
    if 'sr87_sr86' in sequence_tables:
        table=[]
        j=0
        for i in plot_indices:
            if i in sample_indices:
                table.append([i,alldata[i][1][0],2*alldata[i][1][1],alldata[i][4][0],2*alldata[i][4][1],y3[j],2*d_y3[j]])
                j=j+1
            else:
                table.append([i,alldata[i][1][0],2*alldata[i][1][1],alldata[i][4][0],2*alldata[i][4][1],alldata[i][1][0],2*alldata[i][1][1]])
        headerlist=['Cycle','Own ratio','(2s error)','Neptune ratio','(2s error)','Bracketed own','(2s error)']
        headerstring=''
        for h in headerlist:
            headerstring=headerstring+h+'\u0009'
        np.savetxt(output_location+'/ALL_'+ripPath(path)+'_sr87_sr86.csv',table,delimiter='\u0009',fmt="%s",header=headerstring)
        
    #Plot Sr87/Sr86 errors (absolute):    
    if 'sr87_sr86_abserror' in sequence_plots:
        d_y1=[]
        d_y2=[]
        d_y3=[]
        j=0
        for i in plot_indices:
            d_y1.append(2*alldata[i][1][1])
            d_y2.append(2*alldata[i][4][1])
            if i in sample_indices:
                y3.append(true_sr87_sr86[j])
                d_y3.append(2*d_true_sr87_sr86[j])
                j=j+1
        plt.figure()
        plt.errorbar(plot_indices,d_y2,fmt='.',label='Neptune')
        plt.errorbar(plot_indices,d_y1,fmt='.',label='Own Calculation')
        if sample_indices!=[]:
            plt.errorbar(sample_indices,d_y3,fmt='.',label='Standard bracketing')
        plt.title('Sr87/Sr86 absolute error ($2\sigma$) in '+ripPath(path))
        plt.xlabel('No. of measurement')
        plt.ylabel('$\Delta$Sr87/Sr86 (abs.) ($2\sigma$)')
        plt.tight_layout()
        plt.legend()
        plt.savefig(output_location+'/ALL_'+ripPath(path)+'_sr87_sr86_abserror.jpg',dpi=200)
        
    #Plot Sr87/Sr86 errors (relative):    
    if 'sr87_sr86_relerror' in sequence_plots:
        y1=[]
        y2=[]
        y3=[]
        d_y1=[]
        d_y2=[]
        d_y3=[]
        j=0
        for i in plot_indices:
            y1.append(alldata[i][1][0])
            y2.append(alldata[i][4][0])
            d_y1.append(alldata[i][1][1])
            d_y2.append(alldata[i][4][1])
            if i in sample_indices:
                y3.append(true_sr87_sr86[j])
                d_y3.append(d_true_sr87_sr86[j])
                j=j+1
        plt.figure()
        plt.errorbar(plot_indices,2*np.array(d_y2)/np.array(y2),fmt='.',label='Neptune')
        plt.errorbar(plot_indices,2*np.array(d_y1)/np.array(y1),fmt='.',label='Own Calculation')
        if sample_indices!=[]:
            plt.errorbar(sample_indices,np.array(d_y3)/np.array(y3),fmt='.',label='Standard bracketing')
        plt.title('Sr87/Sr86 relative error ($2\sigma$) in '+ripPath(path))
        plt.xlabel('No. of measurement')
        plt.ylabel('$\Delta$Sr87/Sr86 (rel.) ($2\sigma$)')
        plt.tight_layout()
        plt.legend()
        plt.savefig(output_location+'/ALL_'+ripPath(path)+'_sr87_sr86_relerror.jpg',dpi=200)
        
    #Plot differences in 87/86 between own and Neptune calculation
    if 'sr87_sr86_difference' in sequence_plots:
        y1=[]
        y2=[]
        y3=[]
        j=0
        for i in plot_indices:
            y1.append(alldata[i][1][0])
            y2.append(alldata[i][4][0])
        plt.figure()
        plt.errorbar(plot_indices,np.array(y1)-np.array(y2),fmt='.',label='Own Calculation - NEP')
        plt.title('Sr87/Sr86 differences to Neptune in '+ripPath(path))
        plt.xlabel('No. of measurement')
        plt.ylabel('Sr87/Sr86')
        plt.legend()
        plt.tight_layout()
        plt.savefig(output_location+'/ALL_'+ripPath(path)+'_sr87_sr86_difference.jpg',dpi=200)
        
    #Plot Sr88/Sr86 means:
    if ('sr88_sr86' in sequence_plots) or ('sr88_sr86' in sequence_tables):
        y1=[]
        y2=[]
        y3=[]
        d_y1=[]
        d_y2=[]
        d_y3=[]
        j=0
        for i in plot_indices:
            y1.append(alldata[i][2][0])
            y2.append(alldata[i][5][0])
            d_y1.append(alldata[i][2][1])
            d_y2.append(alldata[i][5][1])
            if i in sample_indices:
                y3.append(true_sr88_sr86[j])
                d_y3.append(d_true_sr88_sr86[j])
                j=j+1
    if 'sr88_sr86' in sequence_plots:
        plt.figure()
        plt.errorbar(plot_indices,y2,2*np.array(d_y2),fmt='.',label='Neptune')
        plt.errorbar(plot_indices,y1,2*np.array(d_y1),fmt='.',label='Own Calculation')
        if sample_indices!=[]:
            plt.errorbar(sample_indices,y3,d_y3,fmt='.',label='Standard bracketing')
        plt.axhline(sr88_sr86_srm987, color='grey', linestyle=':',label='SRM987 (NIST)')
        plt.title('Sr88/Sr86 in '+ripPath(path))
        plt.xlabel('No. of measurement')
        plt.ylabel('Sr88/Sr86')
        plt.tight_layout()
        plt.legend()
        plt.savefig(output_location+'/ALL_'+ripPath(path)+'_sr88_sr86.jpg',dpi=200)
    if 'sr88_sr86' in sequence_tables:
        table=[]
        j=0
        for i in plot_indices:
            if i in sample_indices:
                table.append([i,alldata[i][2][0],2*alldata[i][2][1],alldata[i][5][0],2*alldata[i][5][1],y3[j],2*d_y3[j]])
                j=j+1
            else:
                table.append([i,alldata[i][2][0],2*alldata[i][2][1],alldata[i][5][0],2*alldata[i][5][1],alldata[i][2][0],2*alldata[i][2][1]])
        headerlist=['Cycle','Own ratio','(2s error)','Neptune ratio','(2s error)','Bracketed own','(2s error)']
        headerstring=''
        for h in headerlist:
            headerstring=headerstring+h+'\u0009'
        np.savetxt(output_location+'/ALL_'+ripPath(path)+'_sr88_sr86.csv',table,delimiter='\u0009',fmt="%s",header=headerstring)
        
    #Plot Sr88/Sr86 errors (relative):    
    if 'sr88_sr86_relerror' in sequence_plots:
        y1=[]
        y2=[]
        y3=[]
        d_y1=[]
        d_y2=[]
        d_y3=[]
        j=0
        for i in plot_indices:
            y1.append(alldata[i][2][0])
            y2.append(alldata[i][5][0])
            d_y1.append(alldata[i][2][1])
            d_y2.append(alldata[i][5][1])
            if i in sample_indices:
                y3.append(true_sr88_sr86[j])
                d_y3.append(d_true_sr88_sr86[j])
                j=j+1
        plt.figure()
        plt.errorbar(plot_indices,2*np.array(d_y2)/np.array(y2),fmt='.',label='Neptune')
        plt.errorbar(plot_indices,2*np.array(d_y1)/np.array(y1),fmt='.',label='Own Calculation')
        if sample_indices!=[]:
            plt.errorbar(sample_indices,np.array(d_y3)/np.array(y3),fmt='.',label='Standard bracketing')
        plt.title('Sr88/Sr86 relative error ($2\sigma$) in '+ripPath(path))
        plt.xlabel('No. of measurement')
        plt.ylabel('$\Delta$Sr88/Sr86 (rel.) ($2\sigma$)')
        plt.tight_layout()
        plt.legend()
        plt.savefig(output_location+'/ALL_'+ripPath(path)+'_sr88_sr86_relerror.jpg',dpi=200)
        
    #Save bracketed results
    #NOTE: Tables always use 1-sigma errors!!!
    output=[]
    j=0
    for i in range(0,len(alldata)):
        if i in sample_indices:
            output.append([i,measurement_types[i],true_sr84_sr86[j],d_true_sr84_sr86[j],true_sr87_sr86[j],d_true_sr87_sr86[j],true_sr88_sr86[j],d_true_sr88_sr86[j],alldata[i][6][6][0],alldata[i][6][6][1],alldata[i][6][3][0],alldata[i][6][3][1],alldata[i][6][1][0],alldata[i][6][1][1]])
            j=j+1
        else:
            output.append([i,measurement_types[i],alldata[i][0][0],alldata[i][0][1],alldata[i][1][0],alldata[i][1][1],alldata[i][2][0],alldata[i][2][1],alldata[i][6][6][0],alldata[i][6][6][1],alldata[i][6][3][0],alldata[i][6][3][1],alldata[i][6][1][0],alldata[i][6][1][1]])
    headerlist=['Index','Meas. Type','Sr84/86','(error)','Sr87/86','(error)','Sr88/86','(error)','Sr88','(error)','Rb85','(error)','Kr83','(error)']
    headerstring=''
    for h in headerlist:
        headerstring=headerstring+h+'\u0009'
    np.savetxt(output_location+'/neptune_auswertung_output_bracketing.csv',output,delimiter='\u0009',fmt="%s",header=headerstring)
    
    #Save unbracketed results
    #NOTE: Tables always use 1-sigma errors!!!
    output=[]
    j=0
    for i in range(0,len(alldata)):
        output.append([i,measurement_types[i],alldata[i][0][0],alldata[i][0][1],alldata[i][1][0],alldata[i][1][1],alldata[i][2][0],alldata[i][2][1],alldata[i][6][6][0],alldata[i][6][6][1],alldata[i][6][3][0],alldata[i][6][3][1],alldata[i][6][1][0],alldata[i][6][1][1]])
    headerlist=['Index','Meas. Type','Sr84/86','(error)','Sr87/86','(error)','Sr88/86','(error)','Sr88','(error)','Rb85','(error)','Kr83','(error)']
    headerstring=''
    for h in headerlist:
        headerstring=headerstring+h+'\u0009'
    np.savetxt(output_location+'/neptune_auswertung_output_nobracketing.csv',output,delimiter='\u0009',fmt="%s",header=headerstring)
    
    #Save Neptune results
    #NOTE: Tables always use 1-sigma errors!!!
    output=[]
    j=0
    for i in range(0,len(alldata)):
        output.append([i,measurement_types[i],alldata[i][3][0],alldata[i][3][1],alldata[i][4][0],alldata[i][4][1],alldata[i][5][0],alldata[i][5][1]])
    headerlist=['Index','Meas. Type','Sr84/86','(error)','Sr87/86','(error)','Sr88/86','(error)']
    headerstring=''
    for h in headerlist:
        headerstring=headerstring+h+'\u0009'
    np.savetxt(output_location+'/neptune_auswertung_output_neptune.csv',output,delimiter='\u0009',fmt="%s",header=headerstring)

    #Save current Python script
    copyfile(current_location,output_location+'/RESTREX.py')



###############################################################################


#Only change location at first, and run
#The rest should be fine for most applications :-)
#Don't forget to change the directory fields in the beginning when moving the file!

location='Desktop/RESTREX_Starterkit/Example_Data'

evaluateSequence(path=location,measurement_types=[],
                 blank_correction=True,
                 sequence_plots=['sr84_sr86','sr87_sr86','sr88_sr86'],
                 sequence_tables=[],
                 plot_selection='standards_samples',
                 show_neptune=False,
                 measurement_plots=['kr83','rb85','sr84','sr88'],
                 measurement_tables=['sr88'],
                 sigma_clipping=2)
