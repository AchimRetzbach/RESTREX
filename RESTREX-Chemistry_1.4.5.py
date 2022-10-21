# -*- coding: utf-8 -*-
"""
Created on Fri Dec 17 13:59:53 2021

@author: Achim Retzbach

@version: 1.4.5
"""


import numpy as np
import matplotlib.pyplot as plt
from os import listdir
from os import remove
from shutil import copyfile
from tabulate import tabulate
import scipy.odr as odr


plt.rcParams['axes.titley'] = 1.03
plt.rcParams['axes.formatter.limits'] = (-3,4)
plt.rcParams['figure.figsize'] = plt.rcParamsDefault['figure.figsize']
#plt.rcParams["figure.figsize"] = (12,8)


#Fields / fixed variables
current_location='Desktop/RESTREX_Starterkit/RESTREX-Chemistry_1.4.5.py'
output_location='Desktop/RESTREX_Starterkit/Output_Chemistry'


def linear(B,x):
    '''Linear function.'''
    return B[0]*x+B[1]


def proportional(B,x):
    '''Proportional function.'''
    return B[0]*x


def stringToFloat(s):
    '''
    Simple function to convert a single input string into a float number. If
    this is not possible, the rror is caught and 0 is output instead. (i.e. it
    'bulldozes' nans). Be aware that this error catching can lead to serious
    difficulties if the input is not as intended.

    Parameters
    ----------
    s : string
        String that is tried to be converted. No, not to Islam, to integer.

    Returns
    -------
    out : float
        Float number representing the string from before (if possible,
        otherwise zero)

    '''
    try:
        out=float(s)
    except:
        out=0
    return out


def sliceArray(array,slices):
    '''
    Simple function to slice an array according to a list of slice indices.

    Parameters
    ----------
    array : array of arbitrary type
        The array to slice.
    slices : array of two-element integer-arrays (i.e. tuples) / integers
        Single integers in this array will be interpreted as indices at which
        data from array are sliced out. Tuples will be interpreted as intervals
        of indices between which data from array are sliced out. Note that the
        former index within such a tuple is the first to be used while the
        latter is the first index to not be used anymore.

    Returns
    -------
    output : array of arbitrary type
        Sliced version of input, according to the parameters specified.
    '''
    output=[]
    for i in range(0,len(slices)):
        if type(slices[i])==int:
            output.append(array[slices[i]])
        elif len(slices[i])==2:
            for element in array[slices[i][0]:slices[i][1]]:
                output.append(element)
    return output
     

def getData(filepath,outputs,skip_rows=1,transpose_input=True):
    '''
    Function to get specified arrays from a previously generated raw-text csv-file.
    (Just open the Qtegra-.xml in OpenOffice and export as .csv with tabulator separator again.)

    Parameters
    ----------
    filepath : string
        Path of the input csv file
    outputs : string-array
        Abbreviation of the isotope to export. Possible options: 
            '44Ca'
            '83Kr'
            '85Rb'
            '88Sr'
    skip_rows: int, optional
        Integer that is given over to numpy.loadtxt as skiprows. Determines how
        many header rows of the file to read in are ignored. For some QTegra
        output files this seems to be 2, for some 1. Look up if not sure.
        The default is 1.
    transpose_input : boolean, optional
        If the input of a category is given in rows, set this to False.
        If its in columns (as usual), set this to true. The default is True.

    Returns
    -------
    output : array of string/float-arrays
        Array which contains all data specified in outputs, ordered from lowest to highest mass.
        The last array is always the label of the data (e.g. BLK, STD, A1, ...)
    '''
    #Read-in of data (needs to be converted to simple csv list before!)
    daten = np.loadtxt(filepath,dtype='str',skiprows=skip_rows,delimiter='\u0009')
    if transpose_input==True:
        daten=np.matrix.transpose(daten)
    index=daten[0]
    category=daten[1]
    alllabels=daten[2] #shortened version below!
    runs=daten[3]
    sampletype=daten[4]
    dilutionfactor=daten[5]
    comment=daten[6]
    ca44=daten[7]
    kr83=daten[8]
    rb85=daten[9]
    sr88=daten[10]
    
    #Sorting data into array for concentration/intensity, according to elements
    #WARNING: Concentrations for other isotopes than 88Sr are string values!!!
    label=[]
    ca44conc=[]
    d_ca44conc=[]
    ca44int=[]
    d_ca44int=[]
    kr83conc=[]
    d_kr83conc=[]
    kr83int=[]
    d_kr83int=[]
    rb85conc=[]
    d_rb85conc=[]
    rb85int=[]
    d_rb85int=[]
    sr88conc=[]
    d_sr88conc=[]
    sr88int=[]
    d_sr88int=[]
    for i in range(0,int(len(index)/4)):
        label.append(alllabels[4*i])
        ca44int.append(stringToFloat(ca44[4*i]))
        d_ca44int.append(stringToFloat(ca44[4*i+1]))
        ca44conc.append(ca44[4*i+2])
        d_ca44conc.append(ca44[4*i+3])
        kr83int.append(stringToFloat(kr83[4*i]))
        d_kr83int.append(stringToFloat(kr83[4*i+1]))
        kr83conc.append(kr83[4*i+2])
        d_kr83conc.append(kr83[4*i+3])
        rb85int.append(stringToFloat(rb85[4*i]))
        d_rb85int.append(stringToFloat(rb85[4*i+1]))
        rb85conc.append(rb85[4*i+2])
        d_rb85conc.append(rb85[4*i+3])
        sr88int.append(stringToFloat(sr88[4*i]))
        d_sr88int.append(stringToFloat(sr88[4*i+1]))
        sr88conc.append(stringToFloat(sr88[4*i+2]))
        d_sr88conc.append(stringToFloat(sr88[4*i+3]))
    
    #Generate output
    output=[]
    if ('44ca' in outputs) or ('44Ca' in outputs) or ('Ca44' in outputs) or ('ca44' in outputs):
        output.append(ca44int)
        output.append(d_ca44int)
        output.append(ca44conc)
        output.append(d_ca44conc)
    elif ('83kr' in outputs) or ('83Kr' in outputs) or ('Kr83' in outputs) or ('kr83' in outputs):
        output.append(kr83int)
        output.append(d_kr83int)
        output.append(kr83conc)
        output.append(d_kr83conc)
    elif ('85rb' in outputs) or ('85Rb' in outputs) or ('Rb85' in outputs) or ('rb85' in outputs):
        output.append(rb85int)
        output.append(d_rb85int)
        output.append(rb85conc)
        output.append(d_rb85conc)
    elif ('88sr' in outputs) or ('88Sr' in outputs) or ('Sr88' in outputs) or ('sr88' in outputs):
        output.append(sr88int)
        output.append(d_sr88int)
        output.append(sr88conc)
        output.append(d_sr88conc)
    else:
        print('ERROR: Unknown isotope specified')
        return 1
    output.append(label)
    return output


def correctDilution(uncorrected_data,dilution_factors):
    '''
    Function to correct an array of floats with a custom dilution factor each.

    Parameters
    ----------
    uncorrected_data : float-array
        Input array of float values to correct.
    dilution_factors : float/int-array
        Dilution factors array, one factor for EACH element of uncorrected_data!

    Returns
    -------
    corrected_data : float-array
        Array of input data multiplied with the given dilution factors.
    '''
    corrected_data=[]
    for i in range(0,len(uncorrected_data)):
        corrected_data.append(uncorrected_data[i]*dilution_factors[i])
    return corrected_data
        
        
def applyCorrections(data,blank_indices=[],drift_indices=[],calibration_slices=[],calibration_conc=[],d_calibration_conc=[],plots=[],isotope_label=''):
    '''
    Function to apply several corrections on an uncorrected set of counts and
    perform concentration calibration. Divides up into following steps:
        1. Blank correction (if enabled)
        2. Drift correction (if enabled)
        3. Concentration calibration (if enabled)
        4. Zero correction (i.e. blank correction by intercept, if 1. disabled)
        5. Plotting results (optional)
        6. Return of corrected counts and corresponding concentrations

    Parameters
    ----------
    data : array of at least two float-arrays
        Input of the raw data in counts. The first element is expected to be
        the counts themselves, the second array to be their errors. All further
        elements of data are ignored. Supposed to work with the getData()
        output.
    blank_indices : array of ints or int-tuples, optional
        Determines where the blanks are within data. Has to be an array with
        the elements being the indices of blank measurements within the data.
        Here no slice notation, just simple integers accepted!
        If the input is [], a different fit routine is used: Instead of fitting
        a proportional function to calibration vs. concentration, a linear 
        function (with y intercept) is used. The blank correction is then
        determined from the intercept and subtracted after the fit. (= 'zero
        correction'). The default is [].
    drift_indices :  int-array, optional
        Indices of the drift measurements within the dataset. Here no slice
        notation, just simple integers accepted! If left empty, no drift
        correction is performed. The default is [].
    calibration_slices : array of ints or int-tuples, optional
        Determines where the known-concentration measurements are within data.
        Has to have the format used in sliceArray() and specify as many
        elements as there are values within calibration_conc. If left empty, no
        calibration is performed and the last two outputs are canceled.
        The default is [].
    calibration_conc : float/int-array
        The concentrations corresponding to the known-concentration
        measurements specified in calibration_slices.
    d_calibration_conc : float/int-array, optional
        The concentrations corresponding to the known-concentration
        measurements specified in calibration_slices. If none are given, zero
        values are generated. Note that this could lead to a significantly
        lower fit quality because of ODR issues! The default is [].
    plots :  string or string-array, optional
        Whether to display (and save) plots to better assess the quality of the
        fit. Possible options for plots here are:
            'drift'
            'calib'
            'calib_log'
            'calib_counts'
        The default is [].
    isotope_label : string, optional
        A label for the currently corrected isotope. Is used in the plot title
        and for the names of the plot files. The default is ''.

    Returns
    -------
    counts_corr : float-array
        Count array after the drift correction.
    d_counts_corr : float-array
        Error of count array after the drift correction.
    concs_corr : float-array, optional
        Concentration array after the drift correction.
    d_concs_corr : float-array, optional
        Error of concentration array after the drift correction.

    '''
    
    #Read-in
    counts=data[0]
    d_counts=data[1]
    
    #Blank correction
    if blank_indices!=[]:
        print('---- Starting blank correction')
        blank_counts=[]
        d_blank_counts=[]
        for i in blank_indices:
            blank_counts.append(counts[i])
            d_blank_counts.append(d_counts[i])
        mean_blank=np.mean(blank_counts)
        quadsum=0
        for d in d_blank_counts:
            quadsum=quadsum+d**2
        d_mean_blank=np.sqrt(quadsum)/len(blank_indices)
        counts_blank=[]
        d_counts_blank=[]
        for i in range(0,len(counts)):
            counts_blank.append(np.abs(counts[i]-mean_blank))
            d_counts_blank.append(np.sqrt(d_counts[i]**2+d_mean_blank**2))
    if blank_indices==[]:
        counts_blank=counts
        d_counts_blank=d_counts
            
    #Drift correction
    if drift_indices!=[]:
        print('---- Starting drift correction')
        drift_indices=np.sort(drift_indices)
        previousfactor=1
        previouserror=0
        counts_drift=[]
        d_counts_drift=[]
        for k in range(0,drift_indices[0]):
            counts_drift.append(counts_blank[k])
            d_counts_drift.append(d_counts_blank[k])
        for i in range(0,len(drift_indices)):
            currentfactor=counts_blank[drift_indices[i]]/counts_blank[drift_indices[0]]
            currenterror=currentfactor*np.sqrt((d_counts_blank[drift_indices[i]]/counts_blank[drift_indices[i]])**2+(d_counts_blank[drift_indices[0]]/counts_blank[drift_indices[0]])**2)
            for j in range(drift_indices[i-1],drift_indices[i]):
                interpolatedfactor=((currentfactor*(j-drift_indices[i-1])/(drift_indices[i]-drift_indices[i-1])+previousfactor*(drift_indices[i]-j)/(drift_indices[i]-drift_indices[i-1])))
                d_interpolatedfactor=np.sqrt((currenterror*(j-drift_indices[i-1])/(drift_indices[i]-drift_indices[i-1]))**2+(previouserror*(drift_indices[i]-j)/(drift_indices[i]-drift_indices[i-1]))**2)
                counts_drift.append(counts_blank[j]/interpolatedfactor)
                d_counts_drift.append(np.sqrt((d_interpolatedfactor/interpolatedfactor)**2+(d_counts_blank[j]/counts_blank[j])**2)*counts_blank[j]/interpolatedfactor)
            previousfactor=currentfactor
            previouserror=currenterror
        lastdriftcount=counts_blank[j+1]/(counts_blank[drift_indices[len(drift_indices)-1]]/counts_blank[drift_indices[0]])
        counts_drift.append(lastdriftcount)
        d_counts_drift.append(lastdriftcount*np.sqrt((d_counts_blank[j+1]/counts_blank[j+1])**2+(d_counts_blank[drift_indices[len(drift_indices)-1]]/counts_blank[drift_indices[len(drift_indices)-1]])**2+(d_counts_blank[drift_indices[0]]/counts_blank[drift_indices[0]])**2))
        for k in range(drift_indices[len(drift_indices)-1]+1,len(counts_blank)):
            counts_drift.append(counts_blank[k])
            d_counts_drift.append(d_counts_blank[k])
    if drift_indices==[]:
        counts_drift=counts_blank
        d_counts_drift=d_counts_blank
        
    #Concentration calibration
    if calibration_slices!=[]:
        print('---- Starting concentration calibration')
        calibrationcounts=sliceArray(counts_drift,calibration_slices)
        d_calibrationcounts=sliceArray(d_counts_drift,calibration_slices)
        if blank_indices!=[]:
            fitmodel=odr.Model(proportional)
            beta=[1e-6]
        if blank_indices==[]:
            fitmodel=odr.Model(linear)
            beta=[1e-6,0]
        fitdata=odr.RealData(calibrationcounts,calibration_conc,sx=d_calibrationcounts,sy=d_calibration_conc)
        fitodr=odr.ODR(fitdata,fitmodel,beta0=beta)
        fitoutput=fitodr.run()
        print('======= Fit Output =======')
        fitoutput.pprint()
        print('==========================')
        
    #Recalibrate to zero-concentration counts (if activated)
    if blank_indices==[] and calibration_slices!=[]:
        print('---- Starting zero recalibration')
        zero_counts=-fitoutput.beta[1]/fitoutput.beta[0]
        d_zero_counts=zero_counts*np.sqrt((fitoutput.sd_beta[0]/fitoutput.beta[0])**2+(fitoutput.sd_beta[1]/fitoutput.beta[1])**2)
        counts_zero=[]
        d_counts_zero=[]
        for i in range(0,len(counts_drift)):
            counts_zero.append(np.abs(counts_drift[i]-zero_counts))
            d_counts_zero.append(np.sqrt(counts_drift[i]**2+d_zero_counts**2))
            
    #Get new concentrations:
    if calibration_slices!=[]:
        print('---- Getting new concentrations')
        concs_corr=[]
        d_concs_corr=[]
        if blank_indices!=[]:
            for i in range(0,len(counts_drift)):
                newconc=proportional(fitoutput.beta,counts_drift[i])
                concs_corr.append(newconc)
                d_concs_corr.append(np.sqrt((d_counts_drift[i]*fitoutput.beta[0])**2+(counts_drift[i]*fitoutput.sd_beta[0])**2))
        if blank_indices==[]:
            for i in range(0,len(counts_zero)):
                newconc=linear(fitoutput.beta,counts_drift[i])
                concs_corr.append(newconc)
                d_concs_corr.append(np.sqrt((d_counts_drift[i]*fitoutput.beta[0])**2+(fitoutput.sd_beta[0]*counts_drift[i])**2+fitoutput.sd_beta[1]**2))
    
        
    #Plots
    if 'drift' in plots:
        y1=[]
        for i in range(0,len(drift_indices)):
            y1.append(counts_blank[drift_indices[i]])
        plt.figure()
        plt.errorbar(range(0,len(counts_drift)),counts_drift,d_counts_drift,fmt='.',color='green',label='Drift-corrected Counts')
        plt.errorbar(range(0,len(counts_blank)),counts_blank,d_counts_blank,fmt='.',color='blue',label='Uncorrected Counts')
        plt.plot(drift_indices,y1,linestyle='none',marker='o',color='red',label='Drifts')
        plt.axhline(counts_blank[drift_indices[0]], color='grey', linestyle='-',label='Normalized Drift Level')
        plt.title('Drift Correction ('+isotope_label+')')
        plt.xlabel('Index')
        plt.ylabel('Counts')
        plt.legend()
        plt.savefig(output_location+'/drift'+isotope_label+'.jpg',dpi=200)
        
    #Plot calibration
    if 'calib' in plots:
        plt.figure()
        plt.errorbar(calibrationcounts,calibration_conc,d_calibration_conc,d_calibrationcounts,fmt='.',label='Data (baseline-corrected)')
        plt.errorbar(sliceArray(counts,calibration_slices),calibration_conc,d_calibration_conc,sliceArray(d_counts,calibration_slices),fmt='.',label='Data (raw)')
        x=np.linspace(min(calibrationcounts),calibrationcounts[len(calibrationcounts)-1]*1.1,100)
        if blank_indices!=[]:
            plt.plot(x,proportional(fitoutput.beta,x),label='Fit')
        if blank_indices==[]:
            plt.plot(x,linear(fitoutput.beta,x),label='Fit')
        plt.title('Calibration curve ('+isotope_label+')')
        plt.xlabel('Counts')
        plt.ylabel('Concentration [ppb]')
        plt.axhline(y=0, color='grey', linestyle=':')
        plt.legend()
        plt.savefig(output_location+'/calibration'+isotope_label+'.jpg',dpi=200)
        
    #Plot logarithmic calibration
    if 'calib_log' in plots:
        plt.figure()
        plt.errorbar(calibrationcounts,calibration_conc,d_calibration_conc,d_calibrationcounts,fmt='.',label='Data (baseline-corrected)')
        plt.errorbar(sliceArray(counts,calibration_slices),calibration_conc,d_calibration_conc,sliceArray(d_counts,calibration_slices),fmt='.',label='Data (raw)')
        x=np.linspace(min(calibrationcounts),calibrationcounts[len(calibrationcounts)-1]*1.1,100)
        if blank_indices!=[]:
            plt.plot(x,proportional(fitoutput.beta,x),label='Fit')
        if blank_indices==[]:
            plt.plot(x,linear(fitoutput.beta,x),label='Fit')
        plt.title('Calibration curve ('+isotope_label+')')
        plt.xlabel('Counts')
        plt.ylabel('Concentration [ppb]')
        plt.legend()
        plt.xscale('log')
        plt.yscale('log')
        plt.savefig(output_location+'/calibration_log'+isotope_label+'.jpg',dpi=200)
        
    #Plot counts
    if 'calib_counts' in plots:
        plt.figure()
        plt.errorbar(counts,counts_blank,d_counts_blank,d_counts,marker='.',linestyle='none',label='Blank correction only')
        plt.errorbar(counts,counts_drift,d_counts_drift,d_counts,marker='.',linestyle='none',label='Blank+Drift correction')
        plt.plot(np.linspace(0,max(counts)*1.1,2),np.linspace(0,max(counts)*1.1,2),color='grey',label='1:1 line')
        plt.title('Own vs. Qtegra counts ('+isotope_label+')')
        plt.xlabel('QTegra Counts')
        plt.ylabel('Own Counts')
        plt.legend()
        plt.grid()
        plt.savefig(output_location+'/calibration_counts'+isotope_label+'.jpg',dpi=200)
    
    #Outputs
    print('---- Ferdsch.')
    if calibration_slices==[]:
        return counts_drift,d_counts_drift
    elif blank_indices!=[]:
        return counts_drift,d_counts_drift,concs_corr,d_concs_corr
    elif blank_indices==[]:
        return counts_zero,d_counts_zero,concs_corr,d_concs_corr
        
        
def evaluateChem(filepath,isotopes=['88sr'],skip_rows=1,recalibrate=True,blank_indices=[],calibration_conc=[],d_calibration_conc=[],calibration_slices=[], drift_indices=[],
                custom_xs=[],custom_slices=[],custom_labels=[],probe_slices=[],dilution_factors=[],plots='none',table=True, clear_output_destination=True):  
    '''
    Main function to take the data given by getData and plot it/ tabulate it/ 
    calculate column efficiency.
    Uses calibrateData() and correctDrift() in order to perform a custom
    concentration calibration with blank correction and drift correction first.

    Parameters
    ----------
    filepath : string
        Path of input file, is given over to getData().
    isotopes: string/string-array
        Array of string literals identifying the isotope to inspect.
        Possible options are:
            '44ca'/'44Ca'/'Ca44'/'ca44'
            '83kr'/'83Kr'/'Kr83'/'kr83'
            '85rb'/'85Rb'/'Rb85'/'rb85'
            '88sr'/'88Sr'/'Sr88'/'sr88'
        The default is 'sr88'.
    skip_rows: int, optional
        Integer that is given over to numpy.loadtxt as skiprows. Determines how
        many header rows of the file to read in are ignored. For some QTegra
        output files this seems to be 2, for some 1. Look up if not sure.
        Is given over to getData(). The default is 1.
    recalibrate : boolean, optional
        Whether to perform a new calibration with blank correction (True) or
        not (False). The default is True.
    blank_indices : array of ints or int-tuples
        Determines where the blanks are within data. Has to be an array with
        the elements being the indices of blank measurements within the data.
        Here no slice notation, just simple integers accepted!
        If the input is [], a different fit routine is used: Instead of fitting
        a proportional function to calibration vs. concentration, a linear 
        function (with y intercept) is used. The blank correction is then
        determined from the intercept and subtracted after the fit. (= 'zero
        correction'). 
        If the input is [] and the input for drift_indices is [], too, the
        program automatically scans all labels for 'BLK' or 'HNO3 and generates
        the array from that.
        Is given over to applyCorrections(). The default is [].
    calibration_conc : float/int-array
        The concentrations corresponding to the known-concentration
        measurements specified in calibration_slices. If none are given, the
        program is able to find them according to their labels: Each label that
        contains 'Sr', and a number/expression in front of 'ppm'/'ppb'/'ppt'
        is scanned for a usable concentration in number/expression.
        The default is [].
    d_calibration_conc : float/int-array, optional
        The concentrations corresponding to the known-concentration
        measurements specified in calibration_slices. If none are given, ten
        percent rel. are assumed. Note that this could lead to a significantly
        lower fit quality because of ODR issues! The default is [].
    calibration_slices : array of ints or int-tuples, optional
        Determines where the known-concentration measurements are within data.
        Has to have the format used in sliceArray() and specify as many
        elements as there are values within calibration_conc. If none are
        given, the program is able to find them according to their labels: Each
        label that contains 'Sr', and a number/expression in front of 'ppm'/
        'ppb'/'ppt' is taken as a standard then. The default is [].
    drift_indices : int-array, optional
        Indices of the drift measurements within the dataset. Here no slice
        notation, just simple integers accepted! If left empty, no drift
        correction is performed.
        If the input is [] and the input for blank_indices is [], too, the
        program automatically scans all labels for 'Drift' and generates
        the array from that.
        Is given over to applyCorrections(). The default is [].
    custom_xs : array of arrays of arrays of float/int-values, optional
        X-axis data for custom_conc / custom_counts -plots. One overall array
        containing individual arrays for for each different figure containing
        individual arrays for each different data series.
        The default is [].
    custom_slices : array of arrays of arrays of ints or int-tuples, optional
        Slices of the data to plot in custom plot. One overall array containing
        individual arrays for each different figure contaning individual arrays
        for each different data series containing slices (i.e. notation of
        sliceArray()). The default is [].
    custom_labels : array of arrays of {3 strings OR 3 strings + string array},
        opional
        Labels/texts to show in custom plot. One overall array conatining
        individual arrays for each figure in turn containing a title, an
        x-label and a y-label. After the y-label, another array of strings can
        be given, which contains labels for each measurement series in the
        corresponding plot. This will cause a legend to render, obviously.
        The default is [].
    probe_slices : array of ints or int-tuples
        Slices of indices corresponding to probe measurements, in the notation
        used in sliceArray(). Used to cut out anything that is not a probe.
    dilution_factors : float/int-array, optional
        Array of dilution factors to give over to correctDilution() in case a
        dilution is needed. Leave empty if no dilution wished. In that case, do
        NOT specify steps in which dilution correction is needed later on.
        The default is [].
    plots : string or string-array, optional
        Specifies which step (string) or steps (string-array) to plot. Is also
        given over to  calibrateData() and correctDrift().
        Possible options are:
            'calib'
            'calib_log'
            'calib_counts'
            'drift'
            'wash_corr'
            'wash_uncorr'
            'wash_perc'
            'wash_cumul'
            'wash_counts'
            'elution_corr'
            'elution_uncorr'
            'elution_perc'
            'elution_cumul'
            'elution_counts'
            'custom_conc'
            'custom_counts'
        The option 'elution_cumul' includes an efficiency calculation!
        The default is 'none'.
    table : boolean, optional
        If true, prints a table with content corresponding to the other
        parameters and saves it under 'column_auswertung_output/current_table'.
        The default is True.
    clear_output_destination : boolean, optional
        If True, clears all files in output_destination before generating own
        output there. The default is True.

    Returns
    -------
    None.

    '''
    #Remove zombie output
    if clear_output_destination==True:
        for file in listdir(output_location):
            remove(output_location+'/'+file)
    
    alldata=[]
    for isotope in isotopes:
        #Unpacking input
        print('Getting isotope ',isotope,' data for file: ',filepath)
        currentraw=getData(filepath,outputs=isotope,skip_rows=skip_rows,transpose_input=True)
        labels=currentraw[4]
        
        #Generating blank_indices, drift_indices if neither given
        if blank_indices==[] and drift_indices==[] and recalibrate==True:
            print('--Starting blank/drift index fillup')
            for i,l in enumerate(labels):
                if 'blank' in l or 'Blank' in l or 'HNO3' in l or 'HNO_3' in l:
                    blank_indices.append(i)
                    print('---- Found blank at index',i)
                if 'drift' in l or 'Drift' in l or 'DRIFT' in l:
                    drift_indices.append(i)
                    print('---- Found drift at index',i)
                    
        #Generating calibration_slices, calibration_conc if not given
        if calibration_slices==[] and calibration_conc==[] and recalibrate==True:
            math=['0','1','2','3','4','5','6','7','8','9',0,1,2,3,4,5,6,7,8,9,'+','-','*','/','.',',']
            print('--Starting calibration fillup')
            for i,l in enumerate(labels):
                if ('ppt' in l or 'ppb' in l or 'ppm' in l) and ('Sr' in l) and ('Drift' in l or 'drift' in l or 'DRIFT' in l)==False:
                    calibration_slices.append(i)
                    print('---- Found calibration at index',i)
                    print(l)
                    for j in range(0,len(l)):
                        if l[j:j+3]=='ppt' or l[j:j+3]=='ppb' or l[j:j+3]=='ppm':
                            concstring=''
                            for k in range(0,j):
                                if l[j-k-1] in math:
                                    print(l[j-k-1])
                                    concstring=str(l[j-k-1])+concstring
                                    print('ok')
                                if j-k==1:
                                    print('now',concstring)
                                    calibration_conc.append(eval(concstring))
                                    print('------ Concentration: ',concstring)
                                    print('done')
                                    break
                                
        #Generating d_calibration_conc if not given:
        print('--Starting calibration error fillup')
        if d_calibration_conc==[]:
            for c in calibration_conc:
                d_calibration_conc.append(0.1*c)
                    
        #Recalibration / recalibration skipping
        if probe_slices==[]:
            print('-- Generating probe slices')
            probe_slices=[[0,len(currentraw[0])]]
        if recalibrate==True:
            print('-- Starting recalibration for isotope '+isotope)
            allcounts,d_allcounts,allconcs,d_allconcs=applyCorrections(currentraw, blank_indices, drift_indices, calibration_slices, calibration_conc, d_calibration_conc, plots, isotope_label=isotope)
            probecounts=sliceArray(allcounts,probe_slices)
            d_probecounts=sliceArray(d_allcounts,probe_slices)
            probeconc=sliceArray(allconcs,probe_slices)
            d_probeconc=sliceArray(d_allconcs,probe_slices)
        if recalibrate==False:
            allcounts=currentraw[0]
            d_allcounts=currentraw[1]
            allconcs=currentraw[2]
            d_allconcs=currentraw[3]
            probecounts=sliceArray(currentraw[0],probe_slices)
            d_probecounts=sliceArray(currentraw[1],probe_slices)
            for k in range(0,len(currentraw[2])):
                currentraw[2][k]=stringToFloat(currentraw[2][k])
                currentraw[3][k]=stringToFloat(currentraw[3][k])
            probeconc=sliceArray(currentraw[2],probe_slices)
            d_probeconc=sliceArray(currentraw[3],probe_slices)
        probelabels=sliceArray(currentraw[4],probe_slices)
        
        #Collecting data from multiple isotopes
        alldata.append([probecounts,d_probecounts,probeconc,d_probeconc])
        '''alldata currently unused, just exists in case someone wants to
        program cross-isotope plots. probelabels neither.'''
                                           
        #Dilution correction
        if dilution_factors!=[]:
            sr88corr=correctDilution(probeconc, dilution_factors)
            d_sr88corr=correctDilution(d_probeconc, dilution_factors)
        
        #Printing/saving table
        if table==True:
            table1=np.matrix.transpose(np.array([range(0,len(currentraw[4])),np.array(currentraw[4]),np.array(currentraw[0]),np.array(currentraw[1]),np.array(currentraw[2]),np.array(currentraw[3])]))
            header1=['Index','Label','QTegra Counts','(error)','QTegra Concentration','(error)']
            if recalibrate==True:
                table1=np.matrix.transpose(np.array([range(0,len(currentraw[4])),np.array(currentraw[4]),np.array(currentraw[0]),np.array(currentraw[1]),np.array(currentraw[2]),np.array(currentraw[3]),allcounts,d_allcounts,allconcs,d_allconcs]))
                header1=['Index','Label','QTegra Counts','(error)','QTegra Concentration','(error)','Own Counts.','(error)','Own Conc.','(error)']
            if dilution_factors!=[]:
                table1=np.matrix.transpose(np.array([range(0,len(currentraw[4])),np.array(currentraw[4]),np.array(currentraw[0]),np.array(currentraw[1]),np.array(currentraw[2]),np.array(currentraw[3]),sr88corr,d_sr88corr]))
                header1=['Index','Label','QTegra Counts','(error)','QTegra Concentration','(error)','Corr. Conc.','(error)']
            print(tabulate(table1,headers=header1))
            headers=''
            for s in header1:
                headers=headers+s+'\u0009'
            footers=np.array2string(np.array(['Isotope:',isotope,'Calibration data:',calibration_slices,calibration_conc,d_calibration_conc]))
            np.savetxt(output_location+'/column_auswertung_output_'+isotope+'.csv',table1,delimiter='\u0009',fmt="%s",header=headers,footer=footers)
            
        #Plotz
        if 'custom_counts' in plots or 'custom_conc' in plots:
            symbols=['v','o','<','>','*','s']
            while len(symbols)<len(custom_xs[0]):
                symbols.append(symbols)
            
        if 'conc_corr' in plots and dilution_factors!=[]:
            for s in range(0,len(custom_slices)):
                plt.figure()
                for j in range(0,len(custom_xs[s])):
                    if len(custom_labels[s])==4:
                        plt.errorbar(custom_xs[s][j],sliceArray(sr88corr,custom_slices[s][j]),sliceArray(d_sr88corr,custom_slices[s][j]),label=custom_labels[s][3][j],fmt=symbols[j],markersize=9)
                        plt.legend()
                    else:
                        plt.errorbar(custom_xs[s][j],sliceArray(sr88corr,custom_slices[s][j]),sliceArray(d_sr88corr,custom_slices[s][j]),fmt='.')
                plt.title(custom_labels[s][0]+' (Dilution-corrected counts, '+isotope+')')
                plt.xlabel(custom_labels[s][1])
                plt.ylabel('counts in a.u. '+custom_labels[s][2])
                plt.savefig(output_location+'/'+custom_labels[s][0]+'_'+isotope+'_conc_corr.jpg',dpi=200)
        
        if 'custom_counts' in plots:
            for s in range(0,len(custom_slices)):
                plt.figure()
                for j in range(0,len(custom_xs[s])):
                    if len(custom_labels[s])==4:
                        plt.errorbar(custom_xs[s][j],sliceArray(probecounts,custom_slices[s][j]),sliceArray(d_probecounts,custom_slices[s][j]),label=custom_labels[s][3][j],fmt=symbols[j],markersize=9)
                        plt.legend()
                    else:
                        plt.errorbar(custom_xs[s][j],sliceArray(probecounts,custom_slices[s][j]),sliceArray(d_probecounts,custom_slices[s][j]),fmt='.')
                plt.title(custom_labels[s][0]+' (Counts, '+isotope+')')
                plt.xlabel(custom_labels[s][1])
                plt.ylabel('counts in a.u. '+custom_labels[s][2])
                plt.savefig(output_location+'/'+custom_labels[s][0]+'_'+isotope+'_counts.jpg',dpi=200)
                
        if 'custom_conc' in plots:
            for s in range(0,len(custom_slices)):
                plt.figure()
                for j in range(0,len(custom_xs[s])):
                    if len(custom_labels[s])==4:
                        print(custom_xs[s][j])
                        print(sliceArray(probeconc,custom_slices[s][j]))
                        print(sliceArray(d_probeconc,custom_slices[s][j]))
                        print(custom_labels[s][3][j])
                        plt.errorbar(custom_xs[s][j],sliceArray(probeconc,custom_slices[s][j]),sliceArray(d_probeconc,custom_slices[s][j]),label=custom_labels[s][3][j],fmt=symbols[j],markersize=9)
                        plt.legend()
                    else:
                        print(custom_xs[s][j])
                        print(sliceArray(probeconc,custom_slices[s][j]))
                        print(sliceArray(d_probeconc,custom_slices[s][j]))
                        plt.errorbar(custom_xs[s][j],sliceArray(probeconc,custom_slices[s][j]),sliceArray(d_probeconc,custom_slices[s][j]),fmt='.')
                plt.title(custom_labels[s][0]+' (Conc., '+isotope+')')
                plt.xlabel(custom_labels[s][1])
                plt.ylabel('concentration in ppb '+custom_labels[s][2])
                plt.savefig(output_location+'/'+custom_labels[s][0]+'_'+isotope+'_concs.jpg',dpi=200)
                        
        #Save current Python script
        copyfile(current_location,output_location+'/RESTREX_Chemistry.py')

        
        
###############################################################################



file='Desktop/RESTREX_Starterkit/Example_Data_Chemistry/20220223_11365_SR1d-SR1j.csv'



#Basic setup
evaluateChem(filepath=file, isotopes=['Sr88'], skip_rows=2, recalibrate=True,
             blank_indices=[], calibration_conc=[], d_calibration_conc=[], 
             calibration_slices=[], drift_indices=[], 
             plots=['calib','calib_log','drift'],
             table=True,clear_output_destination=False)


#Advanced setup with fancy custom plot
custom_xs=[[[1/2*5.67,1/4*5.67,1/8*5.67,1/16*5.67,1/32*5.67,0],[1/2*5.67,1/4*5.67,1/8*5.67,1/16*5.67,1/32*5.67,0],[1/2*5.67,1/4*5.67,1/8*5.67,1/16*5.67,1/32*5.67,0],[1/2*5.67,1/4*5.67,1/8*5.67,1/16*5.67,1/32*5.67,0],[1/2*5.67,1/4*5.67,1/8*5.67,1/16*5.67,1/32*5.67,0]]]
custom_slices=[[[[12,18]],[[29,35]],[[38,44]],[[47,53]],[[63,69]]]]
custom_labels=[['Double elution resolved over steps & mass',
                'Equivalent coral mass [mg]','',
                ['Elution 2', 'Sample 1', 'Wash 1','Sample 2','Wash 2']]]
evaluateChem(filepath=file, isotopes=['Sr88'], skip_rows=2, recalibrate=True,
             blank_indices=[], calibration_conc=[], d_calibration_conc=[], 
             calibration_slices=[], drift_indices=[], 
             custom_xs=custom_xs, custom_slices=custom_slices,
             custom_labels=custom_labels, 
             plots=['custom_conc'],
             table=True,clear_output_destination=False)


#Pro setup with control over all parameters
calibconc=[5,50,500,5000]
calibconcerrors=[0.5,5,50,500]
calibslices=[3,5,8,9]
blanks=[1,10,18,26,35,44,53,61,69]
drifts=[11,19,27,36,45,54,62,70]
custom_xs=[[[-3,-2,-1,1,2,3,4,5,6]],[[0.5,1.0,1.5,2.0,2.5,3.0]]]
custom_slices=[[[28,37,46,[55,61]]],[[[20,26]]]]
custom_labels=[['Wash of 6.12mg sample','ml wash',''],
               ['Elution of 6.12mg sample','ml elution','']]
evaluateChem(filepath=file, isotopes=['Ca44'], skip_rows=2, recalibrate=True,
             blank_indices=blanks, calibration_conc=calibconc,
             d_calibration_conc=calibconcerrors,
             calibration_slices=calibslices, drift_indices=drifts,
             custom_xs=custom_xs, custom_slices=custom_slices,
             custom_labels=custom_labels,
             plots=['calib_log','custom_conc'],
             table=True, clear_output_destination=False)


#Don't forget to change the directory fields in the beginning when moving the file!