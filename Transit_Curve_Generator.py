import time
import numpy as np

#Curve fitting package
import lmfit as lmfit
from lmfit import Model, fit_report

#Important distinction
from lmfit import Parameters
from lmfit import Parameter

#Graph rendering package
import matplotlib.pyplot as matplot

#Alternative to text boxes that allows right alignment
from matplotlib.offsetbox import AnchoredText

#Transit curve calculation package
import batman

#Profiling library
import cProfile

#WARNING : Only works if OpenMP is enabled
#Do not set this value above 1 if OpenMP is not available
BatmansThreads = 1

#Note:
#periastron  is same as periapsis
#It just refers to the periapsis of objects orbiting stars

#If supplied, will be used as initial fitting parameters
#First array element is value, second is uncertainty
Priors = [
    [0.96, 999],  #t0
    [3.14159265, 999],  #per
    [0.118, 999],  #rp
    [8.0, 999],  #a
    [83.7, 999],  #inc
    [0.0, 999],  #ecc
    [1.0, 999],  #w
    [-0.9, 999],  #u1
    [-0.9, 999],  #u2
    [1, 999]  #ScalingMultiplier
]
ParamNames = 't0', 'per', 'rp', 'a', 'inc', 'ecc', 'w', 'u1', 'u2'
PriorsDict = dict()
#Zip function 'links' objects together like an array of arrays
for Name, Prior in zip(ParamNames, Priors):
    PriorsDict[Name] = Prior

PolynomialOrder = 3

#NOTE: This program assumes only 2 limb-darkening values, more are possible with modification

#Visual modifier only
DataPointRenderSize = 1
#Suggested Values:
#1, apears as point cloud, easy to see function line
#2, larger than error bar, easy to see points when error value is high

StartTimes = []
EndTimes = []
ProgramStartTime = time.time()


def CheckTime(ID, IsStart):
    #used to check how long a section of code took to run
    #Purely for debug

    #Can now support multiple checked blocks
    #Specify "ID" (int) of a block when starting and ending recording to separate them

    global StartTimes
    global EndTimes
    if (IsStart):
        if (len(StartTimes) <= ID):
            StartTimes.append(time.time())
        else:
            StartTimes[ID] = time.time()
    else:
        #Is end
        print("Block Time ID :", str(ID), ": Time :",
              str(time.time() - StartTimes[ID]))
        EndTimes.append(time.time())


def ReplaceZerosInArrayWithLowestValue(DataArray):

    #Replaces zeros in an array with the lowest value present in the array
    #This is normally used on an array of error values, as values of '0' leed to calculation errors and are improbable
    #EX: (Input - Output) / Error = Divide By Zero Error

    FixedArray = DataArray

    LowestArrayValue = GetArrayBounds(DataArray)[0]

    for Count in range(len(FixedArray)):
        if (FixedArray[Count] == 0):
            #print("'0' Array Value At Index : " + str(Count) + " : Replacing With Lowest Array Value")

            FixedArray[Count] = LowestArrayValue

    return (FixedArray)


def CalculateChiSqr(DataX, DataY, Params, DataERROR):

    TransiteParams = ConvertFitParametersToTransitParameters(Params)

    global BatmansThreads
    m = batman.TransitModel(TransiteParams, DataX,nthreads = BatmansThreads)
    flux = m.light_curve(TransiteParams)
    flux = ApplyPolyMultiplier(DataX,flux,Params)

    CheckedOptimizedChiSqr = 0

    DataIncludedErrorBars = GetArrayIsNotNone(DataERROR)

    #sumation of ((data_i - model_i) / uncertainty_i)^2
    if (DataIncludedErrorBars):
        #Previousely was : CheckedOptimizedChiSqr = (((((DataY-flux)**2))/(DataERROR))).sum()
        CheckedOptimizedChiSqr = (((DataY - flux) / DataERROR)**2).sum()
    else:
        CheckedOptimizedChiSqr = (((DataY - flux))**2).sum()
    


    #Removing this section decreases run time, because it increases the ammount of cycles nelder runs when it is included.
    #Previouse comments about this section increasing run time because it was slow were wrong, this section is not to blame. I was wrong.
    ParameterValues = Params.valuesdict()
    global PriorsDict

    for ParamName in PriorsDict.keys():
        if PriorsDict[ParamName][1] is not None:
            CheckedOptimizedChiSqr += ((PriorsDict[ParamName][0] - ParameterValues[ParamName]) /PriorsDict[ParamName][1])**2
    
    

    return (CheckedOptimizedChiSqr)


def ReturnChiModiferOfParameterPrior(Param, Prior, PriorERROR):
    #print(((Param-Prior)/PriorERROR)**2)
    return (((Param - Prior) / PriorERROR)**2)


def Clamp(Value, Min, Max):
    #Should replace with external package

    if (Value < Min):
        return (Min)
    if (Value > Max):
        return (Max)
    return (Value)


def GetArrayBounds(DataArray):
    DataArray = np.array(DataArray, copy=False)
    return ([DataArray.min(), DataArray.max()])


NelderEvaluations = 0


'''
def CustomChiSqrInputFunction(Params, DataX, DataY, DataERROR, Priors):

    ParameterValues = Params.valuesdict()
    
    ''
    FittingTransityFunctionParams = batman.TransitParams()

    #Maybe add function to copy params from parameter values dict, 
    FittingTransityFunctionParams.t0 = ParameterValues["t0"]                        #time of inferior conjunction
    FittingTransityFunctionParams.per = ParameterValues["per"]                       #orbital period
    FittingTransityFunctionParams.rp = ParameterValues["rp"]                       #planet radius (in units of stellar radii)
    FittingTransityFunctionParams.a = ParameterValues["a"]                        #semi-major axis (in units of stellar radii)
    FittingTransityFunctionParams.inc = ParameterValues["inc"]                      #orbital inclination (in degrees)
    FittingTransityFunctionParams.ecc = ParameterValues["ecc"]                       #eccentricity
    FittingTransityFunctionParams.w = ParameterValues["w"]                        #longitude of periastron (in degrees)
    FittingTransityFunctionParams.limb_dark = "quadratic"        #limb darkening model
    FittingTransityFunctionParams.u = [ParameterValues["u1"], ParameterValues["u2"]]      #limb darkening coefficients [u1, u2]
    ''


    global NelderEvaluations
    NelderEvaluations += 1

    #Will try to minimize returned value
    return (CalculateChiSqr(DataX, DataY, Params, DataERROR))
'''

LBMIterations = 0
def LmfitInputFunction(Params, DataX, DataY, DataERROR, Priors, IsNelder):

    if(IsNelder):
        global NelderEvaluations
        NelderEvaluations+=1
    else:
        global LBMIterations
        LBMIterations+=1

    ParameterValues = Params.valuesdict()


    FittingTransityFunctionParams = ConvertFitParametersToTransitParameters(ParameterValues)


    global BatmansThreads
    Flux = batman.TransitModel(FittingTransityFunctionParams, DataX,nthreads = BatmansThreads).light_curve(FittingTransityFunctionParams)
    Flux = ApplyPolyMultiplier(DataX, Flux, ParameterValues)


    ReturnChiArray = abs(DataY - Flux)


    if (not DataERROR is None):
        ReturnChiArray /= DataERROR

    global PriorsDict
    FoundValidPrior = False
    ModifiedPriorValues = []
    for ParamName in PriorsDict.keys():
        if PriorsDict[ParamName][1] is not None:
            ModifiedPriorValues.append(abs((PriorsDict[ParamName][0] - ParameterValues[ParamName]) / PriorsDict[ParamName][1]))

            FoundValidPrior = True
            #print(str((abs((PriorsDict[ParamName][0] - ParameterValues[ParamName])/PriorsDict[ParamName][1]))))

    FoundValidPrior = False

    if FoundValidPrior:
        ReturnChiArray = np.concatenate((ReturnChiArray, ModifiedPriorValues), axis=0)

    #Debug logging
    #If initial params are '0' minor changes will not affect the result enough for proper fitting
    #Check 'a' parameter is set properly, do not leave as initialized value
    #print(Params.valuesdict(),str(DataY-Flux))

    for i in ReturnChiArray:
        if(i < 0):
            print(i)

    #if(IsNelder):
    #    return (np.sum(ReturnChiArray))
    #else:
    #    return (ReturnChiArray)


    #Draws graph after each fit, only for debugging, very slow
    global DrawProgressiveFitGraph
    if(DrawProgressiveFitGraph):
        ContinouseDrawGraph(DataX, DataY, ParameterValues)


    return (ReturnChiArray)


def OptimizeFunctionParameters(DataX, DataY, DataERROR, Priors, UseLBM, StartingParameters):

    Bounds = GetArrayBounds(DataX)
    MinX = Bounds[0]
    MaxX = Bounds[1]

    global PolynomialOrder

    InputParams = lmfit.Parameters()

    if (StartingParameters is not None):
        #Lmfit version
        InputParams.add("t0",value=StartingParameters["t0"],min=MinX,max=MaxX)  #Max?
        InputParams.add("per",value=StartingParameters["per"],min=0.0,max=MaxX-MinX)
        InputParams.add("rp", value=StartingParameters["rp"], min=0, max=1.0)
        InputParams.add("a", value=StartingParameters["a"], min=1.0,max=90)  #In ratio of stelar radii
        InputParams.add("inc",value=StartingParameters["inc"],min=60,max=90) # Degrees from "top" of star to orbit plane
        InputParams.add("ecc",value=StartingParameters["ecc"],min=0.0,max=1.0)
        InputParams.add("w",value=StartingParameters["w"],min=0.0,max=360.0) # from 0-360 or -180 to 180. look at batman docs
        InputParams.add("u1",value=StartingParameters["u1"],min=0.0,max=1.0)
        InputParams.add("u2",value=StartingParameters["u2"],min=0.0,max=1.0)
                   
    elif (Priors is not None):
        #Lmfit version
        InputParams.add("t0", value=Priors[0][0], min=MinX,max=MaxX)  #Max?
        InputParams.add("per", value=Priors[1][0], min=0.0, max=MaxX)
        InputParams.add("rp", value=Priors[2][0], min=0, max=10.0)
        InputParams.add("a", value=Priors[3][0], min=1.0,max=90)  #What should Max Bound be?
        InputParams.add("inc", value=Priors[4][0], min=60, max=90)
        InputParams.add("ecc", value=Priors[5][0], min=0.0, max=1.0)
        InputParams.add("w", value=Priors[6][0], min=0.0, max=360.0)
        InputParams.add("u1", value=Priors[7][0], min=-1.0, max=1.0)
        InputParams.add("u2", value=Priors[8][0], min=-1.0, max=1.0)

    else:
        #Backup 
        InputParams.add("t0", value=MaxX/2, min=MinX, max=MaxX)  #Max?
        InputParams.add("per", value=MaxX/2, min=0.0, max=MaxX)
        InputParams.add("rp", value=5, min=0, max=10.0)
        InputParams.add("a", value=40, min=1.0,max=90)  #What should Max Bound be?
        InputParams.add("inc", value=70, min=60, max=90)
        InputParams.add("ecc", value=0.5, min=0.0, max=1.0)
        InputParams.add("w", value=40, min=0.0, max=360.0)
        InputParams.add("u1", value=0.5, min=-1.0, max=1.0)
        InputParams.add("u2", value=0.5, min=-1.0, max=1.0)

    if(PolynomialOrder > 0):
        for PolyIndex in range(0,PolynomialOrder):
            if(PolyIndex == PolynomialOrder-1):
                InputParams.add(("PolyVal" + str(PolyIndex)), value=1, min=-1000, max=1000, vary = True)
            else:
                InputParams.add(("PolyVal" + str(PolyIndex)), value=0, min=-1000, max=1000, vary = True)

    #InputParams.add(("PolyVal1"), value=0, min=0, max=0.01, vary = False)
    #InputParams.add(("PolyVal2"), value=1, min=1, max=1.01, vary = False)

    #print(InputParams)

    #InputParams.add("Test_1", value=0, min=-1000, max=1000)
    #InputParams.add("Test_2", value=0, min=-1000, max=1000)
    #InputParams.add("Test_3", value=0, min=-1000, max=1000)
    #InputParams.add("Test_4", value=0, min=-1000, max=1000)
    #InputParams.add("Test_5", value=0, min=-1000, max=1000)
    #InputParams.add("Test_6", value=0, min=-1000, max=1000)
    #InputParams.add("Test_7", value=0, min=-1000, max=1000)
    #InputParams.add("Test_8", value=0, min=-1000, max=1000)

    OptimizedFunctionToReturn = None

    if (not UseLBM):
        OptimizedFunctionToReturn = lmfit.minimize(
            LmfitInputFunction,
            InputParams,
            args=(DataX, DataY, DataERROR, Priors, True),
            method="nelder",
            calc_covar=True,
            max_nfev=None,
            nan_policy="raise")
    else:

        #   lmfit.minimize(
        #Function
        #Params
        #Arguments (will be passed to function, not modified by or visible to fitting method)
        #Method, defaults to 'Levenberg-Marquardt' referenced by 'leastsq', assigning it specifically should not be necesary [Different from 'least_squares']
        #Wether to calculate uncertainties (if fitting method supports it), should default to true
        #Maximum evaluations, should default to none
        #Weights not included, are being processed in lmfit function instead of solver to allow for inclusion of prior weights
        #Raise issue if NAN found, this shouldn't happen, but I would like to know about it if it does because it means something isn't being processed correctly.

        OptimizedFunctionToReturn = lmfit.minimize(
            LmfitInputFunction,
            InputParams,
            args=(DataX, DataY, DataERROR, Priors, False),
            method="leastsq",
            calc_covar=True,
            max_nfev=None,
            nan_policy="raise")



        #Weight Implementation According To Documentation:
        '''
        weights
        numpy.ndarray (or None) of weighting values to be used in fit.
        If not None, it will be used as a multiplicative factor of the residual array, so that weights*(data - fit) is minimized in the least-squares sense.
        '''

    return (OptimizedFunctionToReturn)


#Main function
def RunOptimizationOnDataInputFile(Priors):

    
    '''
    #Debug Testing
    TestParameters = batman.TransitParams()
    
    UseGoodValues = True

    if(UseGoodValues != True):
        #Bad
        TestParameters.t0 = 0.97779202
        TestParameters.per = 3.14110619
        TestParameters.rp = 0.10529006
        TestParameters.a =0.10529006
        TestParameters.inc = 0.10529006
        TestParameters.ecc = 0.10529006
        TestParameters.w = 0.10529006
        TestParameters.limb_dark = "quadratic"
        TestParameters.u = [-0.03773462, 0.56399941]
    else:
        #Good
        TestParameters.t0 = 0.9771297617522302
        TestParameters.per = 3.1411255188855454
        TestParameters.rp = 0.10448573735254096
        TestParameters.a =1.4315349098124923
        TestParameters.inc = 63.43883728462587
        TestParameters.ecc = 0.9497616962968043
        TestParameters.w = 2.050529211138452
        TestParameters.limb_dark = "quadratic"
        TestParameters.u = [-0.05394222113055169, -0.05394222113055169]

    SamplePoints = np.linspace(0.0, 26.996528, 10000)
    print(SamplePoints)
    global BatmansThreads
    m = batman.TransitModel(TestParameters, SamplePoints,nthreads = BatmansThreads)
    flux = m.light_curve(TestParameters)
    matplot.plot(SamplePoints, flux, "-", label="Optimized Function")

    matplot.show()
    '''



    #np.array used in this starting section is likely ineficient, but not the main time sink for this program

    print("Running. Please wait...\n")

    DataPoints = np.array([[0, 0, 0]])

    #Clear initialized array value, because I don't know how to initialize empty array with bounds
    DataPoints = np.delete(DataPoints, 0, 0)

    FileName = "Data"
    FileType = "txt"

    # Assume CSV format with no header row:
    datafile = np.loadtxt(FileName + "." + FileType, delimiter=',')
    if datafile.shape[1] == 2:
        DataIncludedErrorBars = False
        DataX = datafile[:, 0]
        DataY = datafile[:, 1]
        DataERROR = -np.ones(DataX.shape)
    elif datafile.shape[1] == 3:
        DataIncludedErrorBars = True
        DataX = datafile[:, 0]
        DataY = datafile[:, 1]
        DataERROR = datafile[:, 2]

    #Get data bounds

    #Used in some parameter bounds
    Bounds = GetArrayBounds(DataX)
    MinX = Bounds[0]
    MaxX = Bounds[1]

    #Not currently used
    Bounds = GetArrayBounds(DataY)
    MinY = Bounds[0]
    MaxY = Bounds[1]

    NumberOfDataPoints = len(DataX)

    CheckTime(0, True)
    # - 5.2s
    #First optimization attempt, used to get error values

    #Running this first iteration with 'nelder' generates significantly better initial results than 'leastsqr'
    #Running using 'leastsqr' results in a badly fit graph that is then used to remove outlier data points. This bad graph leads to good data points being thrown out, and the final graph is bad because of it.

    OptimizedFunction = OptimizeFunctionParameters(DataX, DataY, None, Priors, False, None)
    CheckTime(0, False)
    global NelderEvaluations
    #print("Nelder Evaluations :",str(NelderEvaluations),": Data Points :",len(DataX),": Dif :",str(NelderEvaluations/len(DataX)))

    #Extract parameters used
    OptimizedBatmanParams = ExtractTransitParametersFromFittedFunction(OptimizedFunction.params)
    OptimizedParamsDictionary = OptimizedFunction.params
    '''
    #Debugging
    print("HERE")
    print(fit_report(OptimizedFunction))
    time.sleep(999)
    print("HERE")
    '''

    #Generate function based on extracted parameters
    global BatmansThreads
    FirstOptimizedFunction = batman.TransitModel(OptimizedBatmanParams, DataX,nthreads = BatmansThreads).light_curve(OptimizedBatmanParams)

    #Calculate error from diference between first attempt created function, and given values
    if (not DataIncludedErrorBars):
        #I think "abs" should not be used here, the square of the values is being used instead. Not sure why abs affects the result in this case, but it does.
        DataERROR = (DataY * 0 + np.std((DataY - FirstOptimizedFunction)))
        #CheckChiSqr function expects data error as an array, this allows compatibilty with lmfit.fit instead of lmfit.minimize

        #Debug logging
        #print(np.std((DataY-FirstOptimizedFunction)))


    #Disable this to see if (too many)/(good) data points are being removed after the first fit.
    #If this is happeneing the rpiros are liley too restrictive or too far from the actual values.
    RemoveOutliers = True

    if(RemoveOutliers):
        #Remove Outlier values
        NewDataValues = RemoveOutliersFromDataSet(DataX, DataY, OptimizedParamsDictionary)

        DataX = NewDataValues[0]
        DataY = NewDataValues[1]
        NumberOfDataPoints = NewDataValues[2]
        IndexesRemoved = NewDataValues[3]
        if (len(IndexesRemoved) > 0):
            DataERROR = np.delete(DataERROR, IndexesRemoved)

        #Recalculate error
        #If data did not include errors, and outliers have just been removed
        if (not DataIncludedErrorBars and (len(IndexesRemoved) > 0)):

            #Have to remove removed values from returned light values, because can't calculate std of diference betweeen arrays, when those arrays are of diferent lengths
            UpdatedLightValues = np.delete(FirstOptimizedFunction,IndexesRemoved)

            #Recalcualte error values with the outlier values removed
            DataERROR = (DataY * 0 + np.std((DataY - UpdatedLightValues)))



    #Run second time, this time having removed outliers and calculated error values if they were not provided
    CheckTime(1, True)
    SecondOptimizedFunction = OptimizeFunctionParameters(DataX, DataY, DataERROR, Priors, True, OptimizedParamsDictionary)
    CheckTime(1, False)

    DictionaryParams = SecondOptimizedFunction.params
    BatmanParams = ExtractTransitParametersFromFittedFunction(DictionaryParams)

    #From this point forward, "DataIncludedErrorBars" will no longer be used to decide if Error values need to be calculated
    #It just means wether to use the current Error np.array for Chi calcualtions and wether to render points with error bars
    DataIncludedErrorBars = True

    #Debug Fit Report
    #print(fit_report(FinalOptimizedFunction))
    print(DictionaryParams)
    print("\n")

    #Display points with error bars
    if (DataIncludedErrorBars):
        #Disabled for now
        #Concerned the 'yerr' value being set to the DataERROR value is not an acurate representation of the points error

        #matplot.errorbar(DataX, DataY, yerr = DataERROR, fmt ="o", markersize = DataPointRenderSize)
        matplot.errorbar(DataX, DataY, fmt="o", markersize=DataPointRenderSize)
    else:
        matplot.errorbar(DataX, DataY, fmt="o", markersize=DataPointRenderSize)

    print("-OPTIMIZED PARAMETERS-")
    parameter_names = 't0', 'per', 'rp', 'a', 'inc', 'ecc', 'w', 'u1', 'u2'
    for name in parameter_names:
        print(name + ' : ' + str(DictionaryParams[name].value))
    print("\n")

    print("-OPTIMIZED PARAMETER UNCERTAINTY VALUES-")
    for name in parameter_names:
        print(name + ' : ' + str(DictionaryParams[name].stderr))
    print("\n")

    CheckedOptimizedChiSqr = CalculateChiSqr(DataX, DataY, DictionaryParams, DataERROR)
    #print("\n\n",MinX,"   ",MaxX)
    #Rendering only, uses more sample points than input x-values
    SamplePoints = np.linspace(MinX, MaxX, 10000)
    m = batman.TransitModel(BatmanParams, SamplePoints,nthreads = BatmansThreads)
    flux = m.light_curve(BatmanParams)
    flux = ApplyPolyMultiplier(SamplePoints, flux, DictionaryParams)
    matplot.plot(SamplePoints, flux, "-", label="Optimized Function")

    '''
    #Debug Logging START

    StringData = ""

    DebugFlux = batman.TransitModel(BatmanParams,DataX,nthreads = BatmansThreads).light_curve(BatmanParams)

    for i in range(len(DataX)):
        StringData += (str(DebugFlux[i]) + "\n")
    
    print("-----")

    print(str(np.std((DataY-DebugFlux))))

    print("-----")
    
    with open("Output.txt", "w") as File:
        Lines = File.writelines(StringData)
    '''

    print("\n--- Checked Chi Sqr ---")
    print("ChiSqr : " + str(CheckedOptimizedChiSqr))
    print("Number Of Data Points : " + str(NumberOfDataPoints))
    #The value below should be close to '1'
    print("ChiSqr / # Data Points : " + str(CheckedOptimizedChiSqr / NumberOfDataPoints))

    #Debug Logging END

    #Fixed "χ2" rendering issue
    BestChi = "Optimized χ2 : " + str(round(CheckedOptimizedChiSqr, 2))
    ReducedChi = "Reduced χ2 : " + str(round(CheckedOptimizedChiSqr / NumberOfDataPoints, 5))

    #Text box setup
    ChiAnchoredTextBox = AnchoredText((BestChi + "\n" + ReducedChi), loc=4, pad=0.5)
    matplot.setp(ChiAnchoredTextBox.patch, facecolor="Orange", alpha=0.5)
    matplot.gca().add_artist(ChiAnchoredTextBox)
    matplot.legend(loc=2, borderaxespad=0)


    #Enable this to see a graph of the values the priors would generate
    #All priors must be assigned a value for this to work
    #If this graph is too far from the data and has low uncertainty, this can skew the fit, or result in a high chisqr value even if the fit is good, because it does not align with the priors.
    ShowPriorGraph = False
    if(ShowPriorGraph):
        PriorParams = batman.TransitParams()

        PriorParams.t0 = PriorsDict["t0"][0]  #time of inferior conjunction
        PriorParams.per = PriorsDict["per"][0]  #orbital period
        PriorParams.rp = PriorsDict["rp"][0]  #planet radius (in units of stellar radii)
        PriorParams.a = PriorsDict["a"][0]  #semi-major axis (in units of stellar radii)
        PriorParams.inc = PriorsDict["inc"][0]  #orbital inclination (in degrees)
        PriorParams.ecc = PriorsDict["ecc"][0]  #eccentricity
        PriorParams.w = PriorsDict["w"][0]  #longitude of periastron (in degrees)
        PriorParams.limb_dark = "quadratic"  #limb darkening model
        PriorParams.u = [PriorsDict["u1"][0], PriorsDict["u2"][0]]  #limb darkening coefficients [u1, u2]


        SamplePoints = np.linspace(MinX, MaxX, 10000)
        Flux = batman.TransitModel(PriorParams, SamplePoints,nthreads = BatmansThreads).light_curve(PriorParams)
        matplot.plot(SamplePoints, Flux, "-", label="Prior Graph")

    EndTimeRecording()

    print("\nCompleted")

    global NelderEvaluations
    print("Nelder Iterations :",NelderEvaluations)
    NelderEvaluations = 0

    global LBMIterations
    print("LBM Iterations :",LBMIterations)
    LBMIterations = 0

    #Display plotted data
    #Code after this function is called will not be run untill the graph is closed (this behavior can be changed)
    global TestAvergageTimeMode
    if(not TestAvergageTimeMode):
        matplot.show()


def RemoveOutliersFromDataSet(DataX, DataY, Parameters):

    print(Parameters['t0'])

    TestMode = False
    
    #-Only valid if 'TestMode' is active

    #will not halt further execution of the program, instead overlaying the scatter points or heat map underneath the final graph
    OverlayMode = False

    #Show limits values are allowed between
    HighlightBoundsMode = True

    #-

    NewDataX = DataX
    NewDataY = DataY
    NewNumberOfDataPoints = -1

    BatmanParams = ExtractTransitParametersFromFittedFunction(Parameters)

    global BatmansThreads
    TransitModel = batman.TransitModel(BatmanParams, DataX,nthreads = BatmansThreads)
    LightCurve = TransitModel.light_curve(BatmanParams)
    LightCurve = ApplyPolyMultiplier(DataX, LightCurve, Parameters)

    StandardDeviation = (np.std(DataY - LightCurve))

    Differences = abs((DataY - LightCurve))

    MeanDifference = Differences.sum() / len(Differences)

    DiferenceBounds = GetArrayBounds(Differences)

    DifferenceMin = DiferenceBounds[0]
    DifferenceMax = DiferenceBounds[1]

    #Should not be set below 1
    DiferenceLimitMultiplier = 8
    #Conservative 4
    #Reasonable 2
    #High reduction 1

    ValuesRemoved = 0

    Colors = []

    IndexesToRemove = []

    MaxDifferenceAllowed = MeanDifference * 2 * DiferenceLimitMultiplier

    Count = 0
    for Difference in Differences:
        #If diference between the input YValue and the models YValue is greater than the mean diference of all data points * 2 * 'DiferenceLimitMultiplier'
        #Then remove the data point

        #Justification
        #On average, all points should be between 0 * mean and  2 * mean
        #Thus the mean is in between 0 and 2 times any given value
        #So multiplying the mean diference value by 2 should cover all 'average' data points
        #Multiplying this new value by the 'DiferenceLimitMultiplier' increases the valid bounds to increase the error tolerance
        #This allows for fine tuning, as all values will not likely fit between 0, and 2 times the mean

        if ((Difference) > MaxDifferenceAllowed):
            ValuesRemoved+=1
            IndexesToRemove.append(Count)

            if (TestMode):
                Colors.append((0.5, 0, 0, 0.75))
        else:
            if (TestMode):
                Colors.append((0, 0.5, 0, 0.75))
        Count += 1

    IndexesToRemove = np.array(IndexesToRemove)

    if (len(IndexesToRemove) > 0):
        NewDataX = np.delete(NewDataX, IndexesToRemove)
        NewDataY = np.delete(NewDataY, IndexesToRemove)

    NewNumberOfDataPoints = len(NewDataY)

    #Debug visuals
    if (TestMode):

        PlotType = 1

        if (PlotType == 0):
            #Will plot all input data points, green if not removed, red if removed

            #Difference Value
            #matplot.scatter(DataX, Differences, color = Colors)

            #X,Y Value
            matplot.scatter(DataX, DataY, color=Colors)

        else:
            if (PlotType == 1):
                HeatMapColors = []
                Count = 0
                for Difference in Differences:
                    if (not ((Difference) > MaxDifferenceAllowed)):
                        HeatMapColors.append((1, 0, 0, Clamp((1.0 / MaxDifferenceAllowed * (Difference)), 0.0, 1.0)))
                    else:
                        #HeatMapColors.append((1.0,0.647,0,0.9))
                        HeatMapColors.append((0, 0, 0, 0.9))

                matplot.scatter(DataX, DataY, color=HeatMapColors, s=8)



        XBounds = GetArrayBounds(DataX)
        SamplePoints = np.linspace(XBounds[0], XBounds[1], 10000)
        m = batman.TransitModel(BatmanParams, SamplePoints,nthreads = BatmansThreads)
        LightCurve = m.light_curve(BatmanParams) # m.light_curve(BatmanParams) + (c0 + c1*Xals + c2*Xvals**2)
        LightCurve = ApplyPolyMultiplier(SamplePoints, LightCurve, Parameters)

        matplot.plot(SamplePoints, LightCurve, "-", color="blue")

        if (HighlightBoundsMode):
            matplot.plot(SamplePoints, LightCurve + MaxDifferenceAllowed, "-", color="green")
            matplot.plot(SamplePoints, LightCurve - MaxDifferenceAllowed, "-", color="green")

        if (not OverlayMode):
            #Close graph to continue program
            matplot.show()

        #Results Interpretation:

        #Green lines are bounds in which points are allowd, points that do not fall between these lines will be removed
        #Blue line is the fitted line, it's y values are what the points are being compared to
        #Black points have been removed
        #Red points have not been removed, the transparency of their color indicates how close they are to being removed. Dark red ones are close to the limit, almost clear ones are close to their expected values.
        #Note ^ points overlayed on top of eachother will stack their transparencies, resulting in dark coolors even if they are not close to the border. Zoom in on the graph if you wish to acurately see the colors of points that are overlaping one another.

    if (NewNumberOfDataPoints + len(IndexesToRemove) != len(DataY)):
        #New number of datapoints value given is not equal to the original number minus those removed
        #This means there is an issue with the 'RemoveOutliersFromDataSet' function
        print("ERROR : 'RemoveOutliersFromDataSet' function returned an improper value")

        #This check is here because small diferences in the actual number of data values [len(DataX)] and the recorded number of data points [NumberOfDataPoints] values can easilly go unoticed and lead to inacurate ChiSqr values down the line

        #Stop, definitely better way to do this
        time.sleep(99999)

    RemovalPercentage = (100/len(DataX)*ValuesRemoved)
    print("Values Removed By Outlier Rejection : " + str(ValuesRemoved) + " : Percentage Of Total Values : " + str(int(RemovalPercentage*100)/100) + "%")
    if(RemovalPercentage > 7.5):
        print("Warning :",str(int(RemovalPercentage*100)/100),"% of input values removed by outlier rejection function. This is concerningly high. Is the fit being used to remove outliers too inaccurate?")
    return (NewDataX, NewDataY, NewNumberOfDataPoints, IndexesToRemove)


def ExtractTransitParametersFromFittedFunction(Function):

    #Have to manually assign params instead of using 'OptimizedFunction.params' because 'limb_dark' is not assigned by the function

    return (ConvertFitParametersToTransitParameters(Function))


def GetArrayIsNotNone(InputArray):
    #Should (type(InputArray) == "NoneType") be used instead?
    return (InputArray is not None)


def EndTimeRecording():

    global StartTimes
    global EndTimes
    global ProgramStartTime

    if (len(StartTimes) > 0):
        print("\n----- Block Percents Of Total Time -----")
        for i in range(len(StartTimes)):
            #Block times are only valid if ID's were referenced only once
            print("Block Percent Of Total Time : " + str(100.0 / (time.time() - ProgramStartTime) * (EndTimes[i] - StartTimes[i])))

    print("\n----- Total Time -----")
    print("Total Time : " + str(time.time() - ProgramStartTime))


def ConvertFitParametersToTransitParameters(InputParam):

    #All parameters values being sent from one fucntion to another, shall now be sent as dictionaries
    #No exceptions
    #This is to simplify things, as right now there are 3 types of parameters being juggled arround and I am no longer able to keep track of which is which

    FittingTransitFunctionParams = batman.TransitParams()

    FittingTransitFunctionParams.t0 = InputParam["t0"]  #time of inferior conjunction
    FittingTransitFunctionParams.per = InputParam["per"]  #orbital period
    FittingTransitFunctionParams.rp = InputParam["rp"]  #planet radius (in units of stellar radii)
    FittingTransitFunctionParams.a = InputParam["a"]  #semi-major axis (in units of stellar radii)
    FittingTransitFunctionParams.inc = InputParam["inc"]  #orbital inclination (in degrees)
    FittingTransitFunctionParams.ecc = InputParam["ecc"]  #eccentricity
    FittingTransitFunctionParams.w = InputParam["w"]  #longitude of periastron (in degrees)
    FittingTransitFunctionParams.limb_dark = "quadratic"  #limb darkening model
    FittingTransitFunctionParams.u = [InputParam["u1"], InputParam["u2"]]  #limb darkening coefficients [u1, u2]

    return (FittingTransitFunctionParams)

def ApplyPolyMultiplier(XVal, YVal,  Params):

    #POLYVAL RETURNS IN REVERSE ORDER
    #x^9+x^8+x^7...

    global PolynomialOrder

    if(PolynomialOrder > 0):
        
        PolyVals = []
        for PolyVal in range(0,PolynomialOrder):
            #if(PolyVal == PolynomialOrder-1):
            #    PolyVals.append(1)
            #else:
            #    PolyVals.append(0)
            PolyVals.append(Params["PolyVal" + str(PolyVal)])

        #print(PolyVals)
        return(YVal + np.polyval(PolyVals,XVal))
        
    else:
        #Can not apply polyvals
        #Return unmodifed array
        return(YVal)

def ContinouseDrawGraph(XVal, YVal, Parameters):
        BatmanParams = ConvertFitParametersToTransitParameters(Parameters)

        XBounds = GetArrayBounds(XVal)
        SamplePoints = np.linspace(XBounds[0], XBounds[1], 10000)
        m = batman.TransitModel(BatmanParams, SamplePoints,nthreads = BatmansThreads)
        LightCurve = m.light_curve(BatmanParams)
        LightCurve = ApplyPolyMultiplier(SamplePoints, LightCurve, Parameters)

        matplot.cla()
        matplot.scatter(XVal,YVal)
        matplot.plot(SamplePoints, LightCurve, "-", color="orange")

        matplot.pause(0.01)


TestAvergageTimeMode = False

#Debug only, turns on interative matplot drawing mode
DrawProgressiveFitGraph = False

if(DrawProgressiveFitGraph):
    matplot.ion()


if (not TestAvergageTimeMode):
    RunOptimizationOnDataInputFile(Priors)
else:
    MultiCountStartTime = time.time()
    Iterations = 10

    for i in range(0, Iterations):
        RunOptimizationOnDataInputFile(Priors)
        print("\n============")
        print("Percent Completed :", str(100 / Iterations * (i + 1)))
        print("\n============")

    print(("\n") * 30, "============\nFINISHED\n============\n\nAverage Time :", str(int((time.time() - MultiCountStartTime) / Iterations * 100) / 100))