﻿from statistics import stdev
import numpy as np

#Curve fitting package
import lmfit as lmfit
from lmfit import fit_report

#Graph rendering package
import matplotlib.pyplot as matplot

#Alternative to text boxes that allows right alignment
from matplotlib.offsetbox import AnchoredText

#Transit curve calculation package
import batman

#WARNING : Only works if OpenMP is enabled
#Do not set this value above 1 if OpenMP is not available
BatmansThreads = 1

#NOTE: This program assumes only 2 limb-darkening values, more are possible with modification
ParamNames = 't0', 'per', 'rp', 'a', 'inc', 'ecc', 'w', 'u1', 'u2', 'PolynomialOrder'

#Debug Logging
DebugLog = False

class FitData():
    def __init__(self):

        self.Time = None
        self.Flux = None

        self.Uncertainty = None

        self.t0 =     [0.96, 999],  #t0
        self.per=[3.14159265, 999],  #per
        self.rp=[0.118, 999],  #rp
        self.a=[8.0, 999],  #a
        self.inc=[83.7, 999],  #inc
        self.ecc=[0.0, 999],  #ecc
        self.w=[1.0, 999],  #w
        self.u1=[-0.9, 999],  #u1
        self.u2=[-0.9, 999],  #u2

        self.PolynomialOrder = None


#Visual modifier only
DataPointRenderSize = 1
#Suggested Values:
#1, apears as point cloud, easy to see function line
#2, larger than error bar, easy to see points when error value is high

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

def CalculateChiSqr(DataX, DataY, DataError, Priors, Params, LBMMode):

    flux = ConvertParamatersToFluxValues(DataX, Params)

    '''
    if(Output):
        print("---------------------")
        for Val in DataError:
            print(Val)
        print("\n=======\n")

        #for Val in flux:
        #    print(Val)
        #print("\n=======\n")

        #if(DataError is not None):
        #    for Val in DataError:
        #        print(Val)
        print("---------------------")
    '''

    CheckedOptimizedChiSqr = 0

    DataIncludedErrorBars = (DataError is not None)

    #sumation of ((data_i - model_i) / uncertainty_i)
    CheckedOptimizedChiSqr = (DataY - flux)

    #Divide by uncertainty, because large uncetainty means the chi value for that data point is less revelant, and should not affect the overall chi value as much
    #Low uncertainty values are highly relevant to the final chi sum, so if the deviation is 1, but certainy is 0.1, then that 1/0.1 = 10, afffecting the total chi value more
    #Because even a small deivation at that level of uncertainty is relevant
    if (DataIncludedErrorBars):
        #print(stdev(abs(CheckedOptimizedChiSqr)-abs(DataError)))
        CheckedOptimizedChiSqr /= DataError

    #Apply priors by appending them to the returned chi array
    ParameterValues = Params.valuesdict()
    if(Priors is not None):
        global ParamNames

        for ParamName in ParamNames:
            if(ParamName != "PolynomialOrder"):
                #Don't apply prior if no uncertiainty is given, priors without uncertainty will still be used as initial fitting values but cannot be applied to the chisqr value
                if Priors[ParamName][1] is not None:
                    NewValue = (((Priors[ParamName][0] - ParameterValues[ParamName]) / Priors[ParamName][1]))
                    CheckedOptimizedChiSqr = np.append(CheckedOptimizedChiSqr, [NewValue])
    

    #IMPORTANT
    #Normally the chi array needs to be squared, this returns a positive value and is the convention for 'chisqr' tests
    #But the LBM fit method inputs an array of chi values, and then squares them itself as part of it's evaluation
    #So squareing it here will result in LBM squaring it again later, which creates an unrealistically low or high value
    #The solution is to not square the returned value when in LMBMode.
    if(not LBMMode):
        CheckedOptimizedChiSqr = CheckedOptimizedChiSqr**2
    else:
        CheckedOptimizedChiSqr = abs(CheckedOptimizedChiSqr)

    #if(LBMMode):
    #    print(str((CheckedOptimizedChiSqr**2).sum()) + "   /   " + str(len(DataY)) + "   =   " + str((((CheckedOptimizedChiSqr**2).sum())/len(DataY))))
    #else:
    #    print(str((CheckedOptimizedChiSqr).sum()) + "   /   " + str(len(DataY)) + "   =   " + str((((CheckedOptimizedChiSqr).sum())/len(DataY))))

    return (CheckedOptimizedChiSqr)

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
LBMEvaluations = 0

def ParameterEvaluationFunction(Params, DataX, DataY, DataError, Priors, IsNelder):

    if(IsNelder):
        global NelderEvaluations
        NelderEvaluations+=1
    else:
        global LBMEvaluations
        LBMEvaluations+=1

    ReturnChiArray = CalculateChiSqr(DataX, DataY, DataError, Priors, Params, (not IsNelder))

    #Draws graph after each fit, only for debugging, very slow
    #Comment out for normal use
    #if(not IsNelder):
    #    ContinouseDrawGraph(DataX, DataY, Params)

    return (ReturnChiArray)

def OptimizeFunctionParameters(DataX, DataY, DataError, Priors, UseLBM, StartingParameters):

    Bounds = GetArrayBounds(DataX)
    MinX = Bounds[0]
    MaxX = Bounds[1]

    PolynomialOrder = -1

    if(StartingParameters is not None):
        PolynomialOrder = StartingParameters["PolynomialOrder"].value
    else:
        if(Priors is not None):
            PolynomialOrder = Priors["PolynomialOrder"]

    InputParams = lmfit.Parameters()

    if ((StartingParameters is not None) or (Priors is not None)):
        

        UseParameters = StartingParameters is not None

        

        #Starting parameters are passed from one fitting attempt to another, they let the fitter pick up where it left off and therefore are used as starting values instead of the priors when they are supplied
        #If these are not given, and priors are provided, the priors will be used,because they are better than an arbitrary guess on where to start the fitting process

        if(UseParameters):
            AccessDict = StartingParameters.valuesdict()
        else:
            PolyExcludedPriors = Priors
            PolyExcludedPriors.pop("PolynomialOrder")

            AccessDict = dict()
            for Key,Val in PolyExcludedPriors.items():
                AccessDict[Key] = Val[0]
        
        InputParams.add("t0",value=AccessDict["t0"],min=MinX-(MaxX-MinX),max=MaxX+(MaxX-MinX))  #Max?
        InputParams.add("per",value=AccessDict["per"],min=0.000001,max=MaxX-MinX)
        InputParams.add("rp", value=AccessDict["rp"], min=0, max=1.0)
        InputParams.add("a", value=AccessDict["a"], min=1.0,max=90)  #In ratio of stellar radii
        InputParams.add("inc",value=AccessDict["inc"],min=60,max=90) # Degrees from "top" of star to orbit plane
        InputParams.add("ecc",value=AccessDict["ecc"],min=0.0,max=1.0)
        InputParams.add("w",value=AccessDict["w"],min=0.0,max=360.0) # from 0-360 or -180 to 180 : think 0-360 based on w*pi/180 
        InputParams.add("u1",value=AccessDict["u1"],min=0.0,max=1.0)
        InputParams.add("u2",value=AccessDict["u2"],min=0.0,max=1.0)
        InputParams.add("PolynomialOrder", value=PolynomialOrder, vary = False)
        
        if(PolynomialOrder != -1):
            for PolyIndex in range(0,PolynomialOrder+1):
                PolyName = ("PolyVal" + str(PolyIndex))
                StartingVal = 0
                if(PolyIndex == 0):
                        StartingVal = 1

                if(UseParameters):
                    #Parametrs already generated mode
                    #if(PolyName in StartingParameters):
                    StartingVal = StartingParameters[PolyName].value
                    InputParams.add(PolyName, value=StartingVal, min=-1000, max=1000, vary = True)
                else:
                    #No initial values
                    InputParams.add(PolyName, value=StartingVal, min=-1000, max=1000, vary = True)

    else:
        #Backup - Will result in bad fit if not given a starting point, be that initial params or priors
        InputParams.add("t0", value=MaxX/2, min=MinX, max=MaxX)  #Max?
        InputParams.add("per", value=MaxX/2, min=0.0, max=MaxX)
        InputParams.add("rp", value=0.5, min=0, max=1.0)
        InputParams.add("a", value=40, min=1.0,max=90)  #What should Max Bound be?
        InputParams.add("inc", value=70, min=60, max=90)
        InputParams.add("ecc", value=0.5, min=0.0, max=1.0)
        InputParams.add("w", value=40, min=0.0, max=360.0)
        InputParams.add("u1", value=0.5, min=-1.0, max=1.0)
        InputParams.add("u2", value=0.5, min=-1.0, max=1.0)
        InputParams.add("PolynomialOrder", value=PolynomialOrder, vary = False)

        if(PolynomialOrder != -1):
            for PolyIndex in range(0,PolynomialOrder+1):
                PolyName = ("PolyVal" + str(PolyIndex))
                StartingVal = 0
                if(PolyIndex == 0):
                    StartingVal = 1
                InputParams.add(PolyName, value=StartingVal, min=-1000, max=1000, vary = True)

    OptimizedFunctionToReturn = None

    if (not UseLBM):
        OptimizedFunctionToReturn = lmfit.minimize(
            ParameterEvaluationFunction,
            InputParams,
            args=(DataX, DataY, DataError, Priors, True),
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
            ParameterEvaluationFunction,
            InputParams,
            args=(DataX, DataY, DataError, Priors, False),
            method="leastsq",
            calc_covar=True,
            max_nfev=None,
            nan_policy="raise")



        #Weight Implementation According To lmfit Documentation:
        '''
        weights
        numpy.ndarray (or None) of weighting values to be used in fit.
        If not None, it will be used as a multiplicative factor of the residual array, so that weights*(data - fit) is minimized in the least-squares sense.
        '''

    return (OptimizedFunctionToReturn)

def RemoveOutliersFromDataSet(DataX, DataY, Parameters):

    TestMode = False
    
    #-Only valid if 'TestMode' is active

    #will not halt further execution of the program, instead overlaying the scatter points or heat map underneath the final graph
    OverlayMode = False

    #Show limits values are allowed between
    HighlightBoundsMode = True

    #-----------------------------------

    NewDataX = DataX
    NewDataY = DataY
    NewNumberOfDataPoints = -1

    LightCurve = ConvertParamatersToFluxValues(DataX, Parameters)

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

        LightCurve = ConvertParamatersToFluxValues(SamplePoints, Parameters)

        '''
        TransitModel = batman.TransitModel(BatmanParams, SamplePoints,nthreads = BatmansThreads)
        #Confirm poly aplication
        LightCurve = TransitModel.light_curve(BatmanParams) # TransitModel.light_curve(BatmanParams) + (c0 + c1*Xals + c2*Xvals**2)
        LightCurve = ApplyPolyMultiplier(SamplePoints, LightCurve, Parameters)
        '''

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
        return

    RemovalPercentage = (100/len(DataX)*ValuesRemoved)
    global DebugLog
    if(DebugLog):
        print("Values Removed By Outlier Rejection : " + str(ValuesRemoved) + " : Percentage Of Total Values : " + str(int(RemovalPercentage*100)/100) + "%")
    if(RemovalPercentage > 7.5):
        print("Warning :",str(int(RemovalPercentage*100)/100),"% of input values removed by outlier rejection function. This is concerningly high. Is the fit being used to remove outliers too inaccurate?")
    return (NewDataX, NewDataY, NewNumberOfDataPoints, IndexesToRemove)

def ExtractTransitParametersFromFittedFunction(Function):

    #Have to manually assign params instead of using 'OptimizedFunction.params' because 'limb_dark' is not assigned by the function

    return (ConvertFitParametersToTransitParameters(Function))

def ConvertFitParametersToTransitParameters(InputParam):

    #All parameters values being sent from one fucntion to another, shall now be sent as param dictionaries
    #No exceptions
    #This is to simplify things, as right now there are 3 types of parameter formats being juggled arround and I am no longer able to keep track of which is which

    FittingTransitFunctionParams = batman.TransitParams()
    ParamValues = InputParam.valuesdict()

    FittingTransitFunctionParams.t0 = ParamValues["t0"]  #time of inferior conjunction
    FittingTransitFunctionParams.per = ParamValues["per"]  #orbital period
    FittingTransitFunctionParams.rp = ParamValues["rp"]  #planet radius (in units of stellar radii)
    FittingTransitFunctionParams.a = ParamValues["a"]  #semi-major axis (in units of stellar radii)
    FittingTransitFunctionParams.inc = ParamValues["inc"]  #orbital inclination (in degrees)
    FittingTransitFunctionParams.ecc = ParamValues["ecc"]  #eccentricity
    FittingTransitFunctionParams.w = ParamValues["w"]  #longitude of periastron (in degrees)
    FittingTransitFunctionParams.limb_dark = "quadratic"  #limb darkening model
    FittingTransitFunctionParams.u = [ParamValues["u1"], ParamValues["u2"]]  #limb darkening coefficients [u1, u2]

    return (FittingTransitFunctionParams)

def ConvertParamatersToFluxValues(XVal, Params):

    global  BatmansThreads

    BatmanParams = ConvertFitParametersToTransitParameters(Params)

    TransitModel = batman.TransitModel(BatmanParams, XVal,nthreads = BatmansThreads)
    LightCurve = TransitModel.light_curve(BatmanParams)
    LightCurve = ApplyPolyMultiplier(XVal, LightCurve, Params)

    return(LightCurve)

def ApplyPolyMultiplier(XVal, YVal,  Params):

    #POLYVAL RETURNS IN REVERSE ORDER
    #x^9+x^8+x^7...

    PolynomialOrder = Params["PolynomialOrder"].value

    if(PolynomialOrder != -1):
        
        PolyVals = []

        for PolyVal in range(0,PolynomialOrder+1):
            PolyVals.append(Params["PolyVal" + str(PolyVal)].value)

        return(YVal * np.polyval(PolyVals,XVal))
    else:
        #Can not apply polyvals
        #Return unmodifed array
        return(YVal)

def ContinouseDrawGraph(XVal, YVal, Parameters):
        #Draws matplot graph of fit wihout halting program
        #Requires matplot interactive mode to be inabled with matplot.ion()
        #Used for when I don't understadn why a fit it failing and want to see the process it is taking

        #BatmanParams = ConvertFitParametersToTransitParameters(Parameters)

        XBounds = GetArrayBounds(XVal)
        SamplePoints = np.linspace(XBounds[0], XBounds[1], 10000)
        Flux = ConvertParamatersToFluxValues(SamplePoints, Parameters)

        #TransitModel = batman.TransitModel(BatmanParams, SamplePoints,nthreads = BatmansThreads)
        #LightCurve = TransitModel.light_curve(BatmanParams)
        #LightCurve = ApplyPolyMultiplier(SamplePoints, LightCurve, Parameters)



        matplot.scatter(XVal,YVal, DataPointRenderSize)
        matplot.plot(SamplePoints, Flux, "-", color="orange")
        matplot.show()


        matplot.pause(0.01)
        matplot.cla()

def CalculateDataUncertainty(DataY, FunctionY):
    #DataError =(DataY*0) +(  (abs((DataY-FirstOptimizedFunction))).sum()/len(DataY))
    #Average = abs((DataY - FunctionY).sum()/len(DataY))
    #DataError =ReplaceZerosInArrayWithLowestValue(DataY*0 + Average)
    #DataError = abs((DataY - FunctionY))
    #DataError = abs(DataY-FunctionY)
    DataError = DataY*0 + stdev((DataY-FunctionY))

    #DataError = DataY*0 + 1
    #print("Calculated Error : " + str(DataError))
    return(DataError)

    #CheckChiSqr function expects data error as an array, this allows compatibilty with lmfit.fit instead of lmfit.minimize
    #So an array is created even when a single scaler value is being used

def FitTransitFromData(InputFitData):

    Priors = dict()

    Priors["t0"] = InputFitData.t0
    Priors["per"] = InputFitData.per
    Priors["rp"] = InputFitData.rp
    Priors["a"] = InputFitData.a
    Priors["inc"] = InputFitData.inc
    Priors["ecc"] = InputFitData.ecc
    Priors["w"] = InputFitData.w
    Priors["u1"] = InputFitData.u1
    Priors["u2"] = InputFitData.u2
    Priors["PolynomialOrder"] = InputFitData.PolynomialOrder

    if(Priors["PolynomialOrder"] is None):
        Priors["PolynomialOrder"] = -1

    PolynomialOrder = Priors["PolynomialOrder"]

    print("Running. Please wait...\n")


    #Get data bounds

    DataX = np.array(InputFitData.Time)
    DataY = np.array(InputFitData.Flux)

    DataError = None
    DataIncludedErrorBars = False

    if(InputFitData.Uncertainty is not None):
        DataError = InputFitData.Uncertainty
        DataIncludedErrorBars = True
    else:
        DataIncludedErrorBars = False

    #Used in some parameter bounds
    Bounds = GetArrayBounds(DataX)
    MinX = Bounds[0]
    MaxX = Bounds[1]

    #Not currently used
    Bounds = GetArrayBounds(DataY)
    MinY = Bounds[0]
    MaxY = Bounds[1]

    NumberOfDataPoints = len(DataX)

    global NelderEvaluations
    NelderEvaluations = 0
    global LBMEvaluations
    LBMEvaluations = 0

    #Running this first iteration with 'nelder' generates significantly better initial results than 'leastsqr' (LBM) when the starting values are innacurate. If good starting values are given this can be swapped to LBM for significantly faster evaluation.
    #Running using 'leastsqr' results in a badly fit graph that is then used to remove outlier data points. This bad graph leads to good data points being thrown out, and the final graph is bad because of it.


    FirstOptimizedFunction = OptimizeFunctionParameters(DataX, DataY, DataError, Priors, False, None)
    
    #print("Nelder Evaluations :",str(NelderEvaluations),": Data Points :",len(DataX),": Dif :",str(NelderEvaluations/len(DataX)))

    #Extract parameters used
    FirstOptimizedParamsDictionary = FirstOptimizedFunction.params


    #Generate flux based on extracted parameters
    FirstOptimizedFunctionFlux = ConvertParamatersToFluxValues(DataX, FirstOptimizedParamsDictionary)

    #ContinouseDrawGraph(DataX, DataY, FirstOptimizedParamsDictionary)

    #Disable this to see if (too many)/(good) data points are being removed after the first fit.
    #If this is happening the priors are likely too restrictive or too far from the actual values.
    RemoveOutliers = True

    #print("R 1   ",DataError.sum()/len(DataError))

    '''
    if(RemoveOutliers):
        #Remove Outlier values
        NewDataValues = RemoveOutliersFromDataSet(DataX, DataY, FirstOptimizedParamsDictionary)

        DataX = NewDataValues[0]
        DataY = NewDataValues[1]
        NumberOfDataPoints = NewDataValues[2]
        IndexesRemoved = NewDataValues[3]

        #Recalculate error
        if (not DataIncludedErrorBars):
            #Data did not include error, so need to calculate it
            if(len(IndexesRemoved) > 0):
                #Outliers have been removed, so need to removed rejected points from FirstOptimizedFunctionFlux before calculating uncertainty
                FirstOptimizedFunctionFlux = np.delete(FirstOptimizedFunctionFlux,IndexesRemoved)

            #Calculate error from diference between first attempt created function, and given values
            DataError = CalculateDataUncertainty(DataY, FirstOptimizedFunctionFlux)
        else:
            #Data had error included
            if(len(IndexesRemoved) > 0):
                #Points have been removed, so DataError array needs to have these indexes removed
                DataError = np.delete(DataError, IndexesRemoved)
    else:
        #Calculate error from diference between first attempt created function, and given values
        if (not DataIncludedErrorBars):
            DataError = CalculateDataUncertainty(DataY, FirstOptimizedFunctionFlux)
    '''

        
    if(RemoveOutliers):
        #Remove Outlier values
        NewDataValues = RemoveOutliersFromDataSet(DataX, DataY, FirstOptimizedParamsDictionary)

        DataX = NewDataValues[0]
        DataY = NewDataValues[1]
        NumberOfDataPoints = NewDataValues[2]
        IndexesRemoved = NewDataValues[3]

        if(len(IndexesRemoved) > 0):
            #Outliers have been removed, so need to removed rejected points from FirstOptimizedFunctionFlux before calculating uncertainty
            FirstOptimizedFunctionFlux = np.delete(FirstOptimizedFunctionFlux,IndexesRemoved)
    

    #Run second time, this time having removed outliers and calculated error values if they were not provided

    #SecondOptimizedFunction = OptimizedFunction
    SecondOptimizedFunction = OptimizeFunctionParameters(DataX, DataY, DataError, Priors, True, FirstOptimizedParamsDictionary)

    DictionaryParams = SecondOptimizedFunction.params
    BatmanParams = ExtractTransitParametersFromFittedFunction(DictionaryParams)
    SecondOptimizedFunctionFlux = ConvertParamatersToFluxValues(DataX, FirstOptimizedParamsDictionary)

    #Uncertainty Calculation
    DataError = CalculateDataUncertainty(DataY, SecondOptimizedFunctionFlux)

    #From this point forward, "DataIncludedErrorBars" will no longer be used to decide if Error values need to be calculated
    #It just means wether to use the current Error np.array for Chi calculations and wether to render points with error bars
    DataIncludedErrorBars = True

    #Debug Fit Report
    global DebugLog
    if(DebugLog):
        print(fit_report(SecondOptimizedFunction))
        print("\n")
        fig = matplot.figure()
        #Display points with error bars
        if (DataIncludedErrorBars):
            #Disabled for now
            #Concerned the 'yerr' value being set to the DataError value is not an acurate representation of the points error

            #matplot.errorbar(DataX, DataY, yerr = DataError, fmt ="o", markersize = DataPointRenderSize)
            matplot.errorbar(DataX, DataY, fmt="o", markersize=DataPointRenderSize)
        else:
            matplot.scatter(DataX, DataY, fmt="o", markersize=DataPointRenderSize)

        print("-OPTIMIZED PARAMETERS-")

        for Key in Priors.keys():
            print(Key + ' : ' + str(DictionaryParams[Key].value))
        if(PolynomialOrder != -1):
            for PolyVal in range(0,PolynomialOrder+1):
                print(("PolyVal" + str(PolyVal)) + ' : ' + str(DictionaryParams["PolyVal" + str(PolyVal)].value))
        print("\n")

        print("-OPTIMIZED PARAMETER UNCERTAINTY VALUES-")
        for Key in Priors.keys():
            print(Key + ' : ' + str(DictionaryParams[Key].stderr))
        if(PolynomialOrder != -1):
            for PolyVal in range(0,PolynomialOrder+1):
                print(("PolyVal" + str(PolyVal)) + ' : ' + str(DictionaryParams["PolyVal" + str(PolyVal)].stderr))
        print("\n")



    DataError = CalculateDataUncertainty(DataY, ConvertParamatersToFluxValues(DataX, DictionaryParams))

    '''
    print("[][][][][][][][][][][][][]")
    for Val in  DataY:
        print(Val)
    print("[][][][][][][][][][][][][]")
    '''
    '''
    print("[][][][][][][][][][][][][]")
    for Val in  ConvertParamatersToFluxValues(DataX, DictionaryParams):
        print(Val)
    print("[][][][][][][][][][][][][]")
    '''
    '''
    print("[][][][][][][][][][][][][]")
    for Val in  DataError:
        print(Val)
    print("[][][][][][][][][][][][][]")
    '''


    CheckedOptimizedChiSqr = CalculateChiSqr(DataX, DataY, DataError, Priors, DictionaryParams, False).sum()
    if(DebugLog):
        print(CheckedOptimizedChiSqr)

    #Rendering only, uses more sample points than input x-values
    SamplePoints = np.linspace(MinX, MaxX, 10000)
    flux = ConvertParamatersToFluxValues(SamplePoints, DictionaryParams)
    matplot.plot(SamplePoints, flux, "-", label="Optimized Function")

    if(DebugLog):
        print("\n--- Checked Chi Sqr ---")
        print("ChiSqr : " + str(CheckedOptimizedChiSqr))
        print("Number Of Data Points : " + str(NumberOfDataPoints))
        #The value below should be close to '1'
        print("ChiSqr / # Data Points : " + str(CheckedOptimizedChiSqr / NumberOfDataPoints))

    #Debug Logging END

    #Fixed "χ2" rendering issue
    BestChi = (round(CheckedOptimizedChiSqr, 2))
    BestChiTxt = "Optimized χ2 : " + str(BestChi)

    ReducedChi = (round(CheckedOptimizedChiSqr / NumberOfDataPoints, 5))
    ReducedChiTxt = "Reduced χ2 : " + str(ReducedChi)

    NoPriorReducedChi = (round((CalculateChiSqr(DataX, DataY, DataError, None, DictionaryParams, False).sum() / NumberOfDataPoints), 5))
    NoPriorReducedChiTxt = "No Priors Reduced χ2 : " + str(NoPriorReducedChi)

    if(DebugLog):
        #Text box setup
        ChiAnchoredTextBox = AnchoredText((BestChiTxt + "\n" + ReducedChiTxt + "\n" + NoPriorReducedChiTxt), loc=4, pad=0.5)
        matplot.setp(ChiAnchoredTextBox.patch, facecolor="Orange", alpha=0.5)
        matplot.gca().add_artist(ChiAnchoredTextBox)
        matplot.legend(loc=2, borderaxespad=0)



    print("Completed")

    if(DebugLog):
        print("\nNelder Evaluations :",NelderEvaluations)
        NelderEvaluations = 0

        print("LBM Evaluations :",LBMEvaluations)
        LBMEvaluations = 0

    '''
    #Display plotted data
    #Code after this function is called will not be run untill the graph is closed (this behavior can be changed)
    global TestAvergageTimeMode
    if(not TestAvergageTimeMode):
        matplot.show()
    '''
    if(DebugLog):
        matplot.show()
    ReturnDictionary = dict();

    for Key,Value in DictionaryParams.valuesdict().items():
        if(Key != "PolynomialOrder"):
            ReturnDictionary[Key] = [DictionaryParams[Key].value, DictionaryParams[Key].stderr]
        else:
            ReturnDictionary[Key] = Value

    #fig.savefig("output_plot.jpg", dpi=800)
    ReturnDictionary["ChiSqr"] = BestChi;
    ReturnDictionary["ReducedChiSqr"] = ReducedChi;

    return(ReturnDictionary)



class TransitCurveFitter:
    def FitTransit(InputFitData):
        #matplot.ion()
        ReturnData = FitTransitFromData(InputFitData);
        return(ReturnData)
