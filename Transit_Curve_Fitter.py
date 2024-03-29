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

#Debug Logging - Will print info to console and draw a matplot figure of the output
DebugLog = False
#Visual modifier only for debug mode
DataPointRenderSize = 1

#Class used to store data to be fitted
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

    CheckedOptimizedChiSqr = 0

    DataIncludedErrorBars = (DataError is not None)

    #(data - model)
    CheckedOptimizedChiSqr = (DataY - flux)

    #Divide by uncertainty, because large uncetainty means the chi value for that data point is less revelant, and should not affect the overall chi value as much
    #Low uncertainty values are highly relevant to the final chi sum, so if the deviation is 1, but certainy is 0.1, then that 1/0.1 = 10, afffecting the total chi value more
    #Because even a small deviation at that level of uncertainty is relevant
    if (DataIncludedErrorBars):
        #(data - model) / uncertainty
        CheckedOptimizedChiSqr /= DataError

    #Apply priors by appending them to the returned chi array
    #This will penalize the lmfit fitting algorythims for deviating from the given priors, weighted by their uncertainty
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
    #But the LBM fit method inputs an array of chi values, and then squares them itself as part of its internal evaluation
    #So squareing it here will result in LBM squaring it again later, which creates an unrealistically low or high value
    #The solution is to not square the returned value when in LMBMode.
    #When requestng a chi value for debug or testing, request it with (LBMMode = False) so the result is sqared before being returned
    if(not LBMMode):
        CheckedOptimizedChiSqr = CheckedOptimizedChiSqr**2
    else:
        CheckedOptimizedChiSqr = abs(CheckedOptimizedChiSqr)

    #Always returns as an array, never a scalar
    return (CheckedOptimizedChiSqr)

def Clamp(Value, Min, Max):
    if (Value < Min):
        return (Min)
    if (Value > Max):
        return (Max)
    return (Value)

def GetArrayBounds(DataArray):
    DataArray = np.array(DataArray, copy=False)
    return ([DataArray.min(), DataArray.max()])

#Used for debugging
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

    #Draws matplot graph after each fit, only for debugging, very slow
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
        #If these are not given, and priors are provided, the priors will be used, because they are better than an arbitrary guess on where to start the fitting process

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

        #   lmfit.minimize
        #   (
        #Function - (evaluation function)
        #Params - (input parameters for fitting)
        #Arguments - (will be passed to function, not modified by or visible to fitting method)
        #Method - (defaults to 'Levenberg-Marquardt' referenced by 'leastsq', assigning it specifically should not be necessary [Different from 'least_squares'])
        #Calculate Uncertainty - (whether to calculate uncertainties (if fitting method supports it), should default to true)
        #Maximum evaluations - (should default to none)
        #Weights would go here - (but are instead being processed chi function instead of solver to allow for inclusion of prior weights)
        #Raise issue if NAN found - (this shouldn't happen, but it should notify the user if this happens because it means something isn't being processed correctly)
        #   )
        OptimizedFunctionToReturn = lmfit.minimize(
            ParameterEvaluationFunction,
            InputParams,
            args=(DataX, DataY, DataError, Priors, False),
            method="leastsq",
            calc_covar=True,
            max_nfev=None,
            nan_policy="raise")



        #Weight Implementation According To lmfit Documentation (Note that lmfit weights represent 'Certainty', while this script work with 'Uncertainty'):
        '''
        weights
        numpy.ndarray (or None) of weighting values to be used in fit.
        If not None, it will be used as a multiplicative factor of the residual array, so that weights*(data - fit) is minimized in the least-squares sense.
        '''

    return (OptimizedFunctionToReturn)

def RemoveOutliersFromDataSet(DataX, DataY, Parameters):

    #-----------------------------------

    #If fitter is removing too many good values from the data enable this to see the remova bounds
    TestMode = False
    
    #TestMode parameters - Only valid if 'TestMode' is active

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
    #Highly Conservative : 8
    #Conservative : 4
    #Reasonable : 2
    #High Reduction : 1

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
                        HeatMapColors.append((0, 0, 0, 0.9))

                matplot.scatter(DataX, DataY, color=HeatMapColors, s=8)



        XBounds = GetArrayBounds(DataX)
        SamplePoints = np.linspace(XBounds[0], XBounds[1], 10000)

        LightCurve = ConvertParamatersToFluxValues(SamplePoints, Parameters)

        matplot.plot(SamplePoints, LightCurve, "-", color="blue")

        if (HighlightBoundsMode):
            matplot.plot(SamplePoints, LightCurve + MaxDifferenceAllowed, "-", color="green")
            matplot.plot(SamplePoints, LightCurve - MaxDifferenceAllowed, "-", color="green")

        if (not OverlayMode):
            #Close graph to continue program
            matplot.show()

        #Results Interpretation:

        #Green lines are bounds in which points are allowed, points that do not fall between these lines will be removed
        #Blue line is the fitted line, its y values are what the points are being compared to
        #Black points have been removed
        #Red points have not been removed, the transparency of their color indicates how close they are to being removed. Dark red ones are close to the limit, almost clear ones are close to their expected values.
        #Note ^ points overlayed on top of each other will stack their transparencies, resulting in dark colors even if they are not close to the border. Zoom in on the graph if you wish to accurately see the colors of points that are overlapping one another.

    if (NewNumberOfDataPoints + len(IndexesToRemove) != len(DataY)):
        #New number of datapoints value given is not equal to the original number minus those removed
        #This means there is an issue with the 'RemoveOutliersFromDataSet' function
        #This check is here because small differences in the actual number of data values [len(DataX)] and the recorded number of data points [NumberOfDataPoints] values can easily go unnoticed and lead to inaccurate ChiSqr values later down the line
        raise Exception("'RemoveOutliersFromDataSet' function returned an improper value. The number of input points and (output points + Removed) were not equal. There is likely a bug in this script.")

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

    #All parameters’ values being sent from one function to another, shall now be sent as param dictionaries
    #No exceptions
    #This is to simplify things, as right now there are 3 types of parameter formats being juggled around and I am no longer able to keep track of which is which

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

    #POLYVAL INPUTS IN REVERSE ORDER, EX:
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
        #Draws matplot graph of fit without halting program
        #Requires matplot interactive mode to be enabled with 'matplot.ion()' before this function is called
        #Used for when I don't understand why a fit it is failing and want to see the process it is taking

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

    #DataError = DataY*0 + stdev((DataY-FunctionY))
    #Not a good guess, used to give fit functions a starting place if uncertainty of the data is not given
    DataError = DataY*0 + stdev(FunctionY)

    return(DataError)

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

    global DebugLog

    if(DebugLog):
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
    
    #Extract parameters used
    FirstOptimizedParamsDictionary = FirstOptimizedFunction.params

    #Generate flux based on extracted parameters
    FirstOptimizedFunctionFlux = ConvertParamatersToFluxValues(DataX, FirstOptimizedParamsDictionary)

    #Disable this to see if (too many)/(good) data points are being removed after the first fit.
    #If this is happening the priors are likely too restrictive or too far from the actual values.
    RemoveOutliers = True

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

            if(DataIncludedErrorBars):
                DataError = np.delete(DataError, IndexesRemoved)
    

    #Run second time, this time having removed outliers and calculated error values if they were not provided

    #SecondOptimizedFunction = OptimizedFunction
    SecondOptimizedFunction = OptimizeFunctionParameters(DataX, DataY, DataError, Priors, True, FirstOptimizedParamsDictionary)

    DictionaryParams = SecondOptimizedFunction.params
    BatmanParams = ExtractTransitParametersFromFittedFunction(DictionaryParams)
    SecondOptimizedFunctionFlux = ConvertParamatersToFluxValues(DataX, FirstOptimizedParamsDictionary)

    #Uncertainty Calculation
    DataError = CalculateDataUncertainty(DataY, SecondOptimizedFunctionFlux)

    #Debug Fit Report
    if(DebugLog):
        print(fit_report(SecondOptimizedFunction))
        print("\n")
        fig = matplot.figure()
        matplot.plot(DataX, DataY)

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

    if(not DataIncludedErrorBars):
        DataError = CalculateDataUncertainty(DataY, ConvertParamatersToFluxValues(DataX, DictionaryParams))

    CheckedOptimizedChiSqr = CalculateChiSqr(DataX, DataY, DataError, Priors, DictionaryParams, False).sum()
    if(DebugLog):
        print(CheckedOptimizedChiSqr)

    if(DebugLog):
        #Rendering only, uses more sample points than input x-values
        SamplePoints = np.linspace(MinX, MaxX, 10000)
        flux = ConvertParamatersToFluxValues(SamplePoints, DictionaryParams)
        matplot.plot(SamplePoints, flux, "-", label="Optimized Function")

    
        print("\n--- Checked Chi Sqr ---")
        print("ChiSqr : " + str(CheckedOptimizedChiSqr))
        print("Number Of Data Points : " + str(NumberOfDataPoints))
        #The value below should be close to '1'
        print("ChiSqr / # Data Points : " + str(CheckedOptimizedChiSqr / NumberOfDataPoints))

    #If file is stripped of custom charachters the "χ2" symbols may not render correctly on the debug plot, this is only a visual defect
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


    if(DebugLog):
        print("Completed")

    if(DebugLog):
        print("\nNelder Evaluations :",NelderEvaluations)
        NelderEvaluations = 0

        print("LBM Evaluations :",LBMEvaluations)
        LBMEvaluations = 0

    if(DebugLog):
        matplot.show()
        #Uncomment to save a copy of the output graph
        #fig.savefig("output_plot.jpg", dpi=800)

    ReturnDictionary = dict()

    for Key,Value in DictionaryParams.valuesdict().items():
        if(Key != "PolynomialOrder"):
            ReturnDictionary[Key] = [DictionaryParams[Key].value, DictionaryParams[Key].stderr]
        else:
            if(Value == -1):
                ReturnDictionary[Key] = None
            else:
                ReturnDictionary[Key] = Value

    
    ReturnDictionary["ChiSqr"] = BestChi
    ReturnDictionary["ReducedChiSqr"] = ReducedChi
    ReturnDictionary["TimeValues"] = DataX
    ReturnDictionary["FluxValues"] = DataY

    return(ReturnDictionary)







#Class with function to be called to create model from data stored in a 'FitData' class
#Will output results as a dictionary of the format:

#'Results' is a dictionary of information about the optimized fit in the following format:
#--------------------------
#t0 : [Value, Uncertainty]
#per : [Value, Uncertainty]
#rp : [Value, Uncertainty]
#a : [Value, Uncertainty]
#inc : [Value, Uncertainty]
#ecc : [Value, Uncertainty]
#w : [Value, Uncertainty]
#u1 : [Value, Uncertainty]
#u2 : [Value, Uncertainty]
#PolynomialOrder : Value (The input polynomial value, will not be modified and therefore has no uncertainty value)
#PolyVal0 : [Value, Uncertainty]
#PolyVal1 : [Value, Uncertainty]
#PolyVal2 : [Value, Uncertainty]
#Polyval# : Will have as many of these parameters as needed for the input polynomial order (Order 0 = 1 variable) (Order 1 = 2 variables) (order 2 = 3 variables) (etc.)
#ChiSqr : Float (Should be close to number of input data points)
#ReducedChiSqr : Float (Should be close to 1)
#TimeValues : [All input Time values, minus any values removed by the outlier rejection system]
#FluxValues : [All input Flux values, minus any values removed by the outlier rejection system]
#--------------------------

class TransitCurveFitter:
    def FitTransit(InputFitData):
        ReturnData = FitTransitFromData(InputFitData)
        return(ReturnData)