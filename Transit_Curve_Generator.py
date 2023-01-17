from cmath import exp, pi, sqrt
from datetime import time
from itertools import count
import time
from turtle import color
import numpy as np
import lmfit as lmfit
import matplotlib.pyplot as matplot

from lmfit import Model, fit_report

#Alternative to text boxes that allows right alignment
from matplotlib.offsetbox import AnchoredText

import batman

#important distinction
from lmfit import Parameters
from lmfit import Parameter

#Only used for testing
import random

#Profiling library
import cProfile



#If supplied, will be used as initial fitting parameters
Priors = [
          [0.96,None], #t0
          [3.14159265,None], #per
          [0.118,None], #rp
          [8.0,None], #a
          [83.7,None], #inc
          [0.0,None], #ecc
          [1.0,None], #w
          [-0.5,None], #u1
          [-0.9,None] #u2
          ]

#Visual modifier only
DataPointRenderSize = 2
#Suggested Values:
#0.4, apears as point cloud, easy to see function line
#2, larger than error bar, easy to see points when error value is high

#   TEST OVERIDE MODE
TestMode = False
TestModeRenderPoints = True

if(TestMode):
    OverideFunctionParams = batman.TransitParams()
    OverideFunctionParams.t0 = 0.96
    OverideFunctionParams.per = 3.1
    OverideFunctionParams.rp = 0.118
    OverideFunctionParams.a = 16.2
    OverideFunctionParams.inc = 87.3
    OverideFunctionParams.ecc = 0.1
    OverideFunctionParams.w = 90.0
    OverideFunctionParams.u = [0.09, 0.28]
    OverideFunctionParams.limb_dark = "quadratic"
    OverideFunctionMin = 0
    OverideFunctionMax = 27.19

    fig, ax = matplot.subplots()

    SamplePoints = np.linspace(OverideFunctionMin,OverideFunctionMax,100000)
    m = batman.TransitModel(OverideFunctionParams, SamplePoints)
    flux = m.light_curve(OverideFunctionParams)
    matplot.plot(SamplePoints,flux, "-", label="Optimized Function")



ProgramStartTime = time.time()
StartTime = 0
def CheckTime(IsStart):
    if(IsStart):
        global StartTime
        StartTime = time.time()
    else:
        #Is end
        StartTime = time.time()-StartTime
        print("Block Time : "  + str(StartTime));

def ReplaceZerosInArrayWithLowestValue(DataArray):
    #Replaces zeros in an array with the lowest value present in the array
    #This is normally used on an array of error values, as values of '0' leed to calculation errors and are improbable
    #EX: (Input - Output) / Error = Divide By Zero Error

        FixedArray = DataArray

        LowestArrayValue = "!"

        for Count in range(len(DataArray)):
            if((LowestArrayValue == "!" or DataArray[Count] < LowestArrayValue) and DataArray[Count] != 0):
                LowestArrayValue = DataArray[Count]

        for Count in range(len(FixedArray)):
            if(FixedArray[Count] == 0):
                print("'0' Array Value At Index : " + str(Count) + " : Replacing With Lowest Array Value")

                FixedArray[Count] = LowestArrayValue

        return(FixedArray)

GlobalCount = 0
CountTime = True
def CalculateChiSqr(DataX, DataY, Params, DataERROR, Priors):

    if(CountTime):
        global GlobalCount
        GlobalCount+=1
        #print(GlobaleCount)

    

    m = batman.TransitModel(Params, DataX)
    flux = m.light_curve(Params)

    CheckedOptimizedChiSqr = 0

    DataIncludedErrorBars = GetArrayIsNotNone(DataERROR)

    #sumation of ((data_i - model_i) / uncertainty_i)^2
    if(DataIncludedErrorBars):
        #Previousely was : CheckedOptimizedChiSqr = (((((DataY-flux)**2))/(DataERROR))).sum()
        CheckedOptimizedChiSqr = (((DataY-flux)/DataERROR)**2).sum()
    else:
        CheckedOptimizedChiSqr = (((DataY-flux))**2).sum()

    if(not Priors == None):
        #Unpackage priors, compare to Params, modify CheckedOptimizedChiSqr based on their comparison

        if(not Priors[0][0] == None and not Priors[0][1] == None):
            CheckedOptimizedChiSqr+=ReturnChiModiferOfParameterPrior(Params.t0, Priors[0][0], Priors[0][1])
        if(not Priors[1][0] == None and not Priors[1][1] == None):
            CheckedOptimizedChiSqr+=ReturnChiModiferOfParameterPrior(Params.per, Priors[1][0], Priors[1][1])
        if(not Priors[2][0] == None and not Priors[2][1] == None):
            CheckedOptimizedChiSqr+=ReturnChiModiferOfParameterPrior(Params.rp, Priors[2][0], Priors[2][1])
        if(not Priors[3][0] == None and not Priors[3][1] == None):
            CheckedOptimizedChiSqr+=ReturnChiModiferOfParameterPrior(Params.a, Priors[3][0], Priors[3][1])
        if(not Priors[4][0] == None and not Priors[4][1] == None):
            CheckedOptimizedChiSqr+=ReturnChiModiferOfParameterPrior(Params.inc, Priors[4][0], Priors[4][1])
        if(not Priors[5][0] == None and not Priors[5][1] == None):
            CheckedOptimizedChiSqr+=ReturnChiModiferOfParameterPrior(Params.ecc, Priors[5][0], Priors[5][1])
        if(not Priors[6][0] == None and not Priors[6][1] == None):
            CheckedOptimizedChiSqr+=ReturnChiModiferOfParameterPrior(Params.w, Priors[6][0], Priors[6][1])
        #Separated in case only 1 u is given, or given out of order
        if(not Priors[7][0] == None and not Priors[7][1] == None):
            CheckedOptimizedChiSqr+=ReturnChiModiferOfParameterPrior(Params.u[0], Priors[7][0], Priors[7][1])
        if(not Priors[8][0] == None and not Priors[8][1] == None):
            CheckedOptimizedChiSqr+=ReturnChiModiferOfParameterPrior(Params.u[1], Priors[8][0], Priors[8][1])


    return(CheckedOptimizedChiSqr)

def ReturnChiModiferOfParameterPrior(Param, Prior, PriorERROR):
    return(((Param-Prior)/PriorERROR)**2)

def Clamp(Value, Min, Max):
    if(Value < Min):
        return(Min)
    if(Value > Max):
        return(Max)
    return(Value)

def GetArrayBounds(DataArray):
    MinValue = "NAN"
    MaxValue = "NAN"
    for i in range(len(DataArray)):
        Value = DataArray[i]

        if(MinValue == "NAN"):
            MinValue = Value
            MaxValue = Value
        else:
            if(Value < MinValue):
                MinValue = Value
            else:
                if(Value > MaxValue):
                    MaxValue = Value

    return([MinValue, MaxValue])

def CopyStringDataToList(String):
    #Using a whitelist because there are more charachters we want to ignore than we want to use
    ValidCharachters = ["0","1","2","3","4","5","6","7","8","9",".","-"]

    DataList = []
    CurrentString = ""
    ValuesFound = 0
    for Count in range(len(String)):
        Char = String[Count]
        
        if Char == ",":
            StringValue = float(CurrentString)
            DataList.append(StringValue)
            CurrentString = ""
            ValuesFound+=1
        else:

            # Variable in Array/List
            #Same as C#:
            #Array/List.contains(Variable)

            if Char in ValidCharachters:
                CurrentString+=str(Char)

    DataList.append(float(CurrentString))

    if(ValuesFound < 3):
        DataList.append(-1)

    return(DataList)

def CustomChiSqrInputFunction(Params, DataX, DataY, DataERROR, Priors):

    ParamaterValues = Params.valuesdict()

    FittingTransityFunctionParams = batman.TransitParams()

    FittingTransityFunctionParams.t0 = ParamaterValues["t0"]                        #time of inferior conjunction
    FittingTransityFunctionParams.per = ParamaterValues["per"]                       #orbital period
    FittingTransityFunctionParams.rp = ParamaterValues["rp"]                       #planet radius (in units of stellar radii)
    FittingTransityFunctionParams.a = ParamaterValues["a"]                        #semi-major axis (in units of stellar radii)
    FittingTransityFunctionParams.inc = ParamaterValues["inc"]                      #orbital inclination (in degrees)
    FittingTransityFunctionParams.ecc = ParamaterValues["ecc"]                       #eccentricity
    FittingTransityFunctionParams.w = ParamaterValues["w"]                        #longitude of periastron (in degrees)
    FittingTransityFunctionParams.limb_dark = "quadratic"        #limb darkening model
    FittingTransityFunctionParams.u = [ParamaterValues["u1"], ParamaterValues["u2"]]      #limb darkening coefficients [u1, u2]

    '''
    print(
       str(FittingTransityFunctionParams.t0) + "\n" + 
        str(FittingTransityFunctionParams.per) + "\n" + 
        str(FittingTransityFunctionParams.rp) + "\n" + 
        str(FittingTransityFunctionParams.a) + "\n" + 
        str(FittingTransityFunctionParams.inc) + "\n" + 
        str(FittingTransityFunctionParams.ecc) + "\n" + 
        str(FittingTransityFunctionParams.w) + "\n" + 
        str(FittingTransityFunctionParams.u[0]) + "\n" + 
        str(FittingTransityFunctionParams.u[1]) + "\n" +
        " "
        )
        '''

    #Will try to minimize returned value
    return(CalculateChiSqr(DataX, DataY, FittingTransityFunctionParams, DataERROR, Priors))

def LmfitInputFunction(x,t0,per,rp,a,inc,ecc,w,u1,u2):

    FittingTransityFunctionParams = batman.TransitParams()

    FittingTransityFunctionParams.t0 = t0                        #time of inferior conjunction
    FittingTransityFunctionParams.per =per                       #orbital period
    FittingTransityFunctionParams.rp = rp                       #planet radius (in units of stellar radii)
    FittingTransityFunctionParams.a = a                        #semi-major axis (in units of stellar radii)
    FittingTransityFunctionParams.inc = inc                      #orbital inclination (in degrees)
    FittingTransityFunctionParams.ecc = ecc                       #eccentricity
    FittingTransityFunctionParams.w = w                        #longitude of periastron (in degrees)
    FittingTransityFunctionParams.limb_dark = "quadratic"        #limb darkening model
    FittingTransityFunctionParams.u = [u1, u2]      #limb darkening coefficients [u1, u2]

    m = batman.TransitModel(FittingTransityFunctionParams, x)
    flux = m.light_curve(FittingTransityFunctionParams)

    return (flux)

def OptimizeFunctionParameters(DataX, DataY, DataERROR, Priors, UseLmfit, StartingParameters):
    WeightedDataErrorArray = []

    DataIncludedErrorBars = GetArrayIsNotNone(DataERROR)

    WeightedDataErrorArray

    if(DataIncludedErrorBars):
        #Remove zeros from array, an error value of '0' leeds to calculation errors and is improbable
        WeightedDataErrorArray = (1.0/ReplaceZerosInArrayWithLowestValue(DataERROR))

    else:
        WeightedDataErrorArray = (1.0+0*DataX)    


    Bounds = GetArrayBounds(DataX)
    MinX = Bounds[0]
    MaxX = Bounds[1]



    InputParams = lmfit.Parameters()
    if(UseLmfit and not StartingParameters == None):
        #Lmfit version
        InputParams.add("t0", value=StartingParameters.t0, min = MinX, max = MaxX) #Max?
        InputParams.add("per", value=StartingParameters.per, min = 0.0, max = MaxX)
        InputParams.add("rp", value=StartingParameters.rp, min = 0, max = 10.0)
        InputParams.add("a", value=StartingParameters.a, min = 1.0, max = 90) #What should Max Bound be?
        InputParams.add("inc", value=StartingParameters.inc, min = 60, max = 90)
        InputParams.add("ecc", value=StartingParameters.ecc, min = 0.0, max = 1.0)
        InputParams.add("w", value=StartingParameters.w, min = 0.0, max = 360.0)
        InputParams.add("u1", value=StartingParameters.u[0], min = -1.0, max = 1.0)
        InputParams.add("u2", value=StartingParameters.u[1], min = -1.0, max = 1.0)
    else:
        #Minimize version

        InitialValue_t0 = 0.0
        InitialValue_per = 0.0
        InitialValue_rp = 0.0
        InitialValue_a = 0.0
        InitialValue_inc = 0.0
        InitialValue_ecc = 0.0
        InitialValue_w = 0.0
        InitialValue_u1 = 0.0
        InitialValue_u2 = 0.0

        if(not Priors == None):
            if(not Priors[0] == None):
                InitialValue_t0 = Priors[0][0]
            if(not Priors[1] == None):
                InitialValue_per = Priors[1][0]
            if(not Priors[2] == None):
                InitialValue_rp = Priors[2][0]
            if(not Priors[3] == None):
                InitialValue_a = Priors[3][0]
            if(not Priors[4] == None):
                InitialValue_inc = Priors[4][0]
            if(not Priors[5] == None):
                InitialValue_ecc = Priors[5][0]
            if(not Priors[6] == None):
                InitialValue_w = Priors[6][0]
            if(not Priors[7] == None):
                InitialValue_u1 = Priors[7][0]
            if(not Priors[8] == None):
                InitialValue_u2 = Priors[8][0]

        InputParams.add("t0", value=InitialValue_t0, min = MinX, max = MaxX) #Max?
        InputParams.add("per", value=InitialValue_per, min = 0.0, max = MaxX)
        InputParams.add("rp", value=InitialValue_rp, min = 0, max = 10.0)
        InputParams.add("a", value=InitialValue_a, min = 1.0, max = 90) #What should Max Bound be?
        InputParams.add("inc", value=InitialValue_inc, min = 60, max = 90)
        InputParams.add("ecc", value=InitialValue_ecc, min = 0.0, max = 1.0)
        InputParams.add("w", value=InitialValue_w, min = 0.0, max = 360.0)
        InputParams.add("u1", value=InitialValue_u1, min = -1.0, max = 1.0)
        InputParams.add("u2", value=InitialValue_u2, min = -1.0, max = 1.0)
        

    LmfitOptimizedFunction = None

    if(not UseLmfit):
        LmfitOptimizedFunction = lmfit.minimize(
        
            CustomChiSqrInputFunction,

            InputParams,
            args=(DataX,DataY,DataERROR, Priors), 
            method = "nelder",
            max_nfev=None)
    else:
        LmfitOptimizedFunction = Model(LmfitInputFunction).fit(
            DataY, x=DataX,

            t0=Parameter("t0", value=StartingParameters.t0, min = MinX, max = MaxX),
            per=Parameter("per", value=StartingParameters.per, min = 0.0, max = MaxX),
            rp=Parameter("rp", value=StartingParameters.rp, min = 0, max = 10.0),
            a=Parameter("a", value=StartingParameters.a, min = 1.0, max = 90),
            inc=Parameter("inc", value=StartingParameters.inc, min = 60, max = 90),
            ecc=Parameter("ecc", value=StartingParameters.ecc, min = 0.0, max = 1.0),
            w=Parameter("w", value=StartingParameters.w, min = 0.0, max = 360.0),
            u1=Parameter("u1", value=StartingParameters.u[0], min = -1.0, max = 1.0),
            u2=Parameter("u2", value=StartingParameters.u[1], min = -1.0, max = 1.0),

            weights=WeightedDataErrorArray,
            max_nfev=None
            )

        #Weight Implementation According To Documentation:
        '''
        weights
        numpy.ndarray (or None) of weighting values to be used in fit.
        If not None, it will be used as a multiplicative factor of the residual array, so that weights*(data - fit) is minimized in the least-squares sense.
        '''


    return(LmfitOptimizedFunction)



#Main function
def RunOptimizationOnDataInputFile(Priors):

    print("Running. Please wait...")

    DataPoints = np.array([[0,0,0]])

    #Clear initialized array value, because don't know how to initialize empty array with bounds
    DataPoints = np.delete(DataPoints,0, 0)


    FileName = "Data"
    FileType = "txt"


    DataIncludedErrorBars = False
    Lines = []

    #r - read
    #w - write
    #a - append data

    with open(FileName+"."+FileType, "r") as File:
        Lines = File.readlines(-1)



    for Count in range(len(Lines)):
        #Prevents empty lines from being run through string converter, which if returned would lead to adding an empty array element, or an error
        if(Lines[Count] != "\n"):
            DataPoints = np.append(DataPoints, [CopyStringDataToList(Lines[Count])],0)


    ListDataX = []
    ListDataY = []
    ListDataERROR = []


    NumberOfDataPoints = len(DataPoints)


    for Count in range(NumberOfDataPoints):
        ListDataX.append(DataPoints[Count][0])

    DataX = np.array(ListDataX)

    for Count in range(NumberOfDataPoints):
        ListDataY.append(DataPoints[Count][1])

    DataY = np.array(ListDataY)

    for Count in range(NumberOfDataPoints):
        ListDataERROR.append(DataPoints[Count][2])

    DataERROR = np.array(ListDataERROR)

    if(DataERROR[0] == -1):
        DataIncludedErrorBars = False
    else:
        DataIncludedErrorBars = True


    #Get data bounds
    #Used in some parameter bounds

    Bounds = GetArrayBounds(DataX)
    MinX = Bounds[0]
    MaxX = Bounds[1]

    Bounds = GetArrayBounds(DataY)
    MinY = Bounds[0]
    MaxY = Bounds[1]


    if(TestMode):
        if(TestModeRenderPoints):
            matplot.errorbar(DataX, DataY, fmt ="o", markersize = DataPointRenderSize)
        matplot.show()
        time.sleep(1000)

    '''
    FittingTransityFunctionParams = batman.TransitParams()

    #Initialize params with placeholder values, these will be overridden imediately and will not be used
    #Need to be realistic in order to generate a reasonable step size, uses priors if supplied
    FittingTransityFunctionParams.t0 = 0.96
    FittingTransityFunctionParams.per = 3.1
    FittingTransityFunctionParams.rp = 0.118
    FittingTransityFunctionParams.a = 15.0
    FittingTransityFunctionParams.inc = 87.0
    FittingTransityFunctionParams.ecc = 0.0
    FittingTransityFunctionParams.w = 90.0
    FittingTransityFunctionParams.limb_dark = "quadratic"
    FittingTransityFunctionParams.u = [0.5, 0.1]

    if(not Priors == None):
        if(not Priors[0][0] == None):
            FittingTransityFunctionParams.t0 = Priors[0][0]
        if(not Priors[1][0] == None):
            FittingTransityFunctionParams.per = Priors[1][0]
        if(not Priors[2][0] == None):
            FittingTransityFunctionParams.rp = Priors[2][0]
        if(not Priors[3][0] == None):
            FittingTransityFunctionParams.a = Priors[3][0]
        if(not Priors[4][0] == None):
            FittingTransityFunctionParams.inc = Priors[4][0]
        if(not Priors[5][0] == None):
            FittingTransityFunctionParams.ecc = Priors[5][0]
        if(not Priors[6][0] == None):
            FittingTransityFunctionParams.w = Priors[6][0]
        UParameters = [FittingTransityFunctionParams.u[0],FittingTransityFunctionParams.u[1]]
        if(not Priors[7][0] == None):
            UParameters[0] = Priors[7][0]
        if(not Priors[7][0] == None):
            UParameters[1] = Priors[8][0]
        FittingTransityFunctionParams.u = UParameters

    '''

    #FitModelFunction = batman.TransitModel(FittingTransityFunctionParams, DataXArray)


    #First optimization attempt, used to get error values
    LmfitOptimizedFunction = OptimizeFunctionParameters(DataX, DataY, None, Priors, False, None)

    global CountTime
    CountTime = False

    #Extract parameters used
    OptimizedParams = ExtractParametersFromFittedFunction(LmfitOptimizedFunction)
    #Have to manually assign params instead of using 'LmfitOptimizedFunction.params' because 'limb_dark' is not assigned by the function

    #Generate function based on extracted parameters
    FirstOptimizedFunction = batman.TransitModel(OptimizedParams, DataX)

    #Calculate error from diference between first attempt created function, and given values

    #Should abs be used?
    DataERROR = (DataY*0 + np.std(DataY-FirstOptimizedFunction.light_curve(OptimizedParams)))
    #CheckChiSqr funciton expects data error as an array, this allows compatibilty with lmfit.fit instead of lmfit.minimize

    #data-model
    #stder(data-model)

    #Run second time, using newly calculated error values

    SecondOptimizedFunction = OptimizeFunctionParameters(DataX, DataY, DataERROR, Priors, False, None)

    #Extract paramaters from optimized function
    OptimizedParams = ExtractParametersFromFittedFunction(SecondOptimizedFunction)

    NewDataValues = RemoveOutliersFromDataSet(DataX, DataY, OptimizedParams)

    DataX = NewDataValues[0]
    DataY = NewDataValues[1]
    NumberOfDataPoints = NewDataValues[2]
    IndexesRemoved = NewDataValues[3]
    if (len(IndexesRemoved) > 0):
        DataERROR = np.delete(DataERROR,IndexesRemoved)

    #Run third time, this time using lmfit with newly generated starting parameter values
    #And having removed outliers
    ThirdOptimizedFunction = OptimizeFunctionParameters(DataX, DataY, DataERROR, Priors, False, None)


    #CheckTime(True)
    LmfitOptimizedFunction = ThirdOptimizedFunction
    OptimizedParams = ExtractParametersFromFittedFunction(ThirdOptimizedFunction)
    #CheckTime(False)

    DataIncludedErrorBars = True

    #Debug Fit Report
    print(fit_report(LmfitOptimizedFunction))
    print("\n")

    #Display points with error bars
    if(DataIncludedErrorBars):
        matplot.errorbar(DataX, DataY, yerr = DataERROR, fmt ="o", markersize = DataPointRenderSize)
    else:
        matplot.errorbar(DataX, DataY, fmt ="o", markersize = DataPointRenderSize)

    print("-OPTIMIZED PARAMETERS-")
    parameter_names = 't0', 'per', 'rp', 'a', 'inc', 'ecc', 'w', 'u1', 'u2'
    for name in parameter_names:
        print(name + ' : ' + str(LmfitOptimizedFunction.params[name].value))
    print("\n")

    print("-OPTIMIZED PARAMETER UNCERTAINTY VALUES-")
    for name in parameter_names:
        print(name + ' : ' + str(LmfitOptimizedFunction.params[name].stderr))
    print("\n")

    OptimizedParams = batman.TransitParams()
    OptimizedParams.t0 = LmfitOptimizedFunction.params["t0"].value
    OptimizedParams.per = LmfitOptimizedFunction.params["per"].value
    OptimizedParams.rp = LmfitOptimizedFunction.params["rp"].value
    OptimizedParams.a = LmfitOptimizedFunction.params["a"].value
    OptimizedParams.inc = LmfitOptimizedFunction.params["inc"].value
    OptimizedParams.ecc = LmfitOptimizedFunction.params["ecc"].value
    OptimizedParams.w = LmfitOptimizedFunction.params["w"].value
    OptimizedParams.u = [LmfitOptimizedFunction.params["u1"].value, LmfitOptimizedFunction.params["u2"].value]
    OptimizedParams.limb_dark = "quadratic"



    CheckedOptimizedChiSqr = CalculateChiSqr(DataX,DataY,OptimizedParams, DataERROR, Priors)
    



    SamplePoints = np.linspace(MinX,MaxX,10000)
    m = batman.TransitModel(OptimizedParams, SamplePoints)
    flux = m.light_curve(OptimizedParams)
    matplot.plot(SamplePoints,flux, "-", label="Optimized Function")


    #Debug Logging

    StringData= ""

    for i in range(len(DataX)):
        StringData+=(str(DataY[i]) + "\n" )

    with open("Output.txt", "w") as File:
        Lines = File.writelines(StringData)



    print("\n--- Checked Chi Sqr ---")
    print("ChiSqr : " + str(CheckedOptimizedChiSqr))
    print("Number Of Data Points : " + str(NumberOfDataPoints))
    print(str(1.0/NumberOfDataPoints * (CheckedOptimizedChiSqr)))
    print(str(1.0/NumberOfDataPoints * (CheckedOptimizedChiSqr)-1.0000451473151235))

    #Fixed "χ2" rendering issue
    BestChi = "Optimized χ2 : " + str(round(CheckedOptimizedChiSqr,2))

    ChiAnchoredTextBox = AnchoredText(BestChi
                                    , loc=4, pad=0.5)
    matplot.setp(ChiAnchoredTextBox.patch, facecolor="Orange", alpha=0.5)
    matplot.gca().add_artist(ChiAnchoredTextBox)
    matplot.legend(loc=2, borderaxespad=0)



    EndTimeRecording()

    #Display plotted data
    matplot.show()

def RemoveOutliersFromDataSet(DataX, DataY, Parameters):

    TestMode = True

    NewDataX = DataX
    NewDataY = DataY
    NewNumberOfDataPoints = -1

    TransitModel = batman.TransitModel(Parameters, DataX)
    LightCurve = TransitModel.light_curve(Parameters)

    StandardDeviation = (np.std(DataY-LightCurve))



    Differences = abs((DataY-LightCurve))

    MeanDifference = Differences.sum()/len(Differences)

    DiferenceBounds = GetArrayBounds(Differences)

    DifferenceMin = DiferenceBounds[0]
    DifferenceMax = DiferenceBounds[1]

    DiferenceLimitMax = 4
    #Conservative 6
    #Reasonable 4
    #High reduction 3.5

    Colors = []

    IndexesToRemove = []

    Count = 0
    for Difference in Differences:
        #If diference between the input YValue and the models YValue is greater than the mean diference of all data points * 'DiferenceLimitMax'
        #Then remove the data point
        if((Difference-StandardDeviation) > (MeanDifference*DiferenceLimitMax)):
            IndexesToRemove.append(Count)

            if(TestMode):
                Colors.append((0.5,0,0,0.75))
        else:
            if(TestMode):
                Colors.append((0, 0.5,0,0.75))
        Count+=1

    IndexesToRemove =   np.array(IndexesToRemove)

    if(len(IndexesToRemove) > 0):
        NewDataX = np.delete(NewDataX, IndexesToRemove)
        NewDataY = np.delete(NewDataY, IndexesToRemove)

    NewNumberOfDataPoints = len(DataY)

    if(TestMode):
        #Will plot all input data points, green if not removed, red if removed
        CheckTime(True)

        PlotType = 1

        if(PlotType == 0):
            #Difference Value
            #matplot.scatter(DataX, Differences, color = Colors)

            #X,Y Value
            matplot.scatter(DataX , DataY, color = Colors)
        else:
            if(PlotType == 1):
                #GradientColors = (colo)
                HeatMapColors = []
                Count = 0
                for Difference in Differences:
                    if(not( (Difference-StandardDeviation) > (MeanDifference*DiferenceLimitMax))):
                        HeatMapColors.append((1,0,0,Clamp((1.0/(MeanDifference*DiferenceLimitMax) * (Difference-StandardDeviation)),0.0,1.0)))
                    else:
                        HeatMapColors.append((1.0,0.647,0,0.9))

                matplot.scatter(DataX , DataY, color=HeatMapColors)





        CheckTime(False)

        EndTimeRecording()

        matplot.show()

    return(NewDataX,NewDataY,NewNumberOfDataPoints, IndexesToRemove)

def ExtractParametersFromFittedFunction(Function):

    ExtractedParameters = batman.TransitParams()
    ExtractedParameters.t0 = Function.params["t0"].value
    ExtractedParameters.per = Function.params["per"].value
    ExtractedParameters.rp = Function.params["rp"].value
    ExtractedParameters.a = Function.params["a"].value
    ExtractedParameters.inc = Function.params["inc"].value
    ExtractedParameters.ecc = Function.params["ecc"].value
    ExtractedParameters.w = Function.params["w"].value
    ExtractedParameters.u = [Function.params["u1"].value, Function.params["u2"].value]
    ExtractedParameters.limb_dark = "quadratic"

    return(ExtractedParameters)

def GetArrayIsNotNone(InputArray):
    return(InputArray is not None)

def EndTimeRecording():
    print("Block Percent Of Total Time : " + str(100.0/(time.time()-ProgramStartTime) * StartTime))
    print("Total Time : " + str(time.time()-ProgramStartTime))

RunOptimizationOnDataInputFile(Priors)