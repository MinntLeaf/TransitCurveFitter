from cmath import exp, pi, sqrt
from datetime import time
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


#Only used for testing
import random

#Profiling library
import cProfile



#If supplied, will be used as initial fitting parameters
#First array element is value, second is certainty
#If either value is missing from an element, both elements will be ignored
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
#NOTE: This program assumes only 2 limb-darkening values, more are possible with modification

#Visual modifier only
DataPointRenderSize = 1
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
    if(IsStart):
        if(len(StartTimes) <= ID):
            StartTimes.append(time.time())
        else:
            StartTimes[ID] = time.time()
    else:
        #Is end
        print("Block Time ID :",str(ID),": Time :", str(time.time()-StartTimes[ID]));
        EndTimes.append(time.time())

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
CountTime = False
def CalculateChiSqr(DataX, DataY, Params, DataERROR, Priors):

    if(CountTime):
        global GlobalCount
        GlobalCount+=1
        print(GlobalCount)

    

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

    if(Priors is not None):
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
    #Should replace with external package

    if(Value < Min):
        return(Min)
    if(Value > Max):
        return(Max)
    return(Value)

def GetArrayBounds(DataArray):
    DataArray = np.array(DataArray, copy=False)
    return([DataArray.min(), DataArray.max()])

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


    DataIncludedErrorBars = GetArrayIsNotNone(DataERROR)

    WeightedDataErrorArray = None

    if(DataIncludedErrorBars):
        #Remove zeros from array, an error value of '0' leeds to calculation errors and is improbable
        WeightedDataErrorArray = (1.0/ReplaceZerosInArrayWithLowestValue(DataERROR))


    #else:
    #    WeightedDataErrorArray = (1.0+0*DataX)    
    #   Lmfit should interpret a passed value of 'None' as a lack of weights value, which should weight them all equally


    Bounds = GetArrayBounds(DataX)
    MinX = Bounds[0]
    MaxX = Bounds[1]



    InputParams = lmfit.Parameters()
    if(UseLmfit and StartingParameters is not None):
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

        if(Priors is not None):
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
        

    OptimizedFunctionToReturn = None

    if(not UseLmfit):
        OptimizedFunctionToReturn = lmfit.minimize(
        
            CustomChiSqrInputFunction,

            InputParams,
            args=(DataX,DataY,DataERROR, Priors), 
            method = "nelder",
            max_nfev=None)
    else:
        OptimizedFunctionToReturn = Model(LmfitInputFunction).fit(
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


    return(OptimizedFunctionToReturn)



#Main function
def RunOptimizationOnDataInputFile(Priors):

    print("Running. Please wait...\n")

    DataPoints = np.array([[0,0,0]])

    #Clear initialized array value, because don't know how to initialize empty array with bounds
    DataPoints = np.delete(DataPoints,0, 0)


    FileName = "Data"
    FileType = "txt"

    # Assume CSV format with no header row:
    datafile = np.loadtxt(FileName+"."+FileType, delimiter=',')
    if datafile.shape[1]==2:
        DataIncludedErrorBars = False
        DataX = datafile[:,0]
        DataY = datafile[:,1]
        DataERROR = -np.ones(DataX.shape)
    elif datafile.shape[1]==3:
        DataIncludedErrorBars = True
        DataX = datafile[:,0]
        DataY = datafile[:,1]
        DataERROR = datafile[:,2]

          
    #Get data bounds
    #Used in some parameter bounds

    Bounds = GetArrayBounds(DataX)
    MinX = Bounds[0]
    MaxX = Bounds[1]

    #Not currently used
    Bounds = GetArrayBounds(DataY)
    MinY = Bounds[0]
    MaxY = Bounds[1]


    if(TestMode):
        if(TestModeRenderPoints):
            matplot.errorbar(DataX, DataY, fmt ="o", markersize = DataPointRenderSize)
        matplot.show()
        time.sleep(1000)

    CheckTime(0,True)
    # - 5.2s
    #First optimization attempt, used to get error values
    OptimizedFunction = OptimizeFunctionParameters(DataX, DataY, None, Priors, False, None)
    CheckTime(0,False)

    #Extract parameters used
    OptimizedParams = ExtractParametersFromFittedFunction(OptimizedFunction)

    #Generate function based on extracted parameters
    FirstOptimizedFunction = batman.TransitModel(OptimizedParams, DataX)

    #Calculate error from diference between first attempt created function, and given values
    if(not DataIncludedErrorBars):
        #I think "abs" should not be used here, the square of the values is being used instead. Not sure why abs affects the result in this case, but it does.
        DataERROR = (DataY*0 + np.std((DataY-FirstOptimizedFunction.light_curve(OptimizedParams))))
        #CheckChiSqr function expects data error as an array, this allows compatibilty with lmfit.fit instead of lmfit.minimize

        #Debug logging
        #print(np.std((DataY-FirstOptimizedFunction.light_curve(OptimizedParams))))



    #data-model
    #stder(data-model)

    #Run second time, using newly calculated error values

    CheckTime(1,True)
    # - 5.2s
    #Why?
    #Including error values appears to increase run time signoficantly
    SecondOptimizedFunction = OptimizeFunctionParameters(DataX, DataY, DataERROR, Priors, False, None)
    CheckTime(1,False)

    #Extract paramaters from optimized function
    OptimizedParams = ExtractParametersFromFittedFunction(SecondOptimizedFunction)


    #Remove Outlier values
    NewDataValues = RemoveOutliersFromDataSet(DataX, DataY, OptimizedParams)

    DataX = NewDataValues[0]
    DataY = NewDataValues[1]
    NumberOfDataPoints = NewDataValues[2]
    IndexesRemoved = NewDataValues[3]
    if (len(IndexesRemoved) > 0):
        DataERROR = np.delete(DataERROR,IndexesRemoved)



    #Recalculate error
    #If data did not include errors, and outliers have just been removed
    if(not DataIncludedErrorBars and (len(IndexesRemoved) > 0)):

        #Have to remove removed values from returned light values, because can't calculate std of diference betweeen arrays, when those arrays are of diferent lengths
        UpdatedLightValues = np.delete(FirstOptimizedFunction.light_curve(OptimizedParams),IndexesRemoved)

        #Recalcualte error values with the outlier values removed
        DataERROR = (DataY*0 + np.std((DataY-UpdatedLightValues)))


    #Run third time, this time having removed outliers
    CheckTime(2,True)
    # - 3.5s
    #Why shorter than the other OptimizeFunctionParameters() calls?
    ThirdOptimizedFunction = OptimizeFunctionParameters(DataX, DataY, DataERROR, Priors, False, None)
    CheckTime(2,False)


    FinalOptimizedFunction = ThirdOptimizedFunction
    OptimizedParams = ExtractParametersFromFittedFunction(ThirdOptimizedFunction)
    

    #From this point forward, "DataIncludedErrorBars" will no longer be used to decide if Error values need to be calculated
    #It just means wether to use the current Error np.array for Chi calcualtions and wether to render points with error bars
    DataIncludedErrorBars = True



    #Debug Fit Report
    print(fit_report(FinalOptimizedFunction))
    print("\n")

    #Display points with error bars
    if(DataIncludedErrorBars):
        #Disabled for now
        #Concerned the 'yerr' value being set to the DataERROR value is not an acurate representation of the points error

        #matplot.errorbar(DataX, DataY, yerr = DataERROR, fmt ="o", markersize = DataPointRenderSize)
        matplot.errorbar(DataX, DataY, fmt ="o", markersize = DataPointRenderSize)
    else:
        matplot.errorbar(DataX, DataY, fmt ="o", markersize = DataPointRenderSize)

    print("-OPTIMIZED PARAMETERS-")
    parameter_names = 't0', 'per', 'rp', 'a', 'inc', 'ecc', 'w', 'u1', 'u2'
    for name in parameter_names:
        print(name + ' : ' + str(FinalOptimizedFunction.params[name].value))
    print("\n")

    print("-OPTIMIZED PARAMETER UNCERTAINTY VALUES-")
    for name in parameter_names:
        print(name + ' : ' + str(FinalOptimizedFunction.params[name].stderr))
    print("\n")


    CheckedOptimizedChiSqr = CalculateChiSqr(DataX,DataY,OptimizedParams, DataERROR, Priors)
    


    #Rendering only, uses more sample points than input x-values
    SamplePoints = np.linspace(MinX,MaxX,10000)
    m = batman.TransitModel(OptimizedParams, SamplePoints)
    flux = m.light_curve(OptimizedParams)
    matplot.plot(SamplePoints,flux, "-", label="Optimized Function")


    #Debug Logging START

    StringData= ""

    DebugFlux =  batman.TransitModel(OptimizedParams, DataX).light_curve(OptimizedParams)

    for i in range(len(DataX)):
        StringData+=(str(DebugFlux[i]) + "\n" )
    '''
    print("-----")

    print(str(np.std((DataY-DebugFlux))))

    print("-----")
    '''
    with open("Output.txt", "w") as File:
        Lines = File.writelines(StringData)



    print("\n--- Checked Chi Sqr ---")
    print("ChiSqr : " + str(CheckedOptimizedChiSqr))
    print("Number Of Data Points : " + str(NumberOfDataPoints))
    #The value below should be close to '1'
    print("ChiSqr / # Data Points : " + str(CheckedOptimizedChiSqr/NumberOfDataPoints))

    #Debug Logging END


    #Fixed "χ2" rendering issue
    BestChi = "Optimized χ2 : " + str(round(CheckedOptimizedChiSqr,2))

    #Text box setup
    ChiAnchoredTextBox = AnchoredText(BestChi
                                    , loc=4, pad=0.5)
    matplot.setp(ChiAnchoredTextBox.patch, facecolor="Orange", alpha=0.5)
    matplot.gca().add_artist(ChiAnchoredTextBox)
    matplot.legend(loc=2, borderaxespad=0)



    EndTimeRecording()

    print("\nCompleted")

    #Display plotted data
    #Will 'end' program, code after this function is called will not be run (this behavior can be changed)
    matplot.show()

def RemoveOutliersFromDataSet(DataX, DataY, Parameters):

    TestMode = False
    #Only valid if 'TestMode' is active

    #will not halt further execution of the program, instead overlaying the scatter avlues or heat map underneath the final graph
    OverlayMode = False

    #Show limits are allowed between
    HighlightBoundsMode = False

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

    #Should not be set below 1
    DiferenceLimitMultiplier = 2
    #Conservative 4
    #Reasonable 2
    #High reduction 1

    Colors = []

    IndexesToRemove = []

    MaxDifferenceAllowed = MeanDifference*2*DiferenceLimitMultiplier

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

        if((Difference) > MaxDifferenceAllowed):
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

    NewNumberOfDataPoints = len(NewDataY)

    #Debug visuals
    if(TestMode):

        PlotType = 1

        if(PlotType == 0):
            #Will plot all input data points, green if not removed, red if removed


            #Difference Value
            #matplot.scatter(DataX, Differences, color = Colors)

            #X,Y Value
            matplot.scatter(DataX , DataY, color = Colors)

        else:
            if(PlotType == 1):
                HeatMapColors = []
                Count = 0
                for Difference in Differences:
                    if(not( (Difference) > MaxDifferenceAllowed)):
                        HeatMapColors.append((1,0,0,Clamp((1.0/MaxDifferenceAllowed * (Difference)),0.0,1.0)))
                    else:
                        HeatMapColors.append((1.0,0.647,0,0.9))

                matplot.scatter(DataX ,DataY, color=HeatMapColors, s = 8)

        
        if(HighlightBoundsMode):
            XBounds = GetArrayBounds(DataX)

            SamplePoints = np.linspace(XBounds[0],XBounds[1],10000)
            m = batman.TransitModel(Parameters, SamplePoints)
            LightCurve = m.light_curve(Parameters)

            matplot.plot(SamplePoints,LightCurve, "-", color="yellow")
            matplot.plot(SamplePoints,LightCurve + MaxDifferenceAllowed, "-", color="green")
            matplot.plot(SamplePoints,LightCurve - MaxDifferenceAllowed, "-", color="green")
                
                




        if(not OverlayMode):
            matplot.show()

        #Results Interpretation:

        #Green lines are bounds in which points are allowd, points that do not fall between these lines will be removed
        #Yellow line is the fitted line, it's y values are what the points are being compared to
        #Orange points have been removed
        #Red points have not been removed, the transparency of their color indicates how close they are to being removed. Dark red ones are close to the limit, almost clear ones are close to their expected values.
        #Note ^ points overlayed on top of eachother will stack their transparencies, resulting in dark coolors even if they are not close to the border. Zoom in on the graph if you wish to acurately see the colors of points that are overlaping one another.

    if(NewNumberOfDataPoints + len(IndexesToRemove) != len(DataY)):
        #New number of datapoints value given is not equal to the original number minus those removed
        #This means there is an issue with the 'RemoveOutliersFromDataSet' function
        print("ERROR : RemoveOutliersFromDataSet() returned an improper value")

        #This check is here because small diferences in the actual number of data values [len(DataX)] and the recorded number of data points [NumberOfDataPoints] values can easilly go unoticed and lead to inacurate ChiSqr values down the line

        #Stop, definitely better way to do this
        time.sleep(99999)

    return(NewDataX,NewDataY,NewNumberOfDataPoints, IndexesToRemove)

def ExtractParametersFromFittedFunction(Function):

    #Have to manually assign params instead of using 'OptimizedFunction.params' because 'limb_dark' is not assigned by the function

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
    #Should (type(InputArray) == "NoneType") be used instead?
    return(InputArray is not None)

def EndTimeRecording():

    global StartTimes
    global EndTimes
    global ProgramStartTime

    if(len(StartTimes) > 0):
        print("\n----- Block Percents Of Total Time -----")
        for i in range(len(StartTimes)):
            #Block times are only valid if ID's were referenced only once
            print("Block Percent Of Total Time : " + str(100.0/(time.time()-ProgramStartTime) * (EndTimes[i]-StartTimes[i])))
        
    print("\n----- Total Time -----")
    print("Total Time : " + str(time.time()-ProgramStartTime))


RunOptimizationOnDataInputFile(Priors)
