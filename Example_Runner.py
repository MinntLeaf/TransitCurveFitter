import Transit_Curve_Fitter
import numpy as np
import matplotlib.pyplot as matplot
from matplotlib.offsetbox import AnchoredText
import batman

#Helper Functions
#----------------------

def ConvertParamatersToFluxValues(TimeVal, Params):
    #Convert paramaters and time values into model flux values at the given time values
    #X values do not need to be the ones used to generate the model, flux values can be generated from any time value that is valid within the model

    BatmanParams = ConvertFitParametersToTransitParameters(Params)

    TransitModel = batman.TransitModel(BatmanParams, TimeVal)
    LightCurve = TransitModel.light_curve(BatmanParams)
    LightCurve = ApplyPolyMultiplier(TimeVal, LightCurve, Params)
    return(LightCurve)

def ApplyPolyMultiplier(TimeVal, FluxVal,  Params):
    #Applies the polynomial parameter to a models flux values, if the polynomial modifier is 'None' it will return unmodified flux values.

    if(PolynomialOrder is not None):
        PolyVals = []
        for PolyVal in range(0,PolynomialOrder+1):
            PolyVals.append(Results["PolyVal" + str(PolyVal)][0])

        return(FluxVal * np.polyval(PolyVals,TimeVal))
    return(FluxVal)

def ConvertFitParametersToTransitParameters(InputParam):
    #Converts a dictionary of model paramaters into 'batman.TransitParams()' to be used to generate a transit model with the 'batman' package.
    #The Transit_Curve_Fitter script will output a dictionary of values so this can be used to convert them into usable batman paramteters.

    FittingTransitFunctionParams = batman.TransitParams()

    FittingTransitFunctionParams.t0 = InputParam["t0"][0]
    FittingTransitFunctionParams.per = InputParam["per"][0]
    FittingTransitFunctionParams.rp = InputParam["rp"][0]
    FittingTransitFunctionParams.a = InputParam["a"][0]
    FittingTransitFunctionParams.inc = InputParam["inc"][0]
    FittingTransitFunctionParams.ecc = InputParam["ecc"][0]
    FittingTransitFunctionParams.w = InputParam["w"][0]
    FittingTransitFunctionParams.limb_dark = "quadratic"
    FittingTransitFunctionParams.u = [InputParam["u1"][0], InputParam["u2"][0]]

    return (FittingTransitFunctionParams)

#----------------------





#Example Usage
#----------------------

FileName = "Example_Data"
FileType = "txt"

#Copy in data, assuming CSV format with no header row:

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

#Create a new fit data class
NewFitData = Transit_Curve_Fitter.FitData()

#Fill out the properties of the 'FitData' class with Time, Flux, and Uncertainty values
#If Uncertainty is input as None the fitter will estimate it's own uncertainty values
NewFitData.Time = DataX
NewFitData.Flux = DataY
NewFitData.Uncertainty = None

#Fill out initial fitting values and uncertainty values
#These values will be used by the fitter as starting parameters
#These values will be used as priors if uncertainty is given, and affect the chivalues if the final result, if uncertainty is not given they will only be used as starting parameters
NewFitData.t0 = [0.99, 0.01]  #t0
NewFitData.per = [3.14159265, 0.01]  #per
NewFitData.rp =   [0.118, 0.01]  #rp
NewFitData.a =  [8.0, 1]  #a
NewFitData.inc =  [83.7, 1]  #inc
NewFitData.ecc = [0.0, 1]  #ecc
NewFitData.w =  [1.0, 1]  #w
NewFitData.u1 = [-0.9, 1]  #u1
NewFitData.u2 =  [-0.9, 1]  #u2
NewFitData.PolynomialOrder = 2 #Polynomial Order [Set to None if no modifier is needed, higher values will reduce fitting accuracy]

#Run fit on the filled out NewFitData class and store results as 'Results'
Results = Transit_Curve_Fitter.TransitCurveFitter.FitTransit(NewFitData)

#Print out 'Results' values
print("")
for Key,Value in Results.items():
    print(Key,":",str(Value))
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
#PolynomialOrder : Value (Input polynomial value, will not be modified)
#PolyVal0 : [Value, Uncertainty]
#PolyVal1 : [Value, Uncertainty]
#PolyVal2 : [Value, Uncertainty]
#ChiSqr : Float (Should be close to number of input data points)
#ReducedChiSqr : Float (Should be close to 1)
#TimeValues : [All input Time values, minus any values removed by the outlier rejection system]
#FluxValues : [All input Flux values, minus any values removed by the outlier rejection system]
#--------------------------


PolynomialOrder = Results["PolynomialOrder"]

#Copy 'Results' time values, these may be diferent than the input values because some points may have been removed by the outlier rejection system
NewXValues = Results["TimeValues"]
#Calculate the models flux values at each of the time values
NewYValues = ConvertParamatersToFluxValues(NewXValues, Results)

#Plot original time and flux values
matplot.plot(DataX, DataY, "-", label="Input Data")
#Plot models time and flux values
matplot.plot(NewXValues, NewYValues, "-", label="Optimized Function")

#Add text box to show model Chi values
BestChiTxt = "Optimized χ2 : " + str(Results["ChiSqr"])
ReducedChiTxt = "Reduced χ2 : " + str(Results["ReducedChiSqr"])

#Congigure text box
ChiAnchoredTextBox = AnchoredText((BestChiTxt + "\n" + ReducedChiTxt), loc=4, pad=0.5)
matplot.setp(ChiAnchoredTextBox.patch, facecolor="Orange", alpha=0.5)
matplot.gca().add_artist(ChiAnchoredTextBox)
matplot.legend(loc=2, borderaxespad=0)

#Show created plots
matplot.show()

#----------------------