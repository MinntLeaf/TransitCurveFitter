import time
import Transit_Curve_Fitter
import numpy as np

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


NewFitData = Transit_Curve_Fitter.FitData()

NewFitData.Time = DataX
NewFitData.Flux = DataY
NewFitData.Error = None

NewFitData.t0 = [0.96, 999]  #t0
NewFitData.per = [3.14159265, 999]  #per
NewFitData.rp =   [0.118, 999]  #rp
NewFitData.a =  [8.0, 999]  #a
NewFitData.inc =  [83.7, 999]  #inc
NewFitData.ecc = [0.0, 999]  #ecc
NewFitData.w =  [1.0, 999]  #w
NewFitData.u1 = [-0.9, 999]  #u1
NewFitData.u2 =  [-0.9, 999]  #u2
NewFitData.PolynomialOrder = 0

StartTime = time.time()

Iterations = 1

while Iterations > 0:
    print(Transit_Curve_Fitter.TransitCurveFitter.FitTransit(NewFitData))
    Iterations-=1

#print("Average Time : " + (time.strftime-StartTime)/float(Iterations))