# Transit Curve Generator

## Setup and Running

### Python Version

This was built using Python 3.10, although any version of Python 3 should work.

### Windows
#### Windows Prerequisites
 - Python
https://www.python.org/downloads/windows/

 - Microsoft C++ Build Tools
Download from: https://visualstudio.microsoft.com/visual-cpp-build-tools/  
When you run this, you **only** need to select and install:
"Desktop development with C++"

### Build and Activate yoru Python Virtual Environment
`cd` into the folder with the code and run:
```
python -m venv venv
.\venv\Scripts\activate
```

### Linux
#### Linux Prerequisites
You **may** need to run these commands on Linux (Debian) (or their equivalent) if you do not commonly run Python apps:

```
sudo apt install python3-venv
sudo apt install python3-dev
```

### Build and Activate yoru Python Virtual Environment
`cd` into the folder with the code and run:
```
python -m venv venv
source venv/bin/activate
```

## Install requirements
Ensure your virtual environment is active first!
```
pip3 install -r requirements.txt
```

## Running

```
#-Import
import Transit_Curve_Fitter

#-Create new fit data object
NewFitData = Transit_Curve_Fitter.FitData()

#-Input data values, uncertainty is not required
NewFitData.Time = Array of time values
NewFitData.Flux = Array of brightness values
NewFitData.Error = Array of uncertainty values

#-Priors, will be used as the starting point of the fit, and will be factored into the internal chisqr evaluation of the fit.
#-First array value is expectec value, the second is the uncertainty of the prior.
#-If you do not want the internal evaluation of the fitting process to consider it's deviation from these priors, input an uncertainty value of None.
NewFitData.t0 = [#,#] time of inferior conjunction
NewFitData.per = [#,#] orbital period
NewFitData.rp =   [#,#]  planet radius (in units of stellar radii)
NewFitData.a =  [#,#]  semi-major axis (in units of stellar radii)
NewFitData.inc =  [#,#]  orbital inclination (in degrees)
NewFitData.ecc = [#,#]  eccentricity
NewFitData.w =  [#,#]  longitude of periastron (in degrees)
NewFitData.u1 = [#,#] limb darkening coefficient 1
NewFitData.u2 = [#,#] limb darkening coefficients 2

#-None will not include a polynomial modifer, 0 will use a 0th order polynomial (equivalent to a scaler), 1 will use a 1st order polynomial, etc.
#-This should be set to the lowest acceptable value, as higher values (particularly past 3) lead to low accuracy results.
NewFitData.PolynomialOrder = None

#Call the fitting function
Results = Transit_Curve_Fitter.TransitCurveFitter.FitTransit(NewFitData)

#-Results will be returned as a dictionary with 'keys' of the names of each fitting value ("t0", "per", "rp", "u1", etc.), and 'values' of arrays of the form [OptimizedValue, Uncertainty]
#-'ChiSqr' and 'ReducedChiSqr' will be included in the dictionary as scalar values (Results["ReducedChiSqr"] = 1.0534)
#-Note that if the polynomial parameter is enabled (Not set to None) the resulting outputs will be optimized to data normalized by the polynomial modifier
#-The polynomial paramters will be retuned as individual scalers with key names of the form "PolyVal1", "PolyVal2", "PolyVal3", etc.
for Key, Value in Results.items():
   print(Key,":",Value)

#-Example returned dictionary:
t0 : [0.989971529178046, 0.00027734524799205324]
per : [3.1415942249113384, 0.0001419041035711419]
rp : [0.11799460726960415, 0.00014759410019655065]
a : [7.99578637925794, 0.014688035264285106]
inc : [83.70354206018568, 0.014687631251790645]
ecc : [0.004102587326899443, 0.014754478787111891]
w : [0.999027233286216, 0.014762842258588102]
u1 : [9.019502834606286e-07, 0.014812200909992225]
u2 : [4.740935422020698e-10, 0.012500869499059857]
PolynomialOrder : 2
PolyVal0 : [9.907240121265204e-05, 3.035337627251473e-06]
PolyVal1 : [0.0001211945794921121, 8.440997129644654e-05]
PolyVal2 : [0.9990565785802801, 0.0004942603059137387]
ChiSqr : 7776.63
ReducedChiSqr : 1.00008

```

###Example

Some usable data generated using pi is included in this repo as `Data.txt`. It has a period of 'pi' and should return a period of approximately 3.141 when used as input data.

If you want to use different data, replace the data in the file `Data.txt` with your data.

## Updating
Update your requirements.txt file:  
`pip3 freeze > requirements.txt`

Commit changes:  
`git commit -am "Your commit message."`

Push changes to the Internet:  
`git push`
