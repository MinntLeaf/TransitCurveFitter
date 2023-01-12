# Transit Curve Generator

## Setup and Running

### Python Version

This was built using Python 3.10, although any version of Python 3 should work.

### Data

Some usable data genrated using 'pi' is included in this repo as `Data.txt`. It has a period of 'pi' and should return a period of aproximately 3.141 when used as input data.

If you want to use different data, replace the data in the file `Data.txt` with your data.

#### Data Formating

Data should be input in the format (X Value),(Y Value),(Error).

Commas "," are used as value separators.

Numbers "1,2,3" and minus signs "-" will be converted into floating point values. All other charachters will be ignored.

A line break signifies a new (X,Y,Error) entry. Empty lines will be ignored.

An error value does not need to be supplied, however, if an Error value is not given for 1 or more data points, it will be ignored for **all** of them.


Examples:

1,2,-3

1 , 2 ,- 3

1     ,2,-  hello  3

(1), 2, -(3)

1 , [2],-3hello

1    ,2,-3

All will be interpreted as 1,2,-3



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

## Install requirements
Ensure your virtual environment is active first!
```
pip3 install -r requirements.txt
```

## Run
```
python Transit_Curve_Generator.py
```

## Updating
Update your requirements.txt file:  
`pip3 freeze > requirements.txt`

Commit changes:  
`git commit -am "Your commit message."`

Push changes to the Internet:  
`git push`
