# Transit Curve Generator

## Setup and Running

### Python Version

This was built using Python 3.10, although any version of Python 3 should work.

### Data
A copy of some usable data is included in this repo as `Data.txt`.

If you want to use different data, replace the data in the file `Data.txt` with your data.

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
python3 -m venv venv
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
python3 -m venv venv
.\venv\Scripts\activate
```

## Install requirements
Ensure your virtual environment is active first!
```
pip3 install -r requirements.txt
```

## Run
```
python3 Transit_Curve_Generator.py
```

## Updating
Update your requirements.txt file:  
`pip3 freeze > requirements.txt`

Commit changes:  
`git commit -am "Your commit message."`

Push changes to the Internet:  
`git push`
