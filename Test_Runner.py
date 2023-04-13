import Transit_Curve_Fitter

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

print(Transit_Curve_Fitter(DataX, DataY, None, None))