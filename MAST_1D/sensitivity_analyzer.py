
#  Extracts sensitivity analysis data

import os

spath = str(os.path.join(os.pardir, "Sensitivity_Analysis"))
cfolder = str(os.path.join(os.pardir, "Output"))
slist = [1.1, 0.9]
cdict = {}
outputaxis = [""]

cvars = os.listdir(cfolder)
for c in cvars:
    if 'Out' in c:
        file = open(cfolder + '//' + c, 'r')
        file = file.readlines()
        str(file)
        node = file[3].split()
        svalue = node[-1]
        cdict[str(c)]= svalue
        cvariable = c[4:]
        outputaxis.append(cvariable)

for s in slist:
    outputfile = os.path.join(spath, str(s) + "_sensitivityoutput.txt")
    outputfile = open(outputfile, 'w')
    outputstrip = str(outputaxis).strip('[')
    outputstrip = str(outputstrip).strip(']')
    outputfile.write(str(outputstrip) + '\n')
    nvars = os.listdir(spath)
    for n in nvars:
        if str(s) in n and "txt" not in n:
            perclist = []
            folder = str(os.path.join(spath, n))
            ovars = os.listdir(folder)
            name = str(n)
            name = name[7:-4]
            perclist.append(name)
            for o in ovars:
                if 'Out' in o and o in cvars:
                    percdiff = 'nan'
                    file = open(folder + '//' + o, 'r')
                    file = file.readlines()
                    str(file)
                    node = file[3].split()
                    svalue = float(node[-1])
                    cvalue =float(cdict[str(o)])
                    if cvalue != 0:
                        percdiff = (svalue-cvalue)/cvalue*100
                    perclist.append(percdiff)
            percwrite = str(perclist).strip('[')
            percwrite = str(percwrite).strip(']')
            outputfile.write(str(percwrite) + '\n')
    outputfile.close()




