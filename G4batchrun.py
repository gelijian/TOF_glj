import os
import time
import numpy as np

def createTOFin(En):
    """
    generate the input file(TOF.in) for the G4 simulation 
    """
    ftemplate = open("TOFtemplate.in", "r")
    lines = ftemplate.readlines()
    ftofin = open("TOF.in", "w")   
    energyline = lines[12].split()
    lines[12] = "%s %g %s\n"%(energyline[0], En, energyline[2])
    ftofin.writelines(lines)
    ftemplate.close()
    ftofin.close()

def mergeThreads ():
    """ merge the 4 thread results into one file """
    
    print "merge threads starts"
    data = np.zeros((1, 6))
    for threadID in xrange(0, 4):
        filename = "TOF_nt_TOF_t%d.csv"%(threadID)
        dataThread = np.loadtxt(filename, delimiter = ",", skiprows = 9)
        data = np.vstack((data, dataThread))
    savefmt = ["%.5g", "%.5g", "%.5g", "%d", "%d", "%d"]
    np.savetxt("TOFfile.dat", data[1:, :], fmt = savefmt)
    print "merge threads finished"
# create the directory 

directory = "TOFED_RF_no_digitizer_resolution"
elementdir = "G4csvdata"
if (not os.path.exists(directory)):
    os.mkdir(directory)
    os.mkdir(os.path.join(directory, elementdir))
    print "%s has been created!" % (directory)
    print "%s has been created!" % (os.path.join(directory, elementdir))
    
# set the neutron energy list

Enlist = np.arange(3180, 3520, 20)
#np.savetxt(os.path.join(directory, "Enlist"), Enlist, fmt = "%g")

for i in xrange(len(Enlist)):
    start = time.time()
    En = Enlist[i]    
    createTOFin(En)
    
    print "TOFED G4 Simulation for %g keV neutron starts..." % (En)
    os.system("./TOF TOF.in > TOF.out")
    mergeThreads()
    filename = "En_%dkeV" % (En)
    os.system("mv TOFfile.dat %s" % (os.path.join(directory, elementdir, filename)))
    stop = time.time()
    runtime = stop - start
    print "TOFED G4 Simulation for %g keV neutron finished, run time: %g" % (En, runtime)
    print ""
  
