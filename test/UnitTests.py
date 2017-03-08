'''Runs 15 unique test orbits with IMPACT-prop and compares to STK solutions'''

import numpy as np
import csv
import math
import os
import subprocess
import UnitTests


def run_test_cases():

    outdir = "./IP_Test_Cases"                           #Define the output directory
    N = 15                                               #Define number of test cases
    passtest = 0                                         #Define counter for passed tests
    kep, time, satcon, msis, drag = generate_input(N)    #Generate control data for input file
    forces = generate_forces(N)                          #Generate force data for input file
    compile_code()                                       #Compile the propagator

    for i in range(N):  #Loop over test cases

        modify_input(kep, time, satcon, msis, drag, forces[i], i)    #Rewrite input file with test case data
        run_propagator()                                             #Run the propagator with the test case data
        copy_output(i, outdir)                                       #Rename and copy output to a new directory
        MSSE = case_errors()                                         #Define maximum errors for each test case
        passtest += check_results(time, outdir, MSSE[i], i)          #Compare test case results to STK

    print str(passtest)+" of "+str(N)+" tests passed"



def compile_code():
    os.chdir("../src")
    command = "make"
    os.system(command)
    os.chdir("../test")
    print "Compiled propagator"


def generate_input(N):

    #Define orbit cases
    oc = range(N)                 #All orbits
    co = range(2,11)+range(13,N)  #Complex 
    eq = [0, 11]                  #Circular equatorial
    po = [1, 12]                  #Circular polar

    #Wrap orbit cases in a single array
    orbits = np.array([oc, co, eq, po])
                 
    ########### Define test case arrays for orbit elements ##########
    sma = np.zeros(N)     #Semi-major axis (km)
    ecc = np.zeros(N)     #Eccentricity (unitless)
    inc = np.zeros(N)     #Inclination (degrees)
    aop = np.zeros(N)     #Argument of Perigee (degrees)
    raan = np.zeros(N)    #Right Ascension of the Ascending Node (degrees)
    ta = np.zeros(N)      #True Anomaly (degrees)

    #Define semi-major axis for all cases
    for i in orbits[0]:
        sma[i] = 6778.14
        
    #Define complex orbit elements
    for i in orbits[1]:
        ecc[i] = 0.01
        inc[i] = 50.0
        aop[i] = 40.0
        raan[i] = 30.0
        ta[i] = 20.0
        
    #Define circular equatorial elements
    for i in orbits[2]:
        ecc[i] = 1.0e-10

    #Define circular polar elements
    for i in orbits[3]:
        ecc[i] = 1.0e-10
        inc[i] = 90.0

    #Wrap orbit elements in single array
    kep = np.array([sma, ecc, inc, aop, raan, ta])

    ############ Define time variables #############
    utc = "2002 OCT 22 00:00:00"    #Start time (UTC)
    tst = 86400                     #Total simulation time (seconds)
    ts = 10.0                       #Time step (seconds)

    #Wrap time variables in a single array
    time = np.array([utc, tst, ts])

    ########### Define MSIS variables ##############
    f107 = 150.0
    f107A = 150.0
    ap = 0.0
    
    #Wrap MSIS variables in a single array
    msis = np.array([f107, f107A, ap])

    ########### Define Drag variables ##############
    geo = 1
    gsi = 1
    ads = 1
    pitch = 0.0
    yaw = 0.0
    cylinder_radius = 1.0
    cylinder_length = 1.0
    cuboid_height = 1.0
    cuboid_width = 1.0
    surf_mass = 4.479e-26
    sat_mass = 1.0
    proj_area = 2.0e-7
    Cd = 2.2
    cSRP = 1.89

    #Wrap Drag variables in a single array
    drag = np.array([geo, gsi, ads, pitch, yaw, cylinder_radius, cylinder_length, cuboid_height, cuboid_width, \
                         surf_mass, sat_mass, proj_area, Cd, cSRP])
    
    ########### Define satellite control variables ##############
    sca = 24920   #Satellite catalog number
    satens = 1    #Number of ensembles
    BCSD = 0.1    #Ballistic coefficient scale factor standard deviation
    BCM = 1.0     #Ballistic coefficient scale factor mean
    psig = 0.0    #Satellite position 1-sigma standard deviation
    vsig = 0.0    #Satellite velocity 1-sigma standard deviation

    #Wrap satellite control variables in a single array
    satcon = np.array([sca, satens, BCSD, BCM, psig, vsig])

    print "Generated inputs"

    return kep, time, satcon, msis, drag



def generate_forces(N):

    M = 17 #Number of force flags

    #tbg = Two-Body Gravity
    #sph = Spherical Harmonic Gravity Field
    #gdo = Degree and Order of the gravity Field
    #egm = Earth Gravity Model
    #srp = Solar Radiation Pressure
    #drg = Atmospheric Drag
    #dnm = Atmospheric Density Model
    #msi = Dynamic MSIS
    #rsw = Read Space Weather
    #sun = Sun Gravity Perturbations
    #mon = Moon Gravity Perturbations
    #dtl = Density Time Lookup (Interpolation)
    #rko = Runge-Kutta Order
    #ens = Read State Ensembles File
    #kep = Compute Keplerian Elements

    forces = np.array((N, M))
    #Define test case arrays for force flags
                      #tbg sph gdo egm srp drg dnm msi bcs rsw rsm sun mon dtl rko ens kep
    forces = np.array([[1,  0,  5,  1,  0,  0,  1,  1,  0,  0,  0,  0,  0,  0,  4,  0,  0],  #TBG only
                       [1,  0,  5,  1,  0,  0,  1,  1,  0,  0,  0,  0,  0,  0,  4,  0,  0],  #TBG only
                       [1,  0,  5,  1,  0,  0,  1,  1,  0,  0,  0,  0,  0,  0,  4,  0,  0],  #TBG only
                       [1,  0,  5,  1,  0,  0,  1,  1,  0,  0,  0,  1,  0,  0,  4,  0,  0],  #TBG + Sun
                       [1,  0,  5,  1,  0,  0,  1,  1,  0,  0,  0,  0,  1,  0,  4,  0,  0],  #TBG + Moon
                       [1,  0,  5,  1,  0,  0,  1,  1,  0,  0,  0,  1,  1,  0,  4,  0,  0],  #TBG + Sun+Moon
                       [1,  0,  5,  1,  1,  0,  1,  1,  0,  0,  0,  0,  0,  0,  4,  0,  0],  #TBG + SRP
                       [1,  0,  5,  1,  1,  0,  1,  1,  0,  0,  0,  1,  1,  0,  4,  0,  0],  #TBG + SRP+Sun+Moon
                       [1,  1,  5,  1,  0,  0,  1,  1,  0,  0,  0,  0,  0,  0,  4,  0,  0],  #TBG + GDO5
                       [1,  1, 25,  1,  0,  0,  1,  1,  0,  0,  0,  0,  0,  0,  4,  0,  0],  #TBG + GDO25
                       [1,  1, 25,  1,  1,  0,  1,  1,  0,  0,  0,  1,  1,  0,  4,  0,  0],  #TBG + GDO25+SRP+Sun+Moon
                       [1,  0, 25,  1,  0,  1,  1,  1,  0,  0,  0,  0,  0,  0,  4,  0,  0],  #TBG + Drag
                       [1,  0, 25,  1,  0,  1,  1,  1,  0,  0,  0,  0,  0,  0,  4,  0,  0],  #TBG + Drag
                       [1,  0, 25,  1,  0,  1,  1,  1,  0,  0,  0,  0,  0,  0,  4,  0,  0],  #TBG + Drag
                       [1,  1, 25,  1,  1,  1,  1,  1,  0,  0,  0,  1,  1,  0,  4,  0,  0]]) #TBG + Drag+SRP+Sun+Moon+GD25

    print "Generated configuration"

    return forces

    

def modify_input(kep, time, satcon, msis, drag, forces, tcn):

    InputR = open("../data/prop.inp", "r")  #Open input file to read

    #Define search strings
    kepstring = ["Semi-Major Axis (km)", \
                     "Eccentricity", \
                     "Inclination (degrees)", \
                     "Argument of Perigee (degrees)", \
                     "Right-Ascension of the Ascending Node (degrees)", \
                     "True Anomaly (degrees)"]

    timestring = ["Initial UTC Time (Format: 2013 JAN 01 00:00:00)", \
                      "Total Simulation Time (seconds)", \
                      "Timestep (seconds)"]

    satstring = ["Satellite Catalog Number", \
                     "Number of Satellites", \
                     "Satellite Position Distribution 1-sigma Variance", \
                     "Satellite Velocity Distribution 1-sigma Variance"]

    msisstring = ["Daily Observed F10.7 Solar Flux (sfu)", \
                      "81-day Centered F10.7 Solar Flux Average (sfu)", \
                      "Daily Averaged Geomagnetic Index (ap)"]

    dragstring = ["Satellite geometry (See Cd.c for legend)", \
                      "Gas-surface interation model (0=DRIA; 1=CLL)", \
                      "Adsorption model (0=Langmuir; 1=Freundlich)", \
                      "Pitch angle (degrees)", \
                      "Yaw angle (degrees)", \
                      "Cylinder radius (meters)", \
                      "Cuboid or cylinder length (meters)", \
                      "Cuboid height (meters)", \
                      "Cuboid width (meters)", \
                      "Surface material particle mass (kg)", \
                      "Total satellite mass (kg)", \
                      "Projected satellite area (km^2)", \
                      "Constant Satellite Drag Coefficient (if not HF Cd)", \
                      "SRP reflectivity coefficient for satellite surface"]

    configstring = ["Use Two-Body Gravity (Should always be set to 1)", \
                   "Use Spherical Harmonic Gravity Field", \
                   "Gravity Field Degree & Order", \
                   "Earth Gravity Model (0=JGM-3; 1=EGM96)", \
                   "Use Solar Radiation Pressure", \
                   "Use Atmospheric Drag", \
                   "Density Model (0=CIRA72;1=MSIS;2=GITM;3=US1976;4=HASDM)", \
                   "Dynamically Calculate MSIS Atmospheric Properties", \
                   "Read Space Weather Data for MSIS/HASDM", \
                   "Use High-Fidelity Drag Coefficient (Cd) Model", \
                   "Use Response Surface Model (RSM) Cd (GRACE=1; CHAMP=2)", \
                   "Include 3rd Body Gravity Perturbations due to the Sun", \
                   "Include 3rd Body Gravity Perturbations due to the Moon", \
                   "Time interpolation for density (0=nearest; 1=linear)", \
                   "Integration Scheme (4=RK4; 8=RK8; 0=MCPI)", \
                   "Read from State Ensembles File", \
                   "Compute Keplerian Elements each Timestep"]
                      
                      

    strings = [kepstring, timestring, satstring, msisstring, dragstring, configstring]   #Create array for all strings

    text = [] #Create text array to store the input file

    #Read the file
    for line in InputR:
        text.append(line)

    InputR.close() #Close the file

    for iline in range(len(text)):                  #Loop over all lines in the input file
       for i in range(len(strings)):                #Loop over three types of strings
           for j in range(len(strings[i])):         #Loop over the elements of the given string type
               if strings[i][j] in text[iline]:     #Check if the string label is in the line
                   data = text[iline].split("#")
                   if i==0:   #Keplerian elements
                       data[1] = kep[j][tcn]
                   elif i==1: #Time controls
                       data[1] = time[j]
                   elif i==2: #Satellite controls
                       data[1] = satcon[j]
                   elif i==3: #MSIS controls
                       data[1] = msis[j]
                   elif i==4: #Drag controls
                       data[1] = drag[j]
                   else:      #Config controls
                       data[1] = forces[j]
                   text[iline] = data[0]+"# "+str(data[1])+"\n"
      
    InputW = open("../data/prop.inp", "w") #Reopen the file to write

    #Write the new input file
    for iline in range(len(text)):
        InputW.write(text[iline])

    InputW.close() #Close the file


def run_propagator():

    os.chdir("../src")
    command = "./propagator"
    p = subprocess.Popen(command, shell=True, stdout = subprocess.PIPE, stderr=subprocess.STDOUT)
    output = p.communicate()[0]
    os.chdir("../test")


def copy_output(i, outdir):

    if not os.path.exists(outdir):
        os.makedirs(outdir)
    command = "cp ../data/propdata_24920_0.dat "+outdir+"/CProp_TestCase"+str(i)+".dat"
    os.system(command)


def case_errors():
    
    MSSE = np.zeros((15, 6))
    MSSE[0] = np.array([0.05, 0.05, 0.05, 750.0, 750.0, 750.0])     #Circular Equatorial Orbit - TBG only
    MSSE[1] = np.array([0.05, 0.05, 0.05, 750.0, 750.0, 750.0])     #Circular Polar Orbit -      TBG only
    MSSE[2] = np.array([0.05, 0.05, 0.05, 750.0, 750.0, 750.0])     #Complex Orbit -             TBG only
    MSSE[3] = np.array([0.05, 0.05, 0.05, 750.0, 750.0, 750.0])     #Complex Orbit - TBG + Sun
    MSSE[4] = np.array([0.05, 0.05, 0.05, 750.0, 750.0, 750.0])     #Complex Orbit - TBG + Moon
    MSSE[5] = np.array([0.05, 0.05, 0.05, 750.0, 750.0, 750.0])     #Complex Orbit - TBG + Sun + Moon
    MSSE[6] = np.array([5.0e3, 5.0e3, 5.0e3, 3.0e3, 3.0e3, 3.0e3])  #Complex Orbit - TBG + SRP
    MSSE[7] = np.array([5.0e3, 5.0e3, 5.0e3, 3.0e3, 3.0e3, 3.0e3])  #Complex Orbit - TBG + SRP + Sun + Moon
    MSSE[8] = np.array([1.0e4, 1.0e4, 1.0e4, 1.3e4, 1.3e4, 1.3e4])  #Complex Orbit - GDO5
    MSSE[9] = np.array([8.0e3, 8.0e3, 8.0e3, 1.1e4, 1.1e4, 1.1e4])  #Complex Orbit - GDO25
    MSSE[10]= np.array([6.0e3, 6.0e3, 6.0e3, 7.0e3, 7.0e3, 7.0e3])  #Complex Orbit - GDO25 + SRP + Sun + Moon
    MSSE[11]= np.array([1.1e5, 1.1e5, 1.1e5, 1.4e5, 1.4e5, 1.4e5])  #Circular Equatorial Orbit - TBG + Drag
    MSSE[12]= np.array([1.5e6, 1.5e6, 1.5e6, 1.9e6, 1.9e6, 1.9e6])  #Circular Polar Orbit - TBG + Drag
    MSSE[13]= np.array([2.3e5, 2.3e5, 2.3e5, 2.8e5, 2.8e5, 2.8e5])  #Complex Orbit - TBG + Drag
    MSSE[14]= np.array([1.7e5, 1.7e5, 1.7e5, 2.2e5, 2.2e5, 2.2e5])  #Complex Orbit - TBG + Drag + GDO25 + SRP + Sun + Moon
    
    return MSSE


def check_results(time, outdir, MSSE, tcn):

    counter = 0                  #Counter to test each dimension of the satellite state
    success = 0                  #Flag to determine test success

    SumSqError = np.zeros(6)     #Initialize sum square error array

    NTS = int(float(time[1])/float(time[2]) + 1)    #Number of timesteps

    STK = np.zeros((NTS, 7))     #Define STK data array
    CPR = np.zeros((NTS, 7))     #Define the IMPACT-prop data array

    #Read the STK data file
    with open("./STK_Test_Cases/STK_TestCase"+str(tcn)+".csv", "rb") as csvfile:
        STKfile = csv.reader(csvfile, delimiter=',', quotechar='|')
        next(STKfile)
        for i, row in enumerate(STKfile):
            if row != []:
                STK[i] = np.array(map(float,row))
    
    #Read the IMPACT-prop data file
    CPPfile = open(outdir+"/CProp_TestCase"+str(tcn)+".dat", "r")
    next(CPPfile)
    for i, line in enumerate(CPPfile):
        if line != []:
            data = line.strip().split()
            CPR[i] = np.array(map(float, data))
            
            
    CPPfile.close()   #Close the file

    SqError = np.zeros((NTS, 6))  #Define the error array

    for i in range(NTS):    #Loop over number of timesteps
        for j in range(6):  #Loop over state (X, Y, Z, U, V, W)
            if j<3:
               SqError[i][j] = math.pow((STK[i][j+1] - CPR[i][j+1])*1.0e3, 2.0) #Compute square error, convert to meters
            else:
               SqError[i][j] = math.pow((STK[i][j+1] - CPR[i][j+1])*1.0e6, 2.0) #Compute square error, convert to mm/s

    #Compute vector square error as the sum of the square error at each timestep
    for i in range(NTS):
        for j in range(6):
            SumSqError[j] += SqError[i][j]

    #Compare error for each test case to the maximum value allowed (MSSE)
    for i in range(6):
        #print SumSqError[i], MSSE[i]
        if(SumSqError[i] < MSSE[i]):
            counter += 1

    if counter==6:
        success = 1

    if success == 1:
        print "Test Case #"+str(tcn)+" passed"
    else:
        print "Test Case #"+str(tcn)+" failed"
       
    return success

if __name__ == '__main__':
    
    run_test_cases()
