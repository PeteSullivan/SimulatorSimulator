import subprocess
import numpy as np
from scipy.io import mmwrite

def find_equilibrium_positions(NIons, secular_frequencies, ShowInfo=False):
    '''
    Finds the equilibrium positions of ions given trap frequencies
    '''
    wx, wy, wz = secular_frequencies * 10**6

    #run mathematica script to find positions and export them to a file
    subprocess.run([
        "wolframscript", "-script", "Z:\\Users\\Pete\\SimulatorSimulator\\Simulator\\MathematicaScripts\\equil_positions.wl",
        str(NIons), str(wx), str(wy), str(wz), str(ShowInfo)
    ], check=True)
    print("found equilibrium positions")
    #read csv file to get positions
    data = np.genfromtxt('Z:\\Users\\Pete\\SimulatorSimulator\\Simulator\\TempFiles\\equil_positions.csv', delimiter=',', skip_header=0)
    
    return data

def current_equilibrium_positions():
    '''
    returns the most recently calculated equilibrium positions
    '''

    return np.genfromtxt('Z:\\Users\\Pete\\SimulatorSimulator\\Simulator\\TempFiles\\equil_positions.csv', delimiter=',', skip_header=0)

def find_secular_frequencies(trapName, NIons, Vrf, V, ImportComsol=False):
    '''
    Finds the eigensystem from positions and trap frequencies
    '''
    V1, V2, V3 = V
    scriptPath = "Z:\\Users\\Pete\\SimulatorSimulator\\Simulator\\MathematicaScripts\\secular_frequencies_" + trapName + ".wl"
    outputPath = "Z:\\Users\\Pete\\SimulatorSimulator\\Simulator\\TempFiles\\secular_frequencies_" + trapName + ".csv"
    subprocess.run([
        "wolframscript", "-script", scriptPath,
        str(NIons), str(Vrf), str(V1), str(V2), str(V3), str(ImportComsol)
    ], check=True)
    print("found secular frequencies")
    #read csv file to get positions
    data = np.genfromtxt(outputPath, delimiter=',', skip_header=0)
    
    return data

def current_secular_frequencies(trapName):
    '''
    returns the most recently calculated equilibrium positions for a specific trap
    '''

    outputPath = "Z:\\Users\\Pete\\SimulatorSimulator\\Simulator\\TempFiles\\secular_frequencies_" + trapName + ".csv"
    return np.genfromtxt(outputPath, delimiter=',', skip_header=0)

def Jij_to_weights(JDes, NormalModeEigVecs):
    '''
    returns the weights needed to calculate arbitrary Jij calculations'''
    mmwrite("TempFiles\\JDes.mtx", JDes) # Export to MTX format
    mmwrite("TempFiles\\normal_modes.mtx", NormalModeEigVecs)


    scriptPath = "Z:\\Users\\Pete\\SimulatorSimulator\\Simulator\\MathematicaScripts\\Jij_weights.wl"
    outputPath = "Z:\\Users\\Pete\\SimulatorSimulator\\Simulator\\TempFiles\\Jij_weights.csv"
    subprocess.run([
        "wolframscript", "-script", scriptPath,
    ], check=True)
    print("found secular frequencies")
    #read csv file to get positions
    data = np.genfromtxt(outputPath, delimiter=',', skip_header=0)
    
    return data
