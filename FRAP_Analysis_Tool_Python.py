# -*- coding: utf-8 -*-
"""
Created on Thu Nov  8 09:2:15 2022

_______________________________________________________________________________

These program was made for the data analysis of single measurements 
    of FRAP with THP1 (beta G-actin citrin) cells
    
It is the 1st part of 3 single programs:
    - 1_Fitting-Curves-Python.py
    - 2_Analysis_Combine_Resultes.py
    - 3_Mean_Result_Plot.py 
_______________________________________________________________________________

author: Gesa Hölzer, Anna Kemeter
        IAOB, Prof. Eggeling Group
        Friedrich-Schiller-Universität Jena
                  ____________________
    
These program is for the analysis of the recovery-curves,
in order to do that you import csv-data, fit the curve, decide the part of data you want to include and what you want to save.
It automatically makes a folder with your results in it.
For using program 2 you need to save the analysed data.
        
                  
Define in this program :

     - your measurement folder (file_path) & the file you want to analyse (file_name)
     - If you want to exlude data points of the end of the measurement (m)
     - What you want to plot/ save
     - If you take a bleach-control into account (and from to where to fit the bleach function: x1, x2)
                  ____________________
        
FRAP measurement standards in order to use this program as it is:
    
    - Image min. 60 s recovery (~200 cycles)
    - 3 or 10 cycles before bleaching
    - R1 bleach acquisition, R2 bleachcontrol acquisition, R3 background
                  ____________________

"""

# -------------------------------------
# Import Packages
import math
import csv
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy.optimize as spo
from scipy.optimize import curve_fit
import os.path
import glob
import datetime
import re

# -------------------------------------
#Inititialise variables for later
TwoExp, SaveData, SaveTwoExp, CombinedPlot, SaveCombindedPlot, OneExp, SaveOneExp, ExpA, SaveExpA, Bleachcontrol, IncludeBleachcontrol, CopyMicroscopeData, Background, SaveBleachcontrol, SaveTwoExptofile, SaveOneExptofile, SaveExpAtofile, ComparisonRegion  = False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False

TimeAxisIn_ms = False # If the acquisition interval was measured in ms or min it will be set to True later
TimeAxisIn_min = False

x, f = 0, 1

now = datetime.datetime.now()
now = now.strftime("%Y-%m-%d_%H-%M-%S")

"""___________________________________________________________________________________________________________________________
__________________________________THIS BLOCK NEEDS TO BE CUSTOMISED!__________________________________________________________
______________________________________________________________________________________________________________________________"""

# File path the csv file with all the fit data is being saved to (needs to be customised)
csv_file_path = '/Users/AnnaKemeter/Desktop/BA/Tabelle_der_Daten.csv'

# File path to the overall folder with the FRAP data (needs to be customised)
Path = '/Users/AnnaKemeter/Desktop/BA/DataAnna/UsedData/' 

# Specific folder that is being looked in for csv files (needs to be customised) ----- if there aren't multiple data folders, then the last foldername can be used here and the path customized above needs to lead to this folder
Loop = ['2023-03-17, 50popc50chol', '2023-03-29, 50popc50chol', '2023-03-29, 100popc', '2023-04-04, 50popc50chol', '2023-04-14, 50popc50chol', '2023-04-20, 50popc50chol','Ziliangs Data', '2023-05-11, 50popc50chol', '2023-05-11, 100popc, 0.1mol%', '2023-05-16, 100popc, 0.1mol%', '2023-05-24, 100popc, 0.1mol%', '2023-06-01, 100popc, 0.1mol%', '2023-06-08, 100popc, 0.1mol%', '2023-06-19, 100popc, 0.1mol%', '2023-06-23, 100popc, 10nM Atto488', '2023-07-13, 100popc, 10nM Atto488','2023-09-20, 100popc, 0.1mol%']
#Loop = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]

# --------- specifiy file -------------------------------------------------------------------------
x = 11   # Number of the data folder in the loop starting with 0

f = 45  # Choose number of file in folder x (f_min = 1) 

m = 1  # amount of data points (m > 0) excluded from fitting in the end of the curve 

w = 2   # customise x_axis length in Plots, it will be the seconds you measured + w sec. (standard: 2)

""" --------- Uncomment which fits and plots you want to show and save (Ctr+1 or remove "# ") --------------"""

#TwoExp = True       #required to start adjusting the fit
#SaveTwoExp = True        # save the plot picture
#SaveTwoExptofile = True  # save the generated data from the fit to a csv file

OneExp = True             # Show plot of exp. fit 1st Order
#SaveOneExp = True        # save the plot picture
#SaveOneExptofile = True  # save the generated data from the fit to a csv file

#ExpA = True               # Show plot of exp. fit to power of alpha
#SaveExpA = True           # save the plot picture
#SaveExpAtofile = True     # save the generated data from the fit to a csv file

#CombinedPlot = True       # Show plot of all three fits in one figure
#SaveCombindedPlot = True  # save the combined plot picture

""" Bleachcontrol """
#Bleachcontrol = True           #  ROI 2 needs to be non-bleached region -> creates plot of the ROI 1 and ROI 2 data and the possible data curve with a bleaching correction done
#SaveBleachcontrol = True       #  save the plot picture of bleachcontrol
#IncludeBleachcontrol = True    #  bleaching correction of original data as shown in Bleachcontrol is done -> new data curve can be fitted now as before

x1 = 1     # amount of data points to exclude in the beginning; x1 > 0
x2 = m     # amount of data points to exclude in the end; standard: m

""" Subtract Background """
#Background = True     # substracts the data of ROI 3 from the original data

""" Comparison experimental region (ROI 4) """
# ComparisonRegion = True   # option to analyze ROI 4 with a simultaneous bleaching experiment to the ROI 1 for comparison 
                       
# normally the first region is being analyzed 
# either ROI 1 or ROI 4 can be fitted by the different fitting options above and only one at a time


"""________________________________________________________________________________________________________________________
____________________________DONE WITH CUSTOMISATION________________________________________________________________________
___________________________________________________________________________________________________________________________"""

# --- Initialising saving folders ---

ResultFolder = Path + Loop[x]  + '/Result Data/'  # Where to save the result parameters
PlotFolder = Path + Loop[x] + '/Result Plots/'  # Where to save the plots

""" ___________________________________________________________________________
Include the by 'f' determined file 
"""

# --- Creating a list with the names of all csv data in the folder ---
file_paths = glob.glob(Path + Loop[x] + "/" + "*.csv")
filenames = [os.path.basename(file_path) for file_path in file_paths]
#filenames.sort()  #sort my filenames to get the original order
filenames.sort(key=lambda x: [int(d) if d.isdigit() else d for d in re.split(r'(\d+)', x)])

f = f-1     # Start counting from 1 (instead 0)

#file = file_paths[f]
file = Path + Loop[x] + "/" + filenames[f]
file_name =filenames[f]

z = 0       # Initialize counter
for z in range(len(filenames)):
    print_file = filenames[z]
    #print_file = print_file[(len(Path)+len(Folder)+len(Loop[x])+1):]
    print(print_file)
   


    
# Print information for user  
print('')
print('Amount of files in folder:', len(filenames))
print('')
print('use for f an int between 1 and ', len(filenames))
print('')
print('Using file: ', file_name)
print('')



"""____________________________________________________________________________
Import relevant part of the LSM980 csv-table as an array
"""

#calculate the radius of the ROI used for bleaching
if ComparisonRegion == False:
    A_column = pd.read_csv(file, delimiter=",", usecols=[2])
    Area = A_column.to_numpy()
    A = Area[3] #Area of the ROI
    print(A)
    R = math.sqrt(A/math.pi) #Radius of the ROI
    print(R)
else:
    A_column = pd.read_csv(file, delimiter=",", usecols=[15])
    Area = A_column.to_numpy()
    A = Area[3] #Area of the ROI
    print(A)
    R = math.sqrt(A/math.pi) #Radius of the ROI
    print(R)

if ComparisonRegion == False:
    if Bleachcontrol == True or IncludeBleachcontrol == True:
        df = pd.read_csv(file, delimiter=",", usecols=[0, 5, 8, 13], skiprows=1)
    else:
        df = pd.read_csv(file, delimiter=",", usecols=[0, 5, 13], skiprows=1)

else: 
    if Bleachcontrol == True or IncludeBleachcontrol == True:
        df = pd.read_csv(file, delimiter=",", usecols=[0, 17, 8, 13], skiprows=1)
        desired_columns = [0, 17, 8, 13]
        df = df[desired_columns]
    else:
        df = pd.read_csv(file, delimiter=",", usecols=[0, 17, 13], skiprows=1)
        # Define the desired column order
        desired_columns_order = [0, 2, 1]  # Specify the order you want
        # Reorder the columns based on the desired order
        df = df.iloc[:, desired_columns_order]


data = df.to_numpy() # Convert the entire DataFrame to numpy array


if df.columns[0] == "ms":
    TimeAxisIn_ms = True
    
if df.columns[0] == "min":
    TimeAxisIn_min = True

if TimeAxisIn_min == True:  # If the aquisition intervall was measured in ms instead of s
    data[:, 0] = data[:, 0]*60
    
if TimeAxisIn_ms == True:  # If the aquisition intervall was measured in ms instead of s
    data[:, 0] = data[:, 0]/1000

data_length = len(data[:, 0]) # Amount of data points 
#print(data_length)

df_2 = pd.read_csv(file, delimiter=",", usecols=[1])  #import the markers row in the csv.file
first_non_empty_row_index = df_2.index[df_2.notnull().any(axis=1)][0]

b = first_non_empty_row_index-1   #find the amount of times that intensity was measured before the bleaching

"""____________________________________________________________________________
Check if .csv file exist otherwise create a .csv file 
"""


# Check if the CSV file exists
if os.path.exists(csv_file_path):
        print(csv_file_path)
else:
    # If the file doesn't exist, create a new CSV file and write data in there
    with open(csv_file_path, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter='\t')
        writer.writerow(['Fit used', 'Photocorrected', 'Backgroundcorrected', 'Second region for comparison','alpha', 'IE', 'tau_half', 'Error_tau_half', 'tau_2_half', 'Error_tau_half', 'w', 'D', 'Error_D','filename','lipid solution used'])
        print(csv_file_path)
        

"""____________________________________________________________________________
Implementation of background (from an extra measurement in a ROI outside of my membrane) (Anna)
"""

if Background == True:
    
    def linear_background(t, y, b):       # Define fit for the background noise
        return y + b * t  
    
    popt_lin_back, pcov = spo.curve_fit(linear_background, data[x1:-x2, 0], data[x1:-x2, 3], p0=[data[0, 3], 1], method='lm')  # actual Fit, lm=levenberg-marquardt
    popt_lin_back
    #print(popt_lin_back)
    back_ground = linear_background(data[:,0], *popt_lin_back)
    
    data[:,1] = np.subtract(data[:,1], back_ground) # Adapt data of bleached ROI
    data[:,2] = np.subtract(data[:,2], back_ground) # Adapt data of non-bleached ROI
else: 
    data = data #no adaptation
    
    
"""____________________________________________________________________________
Implementation of bleach function (from an extra measurement) 
"""

if Bleachcontrol == True:

    def linear_bleach(t, y, b):
        return y + b * t            # Define fit function 
    
    popt_lin_bl, pcov = spo.curve_fit(linear_bleach, data[x1:-x2, 0], data[x1:-x2, 2], p0=[data[0, 2], 1], method='lm')   # Fit, lm=levenberg-marquardt
    popt_lin_bl
    print(popt_lin_bl)

    bleach_relation = linear_bleach(data[:, 0], *popt_lin_bl) / data[x1, 2] # Calculate photofade factor for each time point (photofading factor)
    bleach_correction = 1/bleach_relation 

    ndata = np.vstack((data[:, 0], np.multiply(data[:, 1], bleach_correction[:]))).transpose() # Adapt data
    
    


    
    """ --- Plot recovery data, bleach-control data and bleachimplement in a graph ---"""
    
    plt.figure(figsize=(10, 5))
    
    plt.scatter(data[:, 0], data[:, 1], marker='.', # Recovery Data 
                label='Measurement', color='blue') 
    plt.scatter(ndata[:, 0], ndata[:, 1], marker='.',
                label='Measurement', color='black') # Bleachimplement
    plt.scatter(data[:, 0], data[:, 2], marker='.',
                label='Measurement', color='orange') # Bleach-control
    plt.plot(data[x1:-x2, 0], linear_bleach(data[x1:-x2, 0], *popt_lin_bl)) # Linear fit of bleachcontrol
    
    plt.xlabel('Time in s')
    plt.grid()
    
    if not os.path.isdir(PlotFolder): # If required, create result plot folder
        os.makedirs(PlotFolder)
    
    if SaveBleachcontrol == True:    
        plt.savefig(PlotFolder+file_name[:-4]+'-Bleachimplement_'+ now + '_' + str(x1) + '_' +str(x2) + '.png', format="png") # Save plot
    plt.show()

   

    if IncludeBleachcontrol == False: # No adaptation
        ndata = data[:, :2]
else:
    ndata = data # No adaptation


"""____________________________________________________________________________
Translate the start of the recovery curve (IA) to the (0,0)-coordinate 
"""
   
ndata = ndata [:,:] - ndata[b, :]

"""____________________________________________________________________________
Normalisation (Min-Max)
"""
Norm = 0
i = 0
for i in range(b):
    Norm = Norm + ndata[i, 1]  # Summing up the values
    i = i+1

norm = Norm/b  # Calculating the average by dividing by the number of data points


ndata = np.hstack((ndata[:, 0], ndata[:, 1]/norm))
ndata = ndata.reshape(2, data_length).transpose() 


"""____________________________________________________________________________
Additional parameters
"""

m = -m # "-" makes the counter start backwards from the last measured value

# Recovery part of data
t = ndata[b:m, 0]  # Time
F = ndata[b:m, 1]  # R1: bleached part
# print(F)


#%%
# ------------------- Fit with 2-component  exponential function-----------------------------------


def exp_decay_two(t, y_2, A_1, t_1, A_2, t_2):
    return y_2 + A_1 * (1 - np.exp(-t / t_1)) + A_2 * (1 - np.exp(-t / t_2))


# ------------ Set Starting Parameters ------------------------
#          = [y_2, A_1,t_1,A_2,t_2]
#guess_two = [0, 5, 2, 9, 26]    # try something out 
guess_two = [0, 23, 4, 9, 26]    # works most of the time


# ------------------- Fitting ---------------------------
if TwoExp == True:
    popt_two, pcov = spo.curve_fit(exp_decay_two, t, F, p0=guess_two, method='lm')
    popt_two
    d_pop_two = np.sqrt(np.diag(pcov))
    
    y_2, A_1, t_1, A_2, t_2 = popt_two
    d_y_2, d_A_1, d_t_1, d_A_2, d_t_2 = d_pop_two



    print()
    print('For the fit function y0 + A1 *(1 - exp(-t/τ1)) + A2 *(1 - exp(-t/τ2)):')
    print()
    print('y0=', y_2, ' +- ', d_y_2)
    print('A1=', A_1, ' +- ', d_A_1)
    print('τ1=', t_1, ' +- ', d_t_1)
    print('A2=', A_2, ' +- ', d_A_2)
    print('τ2=', t_2, ' +- ', d_t_2)
    print()

# ------------ Calculate Residuals -----------------
    
    Residual = abs(F - exp_decay_two(t, *popt_two)) # Difference between each measurement value and fit (in absolute values)
    Res_Mean = np.mean(Residual)

"""____________________________________________________________________________
Calculate additional parameters:
"""
if TwoExp == True:
    IB = 1
    IA = ndata[b, 1]  # = 0 = y_0
    # IB = 1 included with normalisation
        
    IE = A_1 + A_2

    k_1 = 1/t_1
    k_2 = 1/t_2
    t_1_half = - math.log(0.5)*t_1
    d_t_1_half = - math.log(0.5)*d_t_1
    t_2_half = - math.log(0.5)*t_2          # Diffusion half time
    d_t_2_half = - math.log(0.5)*d_t_2
    D_TwoExp = R**2/(4*t_1_half)            # Diffusion coefficient
    d_D_TwoExp = (4*R**2)/(4*t_1_half)**2*d_t_1_half   # error of D because of the fitting error for tau_1 (tau_2 is not included!!)
        
    # I_mobile = (IE - IA)/(IB - IA) = IE
    I_immobile = 1 - IE
        
    # Fractions 
    F_1 = A_1/IE 
    F_2 = A_2/IE 
    F_immobile = I_immobile 
    F_mobile = IE 


    print('IE=', IE)
    print('I_immobile=', I_immobile)
    print('k1=', k_1)
    print('τ1_half=', t_1_half)
    print('k2=', k_2)
    print('τ2_half=', t_2_half)
    print()
    
"""____________________________________________________________________________
Plot data & fit funktion
"""
# 2
if TwoExp == True:
    
    # Plot (ax) and Residuals (Res)
    fig2, (ax, Res) = plt.subplots(2, 1, gridspec_kw={'top': 0.9,
        'height_ratios': [4, 1]}, figsize=(10, 6.5), constrained_layout=True) 
    fig2.suptitle('Fit with Two Exponential Summands', fontsize=20)
    if ComparisonRegion == True:
        fig2.text(0.5, 0.91, file_name + '  Second Region for Comparison', fontsize=10, ha='center')
    else: 
        fig2.text(0.5, 0.91, file_name, fontsize=10, ha='center')
    
    if IncludeBleachcontrol == True and Background == False:
        ax.scatter(ndata[:, 0], ndata[:, 1], marker='.',
                   label='Measurement (photofade corrected)', color='k')
    if Background == True and IncludeBleachcontrol == False:
        ax.scatter(ndata[:, 0], ndata[:, 1], marker='.',
                      label='Measurement (background subtracted)', color='k')
    if IncludeBleachcontrol == True and Background == True:
        ax.scatter(ndata[:, 0], ndata[:, 1], marker='.',
                      label='Measurement (photofade and background corr)', color='k')
    if IncludeBleachcontrol == False and Background == False:
        ax.scatter(ndata[:, 0], ndata[:, 1], marker='.',
                   label='Measurement', color='k')

    ax.plot(t, exp_decay_two(t, *popt_two), 'g-',
            label='Fit: y0 + A1 *(1 - exp(-t/τ1)) + A2 *(1 - exp(-t/τ2)) \n y0=%5.3f, A1=%5.3f, τ1=%5.3f,\n A2=%5.3f, τ2=%5.3f' % tuple(popt_two))
    ax.grid()
    ax.set_ylabel('Intensity (normalized)', fontsize=16)
    #ax.set_xlabel('Time in s', fontsize=14)
    ax.set_xlim(-5, ndata[-1, 0]+w)
    ax.legend(loc='best', fontsize=14)
    

    Res.plot(t, Residual)
    Res.set_ylabel('Residuals', fontsize=16)
    Res.set_xlabel('Time in s', fontsize=16)
    Res.axhline(y=Res_Mean, xmin=0, xmax=60, color='r', linewidth=0.5,
                label='Mean difference: %5.3f' % float(Res_Mean))
    Res.legend(loc='best', fontsize=14)
    Res.set_xlim(-5, ndata[-1, 0]+w)
    Res.grid()

    fig2.align_ylabels((ax, Res)[:])
    fig2.tight_layout()
    fig2.show()

    if SaveTwoExp == True:
        
        if not os.path.isdir(PlotFolder):
            os.makedirs(PlotFolder)
        
        if ComparisonRegion == True:
            if IncludeBleachcontrol ==  True:
                fig2.savefig(PlotFolder+file_name[:-4]+ '--ROI 4'+'---Fit-2nd-Order-Plot_+PhotofadeCorr'+ now + '_' + str(x1) + '_' +str(x2) +'.pdf', format="pdf")   
            else: # different file name
                fig2.savefig(PlotFolder+file_name[:-4]+ '--ROI 4'+'---Fit-2nd-Order-Plot_noCorr'+ now + '_' + str(x1) + '_' +str(x2) +'.pdf', format="pdf")
        else:
            if IncludeBleachcontrol ==  True:
                fig2.savefig(PlotFolder+file_name[:-4]+'---Fit-2nd-Order-Plot_+PhotofadeCorr'+ now + '_' + str(x1) + '_' +str(x2) +'.pdf', format="pdf")   
            else: # different file name
                fig2.savefig(PlotFolder+file_name[:-4]+'---Fit-2nd-Order-Plot_noCorr'+ now + '_' + str(x1) + '_' +str(x2) +'.pdf', format="pdf")
        
    if SaveTwoExptofile == True:
        #include the generated data
        if ComparisonRegion == True: 
            with open(csv_file_path, 'a', newline='') as csvfile:
                    csvwriter = csv.writer(csvfile, delimiter='\t')
                    csvwriter.writerow(['TwoExp', IncludeBleachcontrol, Background, ComparisonRegion, '-', IE, t_1_half, d_t_1 ,t_2_half, d_t_2 ,R , D_TwoExp, d_D_TwoExp , file_name + '  Comparison region', Loop[x]])
                    #csvfile.write('\n')  # Add a newline
        else: 
            with open(csv_file_path, 'a', newline='') as csvfile:
                    csvwriter = csv.writer(csvfile, delimiter='\t')
                    csvwriter.writerow(['TwoExp', IncludeBleachcontrol, Background, ComparisonRegion, '-', IE, t_1_half, d_t_1, t_2_half, d_t_2 ,R , D_TwoExp, d_D_TwoExp , file_name, Loop[x]])
                    #csvfile.write('\n')  # Add a newline

#%%
# ------------------- Fit with 1 exponential function-----------------------------------

# Define Fit Function

def exp_decay_one(t, y_1, A_0, t_0):
    return y_1 + A_0 * (1 - np.exp(-t / t_0))

# Set Starting Parameters
#          = [ B_1, A_0, t_0]
#guess_one = [0, 28, 10.5]   # try something out here 
guess_one = [0, 30, 4]       # works most of the time


# Fit
if OneExp == True:
    popt_one, pcov = curve_fit(exp_decay_one, t, F, p0=guess_one, method='lm')
    #popt_one, pcov = curve_fit(exp_decay_one, t, F, p0=guess_one, method='dogbox') 
    popt_one
    d_popt_one = np.sqrt(np.diag(pcov)) # Standard deviation
    
    y_1, A_0, t_0 = popt_one
    d_y_1, d_A_0, d_t_0 = d_popt_one


    print('For the fit function y0 +  A * (1 - exp(-t/τ)) :')
    print('y0=', y_1, ' +- ', d_y_1)
    print('A=', A_0, ' +- ', d_A_0)
    print('τ=', t_0, ' +- ', d_t_0)
    print()

# Residuals
    Residual_one = abs(F - exp_decay_one(t, *popt_one))
    Res_Mean_one = np.mean(Residual_one)

"""____________________________________________________________________________
Calculate additional parameters:
"""

if OneExp == True:
    IB_1 = 1
    IA_1 = ndata[b, 1]  # = 0 = y_0
    # IB = 1 included with normalisation
    
    IE_1 = y_1 + A_0
    
    k_0 = 1/t_0
    t_0_half = - math.log(0.5)*t_0       # Diffusion half time
    d_t_0_half = - math.log(0.5)*d_t_0
    D_OneExp = R**2/(4*t_0_half)         # Diffusion coefficient
    d_D_OneExp = (4*R**2)/(4*t_0_half)**2*d_t_0_half  # error of D because of the fitting error for tau
    
    # I_mobile = (IE - IA)/(IB - IA) = IE
    I_immobile_0 = 1 - IE_1
    
    # Fractions 
    F_0 = (A_0)/(IE_1)
    F_immobile_0 = I_immobile_0 
    F_mobile_0 = IE_1 


    print('IE=', IE_1)
    print('I_immobile=', I_immobile_0)
    print('k=', k_0)
    print('τ_half=', t_0_half)
    print()

# Create Plot
if OneExp == True:

    fig1, (ax, Res) = plt.subplots(2, 1, gridspec_kw={'top': 0.9,
        'height_ratios': [4, 1]}, figsize=(10, 6.5), constrained_layout=True)
    fig1.suptitle('Fit with One Exponential Summand', fontsize=20)
    if ComparisonRegion == True:
        fig1.text(0.5, 0.91, file_name + '  Second Region for Comparison', fontsize=10, ha='center')
    else: 
        fig1.text(0.5, 0.91, file_name, fontsize=10, ha='center')
    
    if IncludeBleachcontrol == True and Background == False:
        ax.scatter(ndata[:, 0], ndata[:, 1], marker='.', label='Measurement (photofade corrected)', color='k')
    if Background == True and IncludeBleachcontrol == False:
        ax.scatter(ndata[:, 0], ndata[:, 1], marker='.', label='Measurement (background subtracted)', color='k')
    if Background == True and IncludeBleachcontrol == True:
        ax.scatter(ndata[:, 0], ndata[:, 1], marker='.', label='Measurement (photofade and background corr)', color='k')       
    if IncludeBleachcontrol == False and Background == False:
        ax.scatter(ndata[:, 0], ndata[:, 1], marker='.', label='Measurement', color='k')
        
    ax.plot(t, exp_decay_one(t, *popt_one), 'b-', label='Fit: y0 + A * (1 - exp(-t/τ))\n y0=%5.3f, A=%5.3f,τ=%5.3f' % tuple(popt_one))
    ax.grid()
    ax.set_ylabel('Intensity (normalized)', fontsize=16)
    # ax.set_xlabel('Time in s', fontsize=14)
    ax.set_xlim(-5, ndata[-1, 0]+w)
    ax.legend(loc='best', fontsize=14)

    Res.plot(t, Residual_one)
    Res.set_ylabel('Residuals', fontsize=16)
    Res.set_xlabel('Time in s', fontsize=16)
    Res.axhline(y=Res_Mean_one, xmin=0, xmax=60, color='r', linewidth=0.5,
                label='Mean difference: %5.3f' % float(Res_Mean_one))
    Res.legend(loc='best', fontsize=14)
    Res.set_xlim(-5, ndata[-1, 0]+w)
    Res.grid()

    fig1.align_ylabels((ax, Res)[:])
    fig1.tight_layout()
    fig1.show()

    if SaveOneExp == True:
        if not os.path.isdir(PlotFolder):
            os.makedirs(PlotFolder)
        if ComparisonRegion == True:    
            if IncludeBleachcontrol ==  True:
                fig1.savefig(PlotFolder+file_name[:-4]+ '--ROI 4'+'---Fit-1st-Order-Plot_+PhotofadeCorr'+ now + '_' + str(x1) + '_' +str(x2) +'.pdf', format="pdf")
            else:
                fig1.savefig(PlotFolder+file_name[:-4]+ '--ROI 4'+'---Fit-1st-Order-Plot_noCorr'+ now + '_' + str(x1) + '_' +str(x2) +'.pdf', format="pdf")
        else: 
            if IncludeBleachcontrol ==  True:
                fig1.savefig(PlotFolder+file_name[:-4]+'---Fit-1st-Order-Plot_+PhotofadeCorr'+ now + '_' + str(x1) + '_' +str(x2) +'.pdf', format="pdf")
            else:
                fig1.savefig(PlotFolder+file_name[:-4]+'---Fit-1st-Order-Plot_noCorr'+ now + '_' + str(x1) + '_' +str(x2) +'.pdf', format="pdf")
        
    if SaveOneExptofile == True:
        #include the generated data
        if ComparisonRegion == True:
            with open(csv_file_path, 'a', newline='') as csvfile:
                    csvwriter = csv.writer(csvfile, delimiter='\t')
                    csvwriter.writerow(['OneExp', IncludeBleachcontrol, Background, ComparisonRegion, '-', IE_1, t_0_half, d_t_0_half,'-', '-',R , D_OneExp, d_D_OneExp , file_name + '  Comparison region', Loop[x]])
                    #csvfile.write('\n')  # Add a newline
        else:
            with open(csv_file_path, 'a', newline='') as csvfile:
                    csvwriter = csv.writer(csvfile, delimiter='\t')
                    csvwriter.writerow(['OneExp', IncludeBleachcontrol, Background, ComparisonRegion, '-', IE_1, t_0_half, d_t_0_half,'-', '-',R , D_OneExp, d_D_OneExp , file_name , Loop[x]])
                    #csvfile.write('\n')  # Add a newline

#%%
# ------------------- Fit with one exponatial summand to the power of alpha -----------------------------------

# Define Fit Function

def exp_decay_a(t, y_a, A_a, t_a, alpha):
    return y_a + A_a * (1 - np.exp(-(t/t_a)**alpha)) # without setting it here: alpha < 1

# Set Starting Parameters
#        = [y_a, A_a, t_a, alpha]
#guess_a = [0, 20, 10, 1] #try something out 
guess_a = [0, 30, 6, 0.7] #works most of the time


# Fit
if ExpA == True:
    #popt_a, pcov = curve_fit(exp_decay_a, t, F, p0=guess_a, method='dogbox') #try something out (Anna)
    popt_a, pcov = curve_fit(exp_decay_a, t, F, p0=guess_a, bounds=(-5, 500))
    popt_a
    d_popt_a = np.sqrt(np.diag(pcov))
    y_a, A_a, t_a, alpha = popt_a
    d_y_a, d_A_a, d_t_a, d_alpha = d_popt_a


    print('For the fit function y0 + A * (1 - exp((-t/τ)**a) ) :')
    print('y0=', y_a, ' +- ', d_y_a)
    print('A=', A_a, ' +- ', d_A_a)
    print('τ=', t_a, ' +- ', d_t_a)
    print('alpha =', alpha, ' +- ', d_alpha)
    print()

# Residuals
    Residual_a = abs(F - exp_decay_a(t, *popt_a))
    Res_Mean_a = np.mean(Residual_a)

"""____________________________________________________________________________
Calculate additional parameters:
"""
if ExpA == True:
    IB_a = 1
    IA_a = ndata[b, 1]  # = 0 = y_0
    # IB = 1 included with normalisation
    
    IE_a = y_a + A_a
    
    k_a = 1/t_a
    t_a_half = - math.log(0.5)*t_a      # Diffusion half time
    d_t_a_half = - math.log(0.5)*d_t_a
    D_ExpA = R**2/(4*t_a_half)          # Diffusion coefficient
    d_D_ExpA = (4*R**2)/(4*t_a_half)**2*d_t_a_half  # error of D because of the fitting error for tau
        
    
    # I_mobile = (IE - IA)/(IB - IA) = IE
    I_immobile_a = 1 - IE_a
    
    # Fractions 
    F_a = (A_a)/(IE_a)
    F_immobile_a = I_immobile_a 
    F_mobile_a = IE_a


    print('IE=', IE_a)
    print('I_immobile=', I_immobile_a)
    print('k=', k_a)
    print('τ_half=', t_a_half)
    print()

if ExpA == True:
    figA, (ax, Res) = plt.subplots(2, 1, gridspec_kw={'top': 0.9,
        'height_ratios': [4, 1]}, figsize=(10, 6.5), constrained_layout=True)
    figA.subplots_adjust(top=0.88) # Increase top margin
    
    figA.suptitle('Fit with One Exponential Summand to Power of α', fontsize=20)
    if ComparisonRegion == True:
        figA.text(0.5, 0.91, file_name + '  Second Region for Comparison', fontsize=10, ha='center')
    else: 
        figA.text(0.5, 0.91, file_name, fontsize=10, ha='center')

    if IncludeBleachcontrol == True and Background == False:
        ax.scatter(ndata[:, 0], ndata[:, 1], marker='.',
                   label='Measurement (photofade corrected)', color='k')
    if Background == True and IncludeBleachcontrol == False:
        ax.scatter(ndata[:, 0], ndata[:, 1], marker='.',
                   label='Measurement (background subtracted)', color='k')
    if Background == True and IncludeBleachcontrol == True:
        ax.scatter(ndata[:, 0], ndata[:, 1], marker='.',
                   label='Measurement (background and photofade corr)', color='k')
    if IncludeBleachcontrol == False and Background == False:
        ax.scatter(ndata[:, 0], ndata[:, 1], marker='.',
                   label='Measurement', color='k')

    ax.plot(t, exp_decay_a(t, *popt_a), 'y-',
            label='Fit of y0 + A * (1 - exp((-t/τ)**α)): \n y0=%5.3f, A=%5.3f, τ=%5.3f, α=%5.3f' % tuple(popt_a))
    ax.grid()
    ax.set_ylabel('Intensity (normalized)', fontsize=16)
    # ax.set_xlabel('Time in s', fontsize=14)
    ax.set_xlim(-5, ndata[-1, 0]+w)
    ax.legend(loc='best', fontsize=14)

    Res.plot(t, Residual_a)
    Res.set_ylabel('Residuals', fontsize=16)
    Res.axhline(y=Res_Mean_a, xmin=0, xmax=60, color='r', linewidth=0.5,
                label='Mean difference: %5.3f' % float(Res_Mean_a))
    Res.legend(loc='best', fontsize=14)
    Res.set_xlabel('Time in s', fontsize=16)
    Res.set_xlim(-5, ndata[-1, 0]+w)
    Res.grid()

    figA.align_ylabels((ax, Res)[:])
    figA.tight_layout()
    figA.show()

    if SaveExpA == True:
        
        if not os.path.isdir(PlotFolder):
            os.makedirs(PlotFolder)
        
        if ComparisonRegion == True:
            if IncludeBleachcontrol ==  True:
                figA.savefig(PlotFolder+file_name[:-4]+'--ROI 4'+'---Fit-Exp-Alpha-Plot_+PhotofadeCorr' + now + '_' + str(x1) + '_' +str(x2) +'.pdf', format="pdf")
            else:
                figA.savefig(PlotFolder+file_name[:-4]+ '--ROI 4'+'---Fit-Exp-Alpha-Plot_noCorr'+ now + '_' + str(x1) + '_' +str(x2) +'.pdf', format="pdf")
        else: 
            if IncludeBleachcontrol ==  True:
                figA.savefig(PlotFolder+file_name[:-4]+'---Fit-Exp-Alpha-Plot_+PhotofadeCorr' + now + '_' + str(x1) + '_' +str(x2) +'.pdf', format="pdf")
            else:
                figA.savefig(PlotFolder+file_name[:-4]+'---Fit-Exp-Alpha-Plot_noCorr'+ now + '_' + str(x1) + '_' +str(x2) +'.pdf', format="pdf")
    
    if SaveExpAtofile == True:
        #include the generated data
        if ComparisonRegion == True: 
            with open(csv_file_path, 'a', newline='') as csvfile:
                    csvwriter = csv.writer(csvfile, delimiter='\t')
                    csvwriter.writerow(['ExpA', IncludeBleachcontrol, Background, ComparisonRegion, alpha, IE_a, t_a_half, d_t_a_half, '-', '-',R , D_ExpA, d_D_ExpA , file_name + '  Comparison region', Loop[x]])
                    #csvfile.write('\n')  # Add a newline
        else:
            with open(csv_file_path, 'a', newline='') as csvfile:
                    csvwriter = csv.writer(csvfile, delimiter='\t')
                    csvwriter.writerow(['ExpA', IncludeBleachcontrol, Background, ComparisonRegion, alpha, IE_a, t_a_half, d_t_a_half, '-', '-',R , D_ExpA, d_D_ExpA , file_name, Loop[x]])
                    #csvfile.write('\n')  # Add a newline
               
"""____________________________________________________________________________
Combined Plot
"""

if CombinedPlot == True:

    figC, ([ax1, ax2, ax3], [Res1, Res2, Res3]) = plt.subplots(2, 3, figsize=(
        30, 7), gridspec_kw={'height_ratios': [4, 1]}) #,layout='constrained'
    figC.suptitle('Comparison of Different Fitting Functions', fontsize=20)

    # ______axes_______
    ax1.set_ylabel('Intensity (normalized)', fontsize=16)
    ax2.set_ylabel('Intensity (normalized)', fontsize=16)
    ax3.set_ylabel('Intensity (normalized)', fontsize=16)
    Res1.set_ylabel('Residuals', fontsize=16)
    Res2.set_ylabel('Residuals', fontsize=16)
    Res3.set_ylabel('Residuals', fontsize=16)
    Res1.set_xlabel('Time in s', fontsize=16)
    Res2.set_xlabel('Time in s', fontsize=16)
    Res3.set_xlabel('Time in s', fontsize=16)
    ax1.set_xlim(-5, ndata[-1, 0]+w)
    ax2.set_xlim(-5, ndata[-1, 0]+w)
    ax3.set_xlim(-5, ndata[-1, 0]+w)
    Res1.set_xlim(-5, ndata[-1, 0]+w)
    Res2.set_xlim(-5, ndata[-1, 0]+w)
    Res3.set_xlim(-5, ndata[-1, 0]+w)

    figC.align_ylabels((ax1, Res1)[:])

    # ______1______
    ax1.set_title('One Exponential Summand')
    ax1.plot(t, exp_decay_one(t, *popt_one), 'b-',
             label='Fit: y0 + A * (1 - exp(-t/τ)): \n y0=%5.3f, A=%5.3f,τ=%5.3f' % tuple(popt_one))
    if IncludeBleachcontrol == True:
        ax1.scatter(ndata[:, 0], ndata[:, 1], marker='.',
                    label='Measurement (photofade corrected)', color='k')
    else:
        ax1.scatter(ndata[:, 0], ndata[:, 1], marker='.',
                    label='Measurement', color='k')
    ax1.grid(True)
    ax1.legend(loc='best', fontsize=14)

    Res1.plot(t, Residual_one)
    Res1.set_xlabel('Time in s')
    Res1.grid(True)
    Res1.axhline(y=Res_Mean_one, xmin=0, xmax=60, color='r', linewidth=0.5,
                 label='Mean difference: %5.4f' % float(Res_Mean_one))
    Res1.legend(loc='best', fontsize=14)

    # ____2_______
    ax2.set_title('Two Exponential Summands')
    ax2.plot(t, exp_decay_two(t, *popt_two), 'g-',
             label='Fit: y0 + A1 *(1 - exp(-t/τ1)) + A2 *(1 - exp(-t/τ2)) : \n y0=%5.3f, A1=%5.3f, \n τ1=%5.3f, A2=%5.3f, τ2=%5.3f' % tuple(popt_two))  # fit
    if IncludeBleachcontrol == True:
        ax2.scatter(ndata[:, 0], ndata[:, 1], marker='.',
                    label='Measurement (photofade corrected)', color='k')
    else:
        ax2.scatter(ndata[:, 0], ndata[:, 1], marker='.',
                    label='Measurement', color='k')
    ax2.grid(True)
    ax2.legend(loc='best', fontsize=14)

    Res2.plot(t, Residual)
    Res2.set_xlabel('Time in s')
    Res2.axhline(y=Res_Mean, xmin=0, xmax=60, color='r', linewidth=0.5,
                 label='Mean difference: %5.4f' % float(Res_Mean))
    Res2.grid(True)
    Res2.legend(loc='best', fontsize=14)

    # _____alpha______
    ax3.set_title('Exponential Function to Power of α')
    ax3.plot(t, exp_decay_a(t, *popt_a), 'y-',
             label='Fit: y0 + A * (1 - exp((-t/τ)**α)): \n y0=%5.3f, A=%5.3f, \n τ=%5.3f, α=%5.3f' % tuple(popt_a))
    if IncludeBleachcontrol == True:
        ax3.scatter(ndata[:, 0], ndata[:, 1], marker='.',
                    label='Measurement (photofade corrected)', color='k')
    else:
        ax3.scatter(ndata[:, 0], ndata[:, 1], marker='.',
                    label='Measurement', color='k')
    ax3.grid(True)
    ax3.legend(loc='best', fontsize=14)

    Res3.plot(t, Residual_a)
    Res3.set_xlabel('Time in s')
    Res3.grid(True)
    Res3.axhline(y=Res_Mean_a, xmin=0, xmax=60, color='r', linewidth=0.5,
                 label='Mean difference: %5.4f' % float(Res_Mean_a))
    Res3.legend(loc='best', fontsize=14)

    figC.show()

    if SaveCombindedPlot == True:
        if not os.path.isdir(PlotFolder):
            os.makedirs(PlotFolder)
        if IncludeBleachcontrol ==  True:
            figC.savefig(
                PlotFolder+file_name[:-4]+'---Fit-Plots-Combined_+PhotofadeCorr' + now + '_' + str(x1) + '_' +str(x2) + '.pdf', format="pdf")
        else:
            figC.savefig(
                PlotFolder+file_name[:-4]+'---Fit-Plots-Combined_noCorr'+ now + '_' + str(x1) + '_' +str(x2) +'.pdf', format="pdf")


"""____________________________________________________________________________
Dictonary with all Results
"""
if TwoExp == True:
    ResultDict_TwoExp = {'Exp.function 2nd order': {'IE': IE,
                                             'I_immobile': I_immobile,
                                             'y0': [y_2, d_y_2],
                                             'F1 [%]': F_1,
                                             'A1': [A_1, d_A_1],
                                             't1': [t_1, d_t_1],
                                             't1_half': t_1_half,
                                             'k1': k_1,
                                             'F2 [%]': F_2,
                                             'A2': [A_2, d_A_2],
                                             't2': [t_2, d_t_2],
                                             't2_half': t_2_half,
                                             'k2': k_2,
                                             'F_mobile [%]': F_mobile,
                                             'F_immobile [%]': F_immobile}
                  }
    
if OneExp == True:
    ResultDict_OneExp = {'Exp.function 1st order': {'IE': IE_1,
                                             'I_immobile': I_immobile_0,
                                             'y0': [y_1, d_y_1],
                                             'A1': [A_0, d_A_0],
                                             't1': [t_0, d_t_0],
                                             't_half': t_0_half,
                                             'k1': k_0,
                                             'F_mobile [%]': F_mobile_0,
                                             'F_immobile [%]': F_immobile_0}
                  
                  }

if ExpA == True:
    ResultDict_ExpA = {'Exp.function with power to a': {'IE': IE_a,
                                                   'I_immobile': I_immobile_a,
                                                   'y0': [y_a, d_y_a],
                                                   'A': [A_a, d_A_a],
                                                   't0': [t_a, d_t_a],
                                                   'α': [alpha, d_alpha],
                                                   't_half': t_a_half,
                                                   'k1': k_a,
                                                   'F_mobile [%]': F_mobile_a,
                                                   'F_immobile [%]': F_immobile_a}
                  }

##############################################################################

print('___________')
print('done')
