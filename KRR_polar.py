#!/usr/bin/python3
import sys
import numpy as np
import time
import matplotlib.pyplot as plt
import scipy.stats as stats
from scipy.special import erf
import scipy.spatial
import os,argparse
import random



def file_len(fname):
  '''
  This function gives the number of lines of a given text file
  '''
  with open(fname) as f:
      for i, l in enumerate(f):
          pass
  return i + 1


def read_txt(name,array):
  '''
  The following function reads a text file containing alternatively one line of comment and one line of data, and returns a numpy array
  '''
  name = os.path.join('data', name)
  N_points = int(file_len(name)/2)
  print ('There are', N_points, 'points in the file', name)
  data = open(name,'r')
  array = []
  counter = 0
  while (counter < N_points):
    line = data.readline() # Skip line of comment
    if not line:
      break
    line = data.readline() # Line containing the information wanted
    array.append(list(map(float, line.split()[0:])))
    counter = counter+1
  data.close()
  array = np.array(array)
  return array

# The following function reads a numpy file (.npy) containing the quantity of interest stored directly as an array.

def read_npy(name,array):
  '''
  The following function reads a numpy array.
  '''
  name = os.path.join('data', name)
# Load the saved numpy array
  array_reloaded = np.load(name)
# Determine length of the array (corresponds to the number of training (or test) points)
  N_points = len(array_reloaded)
  print ('There are', N_points, 'points in the file', name)
  return array_reloaded

def read_data(name,array):
  '''
  Reads data either as text or numpy format
  '''
  extension = os.path.splitext(name)[1]
  if (extension == '.dat'):
    array = read_txt(name,array)
  elif (extension == '.npy'):
    array = read_npy(name,array)
  else:
    print ('Extension for %s unrecognized !'.format(name))
    sys.exit()
  return array


def shuffle_data(array,index_shuf):
  '''
  Shuffle an array according to the list of indexes given in input
  '''
  array_shuf = []
  for i in index_shuf:
      array_shuf.append(array[i])
  return array_shuf

def select_data(array,N):
  '''
Selects only the first N values of the input array
  '''
  array = array[0:N]
  array = np.array(array)
  return array

def print_error(polar_dfpt, polar_ML,f_,f_polar,f_Raman_predict,f_Raman_dfpt,N_training):
  '''
  Print error on screen and write it to file
  '''
  mae = np.mean(abs(polar_dfpt - polar_ML),axis=0) # MAE
  max_error = np.amax(abs(polar_dfpt - polar_ML),axis=0)
  min_error = np.amin(abs(polar_dfpt - polar_ML),axis=0)
  rmse = np.sqrt(np.mean((polar_dfpt - polar_ML)**2,axis=0)) # RMSE
  std = np.std(polar_dfpt,axis=0) # STD
  rmse_normal = 100*rmse/std
  print ("Error                     |      xx          yy          zz          xy          xz          yz")
  print ("--------------------------|---------------------------------------------------------------------")
  print ("MAE                       |", end=' ')
  for a in mae: print ("{0:.4e}".format(a),end=' ')
  print ('')
  print ("Max                       |",end=' ')
  for a in max_error: print ("{0:.4e}".format(a),end=' ')
  print ('')
  print ("Min                       |",end=' ')
  for a in min_error: print ("{0:.4e}".format(a),end=' ')
  print ('')
  print ("RMSE/STD                  |", end= ' ')
  for a in rmse_normal: print ("{0:.4e}".format(a), end='  ')
  print ('')

  # Write to file the error versus sigma and Lambda
  f_.write( ("{0:.2e}".format(sigma)) + " " + "{0:.2e}".format(Lambda) + " " +("{:d}".format(N_training)) + " " )
  for a in mae: f_.write( ("{0:.5e}".format(a)) + " " )
  for a in rmse_normal: f_.write( ("{0:.5e}".format(a)) + " " )
  f_.write("\n")

  f_polar.write("# num | polar_DFPTxx | predictionxx | polar_DFPTyy | predictionyy | ... ")
  for a in rmse_normal: f_polar.write( ("{0:.5e}".format(a)) + " " )
  f_polar.write("\n")
  for i in range(len(polar_dfpt)):
    f_polar.write(str(i)+" ")
    for comp in range(6):
      f_polar.write(str(polar_dfpt[i,comp]) + " " + str(polar_ML[i,comp]) + " ")
    f_polar.write("\n")

  for i in range(len(polar_dfpt)):
    f_Raman_predict.write(str(polar_ML[i,0]) + " " + str(polar_ML[i,1]) + " " + str(polar_ML[i,2]) + " " + str(polar_ML[i,3]) + " " + str(polar_ML[i,4]) + " " + str(polar_ML[i,5]) + '\n' )

  # Write the corresponding true spectrum (i.e. DFPT on training set) 
  for i in range(len(polar_dfpt)):
          f_Raman_dfpt.write(str(polar_dfpt[i,0]) + " " + str(polar_dfpt[i,1]) + " " + str(polar_dfpt[i,2]) + " "  + str(polar_dfpt[i,3]) + " " + str(polar_dfpt[i,4]) + " " + str(polar_dfpt[i,5]) + '\n' )

  return mae


def compute_Raman(name_raman_dfpt,name_raman_ML):
  '''
  Computes polarizability autocorrelation function and Raman spectrum
  '''
  name_list = "list.dat" # Input file for autocorrelation script
  os.chdir("predictions")
  
  with open(name_list,'w') as f_:
    f_.write(name_raman_dfpt + ' ') # The space is needed so that the autocorrelation script parses the name correctly
  # Calculate autocorrelation and Raman spectrum of dfpt
  os.system("python3 autocorr.py control.autocorr.in")
  
  with open(name_list,'w') as f_:
    f_.write(name_raman_ML + ' ')
  # Calculate autocorrelation and Raman spectrum of prediction
  os.system("python3 autocorr.py control.autocorr.in")
  
  # Move Raman data to the "plots" directory
  os.chdir("..")
  os.rename('predictions/Raman_'+name_raman_dfpt, 'plots/Raman_'+name_raman_dfpt)
  os.rename('predictions/Raman_'+name_raman_ML, 'plots/Raman_'+name_raman_ML)


# Define maximum and minimum values for hyperparameters in case a grid search is desired
Lambda_min = 1.e-8
Lambda_max = 1
sigma_min = 1.e-4
sigma_max = 1
sigma = sigma_min # Current sigma
Lambda = Lambda_min # current Lambda
sigma_opt = sigma # Optimal sigma
Lambda_opt = Lambda # Optimal lambda

def vary_hyper(sigma,Lambda,sigma_min,sigma_max,Lambda_min,Lambda_max):
  '''
  Vary hyperparameters when a grid search is requested
  '''
  if (Lambda<Lambda_max):
    Lambda=Lambda*10**1
    repeat = True
  elif (sigma<sigma_max):
    #sigma = sigma*10**0.25
    sigma = sigma*10**0.5
    Lambda = Lambda_min
    repeat = True
  else:
    repeat = False
  return sigma,Lambda,repeat
  

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description = ''' This program will interpolate polarizability tensors based on their geometrical features (atomic coordinates, atomic density, etc.) 
All specifications are given in the \'control_KRR\' file.
  It requires at least 2 files for the training set:
   - 1 file containing the geometrical information of your system for each structure (i.e., atomic coordinates, atomic densities, etc.)
   - 1 file containing the polarizabilities (or dipoles, etc.) corresponding to the aforementioned structures 
  and 1 file for the test set:
   - 1 file containing the geometrical information of your system for each structure, different from the training set
  After training, the routine will predict the polarizability tensors of the new (test) set of structures.
  The script accepts data as text (\'.dat\' files) or numpy format (\'.npy\' files). In the case of text data, it is expected that there is 1 alternatively 1 line of comment, and 1 line of data, with the elements to be predicted being separated by spaces
  ''',formatter_class=argparse.RawTextHelpFormatter)

  args = parser.parse_args()


  # I) Preparation
  
  baselining = False

  # a) Read control file
  with open("control_KRR",'r') as control_file:
    for line in control_file:
       if "features_training " in line:
         file_features_training = line.split()[-1]
       if "dfpt_training " in line:
         file_dfpt_training = line.split()[-1]
       if "polar_mol_training " in line:
         file_molpol_training = line.split()[-1]
         baselining = True
       if "features_test " in line:
         file_features_test = line.split()[-1]
       if "dfpt_test " in line:
         file_dfpt_test = line.split()[-1]
       if "polar_mol_test " in line:
         file_molpol_test = line.split()[-1]
         baselining = True
       if "length_training " in line:
         N_training = int(line.split()[-1])
       if "length_test " in line:
         N_test = int(line.split()[-1])
       if "grid_search " in line:
         if str(line.split()[-1])=="yes":
           grid_search = True
         elif str(line.split()[-1])=="no":
           grid_search = False
         else:
           print("What do you want in life?\nMake a choice and come back when you are ready.\nExiting (bad keyword for \'grid_search\').")
           sys.exit()
       if "plot " in line:
         if str(line.split()[-1])=="yes":
           plot_Raman=True
         elif str(line.split()[-1])=="no":
           plot_Raman=False
         else:
           print("What do you want in life?\nMake a choice and come back when you are ready.\nExiting (bad keyword for \'plot\').")
           sys.exit()
  
  # b) Read necessary files for the training set
  
  print ("You have chosen", N_training, "data points in the training set")
  

  # Read descriptors
  u_training = []
  u_training = read_data(file_features_training,u_training)
  if (len(u_training) < N_training):
    print ("The number of training points you have asked for is larger than the number of data points available \n Exiting")
    sys.exit()

  # Shuffle data set. The following will produce a list of random indices of length the total size of the training set file
  index_shuf = np.array(list(range(len(u_training))))
  np.random.shuffle(index_shuf)
  
  # Use same indexes for all calculations, this way when increasing N_train, we always add new training points
  #index_shuf = np.genfromtxt("index2.dat", usecols = 0, dtype=int,unpack=True)

  # Print indices used for this run
  f_indices=open("indices_Ntrain"+str(N_training)+'.dat','w')
  for a in index_shuf[0:N_training]: 
    f_indices.write( ("{:d}".format(a)) + "\n" )
  f_indices.close()

  #np.random.shuffle(index_shuf)

  # First shuffle the data... The same shuffling has to be applied to all the training data, not just the descriptors
  u_training = shuffle_data(u_training,index_shuf)
  #...then select the first N_training values, as requested in the control file
  u_training = select_data(u_training,N_training)
  
  # Read training DFPT polarizabilities
  polar_dfpt_training = []
  polar_dfpt_training = read_data(file_dfpt_training,polar_dfpt_training)
  polar_dfpt_training = shuffle_data(polar_dfpt_training,index_shuf)
  polar_dfpt_training = select_data(polar_dfpt_training,N_training)
  mean_dfpt_training = np.mean(polar_dfpt_training,axis=0)
  print ("Mean DFPT polarizability (training): ", mean_dfpt_training)
  if (len(polar_dfpt_training) < N_training):
    print ("The number of training points you have asked for is larger than the number of data points available \n Exiting")
    sys.exit()

  # Read sum of molecular polarizabilities (used for baselining, does not enter kernel)
  if baselining:
    molpol_training = []
    molpol_training = read_data(file_molpol_training,molpol_training)
    molpol_training = shuffle_data(molpol_training,index_shuf)
    molpol_training = select_data(molpol_training,N_training)
    if (len(molpol_training) < N_training):
      print ("The number of training points you have asked for is larger than the number of data points available \n Exiting")
      sys.exit()
    mean_molpol_training = np.mean(molpol_training,axis=0)
    print ('Average sum of molecular polarizabilities (training): ', mean_molpol_training)

  
  # c) Read necessary files for the test set
  
  print ("You chose", N_test, "data points to extrapolate.")
  
  # Read descriptors 
  u_test = []
  u_test = read_data(file_features_test,u_test)
  if (len(u_test) < N_test):
    print ("The number of test points you have asked for is larger than the number of data points available \n Exiting")
    sys.exit()
  # Select first N_test points
  u_test = select_data(u_test,N_test)

  
  # Read the DFPT polarizabilities (only used to calculate the error we make with our predictive model, in principles we do not have access to them)
  polar_dfpt_test = []
  polar_dfpt_test = read_data(file_dfpt_test,polar_dfpt_test)
  polar_dfpt_test = select_data(polar_dfpt_test,N_test)
  mean_dfpt_test = np.mean(polar_dfpt_test,axis=0)
  print ("Mean DFPT polarizability (test): ", mean_dfpt_test)
  if (len(polar_dfpt_test) < N_test):
    print ("The number of test points you have asked for is larger than the number of data points available \n Exiting")
    sys.exit()

  # Read sum of molecular polarizabilities (used for baselining)
  if baselining:
    molpol_test = []
    molpol_test = read_data(file_molpol_test,molpol_test)
    molpol_test = select_data(molpol_test,N_test)
  

  # II) Now the real work will start: start constructing the kernel, weights, etc. for a given Lambda and sigma
 
  # Open files to store errors and predictions 

  folders = ["errors","predictions","plots"]
  for d_ in folders:
    if not os.path.exists(d_):
      os.mkdir(d_)

  f_training=open("errors/error_training_"+str(N_training)+'.dat','w')
  f_test=open("errors/error_Ntest"+str(N_test)+"_Ntrain"+str(N_training)+'.dat','w')
  f_polar_training=open("predictions/predict_polar_training_"+str(N_training)+'.dat','w')
  f_polar_test=open("predictions/predict_polar_test_Ntest"+str(N_test)+"_Ntrain"+str(N_training)+'.dat','w')
  f_Raman_predict_training=open("predictions/polar_ML_training_"+str(N_training)+'.dat','w')
  f_Raman_dfpt_training=open("predictions/polar_dfpt_training_"+str(N_training)+'.dat','w')
  name_raman_dfpt = "polar_dfpt_test_"+str(N_test)+".dat"
  name_raman_ML = "polar_ML_test_Ntest"+str(N_test)+"_Ntrain"+str(N_training)+".dat"
  f_Raman_predict_test=open("predictions/" + name_raman_ML,'w')
  f_Raman_dfpt_test=open("predictions/" + name_raman_dfpt,'w')
  
  error_old = -1 # This will serve to determine the optimal hyperparameter
  
  repeat = True # Says whether to make another prediction with updated hyperparameters
  first = True # First pass to the loop
  while(repeat == True):
    if(grid_search and first):
      first = False
    elif(grid_search and not first):
      sigma,Lambda,repeat = vary_hyper(sigma,Lambda,sigma_min,sigma_max,Lambda_min,Lambda_max)
    else: # Standard values
      sigma = 3e-3 # Good value for the paracetamol molecule with a few hundreds of training points 
      #sigma = 1e-2 
      Lambda = 1e-5 
      repeat = False

    # 1) Now predict the polarizability tensors of each training structure. The predicted value should almost be exact, but some deviation is allowed.
    # The optimal lambda and sigma have to be determined against a validation test.
    
    # Construct the square-exponential kernel matrix, which has dimensions N_training*N_training
    # k(i,j)=exp(-|u_i-u_j|^2/(2*sigma^2)) with u containing geometrical features
    
    print ("Constructing kernel for the training set...")
    
    start = time.time() 
    t1 = time.time()
    
    eucl = scipy.spatial.distance.pdist(u_training[:,:],metric='sqeuclidean') # Dense matrix containing the norms 
    eucl = scipy.spatial.distance.squareform(eucl) # Put the previous matrix back in square form, with redundancies. It is of size N_training*N_training, and its elements are (u_ai-u_bi)**2, where u_ai is the ith feature of the ath configuration
    
    t2 = time.time()
    print('Took ',t2-t1, 'seconds')
    
    # Build the kernel
    k_training = np.exp(-eucl/(2*sigma**2))
    
    t2 = time.time()
    print ("Took ", t2-t1, "seconds")
    
    # Calculate the weights that will be used to "predict" polarizabilities:
    
    print ("Inverting matrix...\n")
    
    L = Lambda*np.identity(N_training)
    inv = np.linalg.inv(k_training+L)
    
    print ("Took ", time.time()-start, "seconds")
    
    print ("Calculating weights...\n")
    
    if (not baselining):
      weights = np.dot(inv,(polar_dfpt_training - mean_dfpt_training))
    else:
      weights = np.dot(inv,(polar_dfpt_training - mean_dfpt_training) - (molpol_training - mean_molpol_training) )
    
    print ("Took ", time.time()-start, "seconds")
    
    print ("Predicting training set polarizabilities (although we have trained on the same points)\n")
    print ("lambda = ", Lambda, ", sigma = ", sigma)
    
    # alpha_ML = alpha_mean_training + sum_1^N w_l k(u,u_l)
    if (not baselining):
      polar_training = np.dot(k_training,weights) + mean_dfpt_training # Add mean only if it was retrieved from the weights
    else:
      polar_training = np.dot(k_training,weights) + mean_dfpt_training + (molpol_training-mean_molpol_training)
  
    # Calculates error
    print_error(polar_dfpt_training, polar_training,f_training,f_polar_training,f_Raman_predict_training,f_Raman_dfpt_training,N_training)
  
    # 2) Predict data for structures the model has not seen before (test set)
    print ("Predict polarizabilities for the test set: be ready, it's starting!")
    
    # Calculate the "kernel". Note that it is a priori not a square matrix here (we don't need to invert it afterwards), since it has dimensions N_test x N_training, with N_test the number of points to be predicted
    
    print ("Constructing \"kernel\"...")
    
    t1 = time.time()
    
    eucl = scipy.spatial.distance.cdist(u_test[:,:],u_training[:,:],metric='sqeuclidean') # Dense matrix containing the norms
    # Build kernel
    k_test = np.exp(-eucl/(2*sigma**2))
    
    t2 = time.time()
    print('Took ',t2-t1, 'seconds')
    
    print ("Predicting values...")
    
    if (not baselining):
      polar_test = np.dot(k_test,weights) + mean_dfpt_training # Add mean only if it was retrieved from the weights
    else:
      polar_test = np.dot(k_test,weights) + mean_dfpt_training + (molpol_test - mean_molpol_training) 
    
    # Print the error on the test set
    
    mae = print_error(polar_dfpt_test, polar_test,f_test,f_polar_test,f_Raman_predict_test,f_Raman_dfpt_test,N_training)
  
    # Find optimal hyperparameters
    if(mae[0] < error_old or error_old == -1): # Check only xx at the moment
      sigma_opt = sigma
      Lambda_opt = Lambda
      error_old = mae[0]
  
  if (grid_search):
    print("Optimal hyperparameters (sigma,lambda) = ({:.2e},{:.2e})".format(sigma_opt,Lambda_opt))

  # Close all files
  f_training.close()
  f_test.close()
  f_polar_training.close()
  f_polar_test.close()
  f_Raman_predict_training.close()
  f_Raman_dfpt_training.close()
  f_Raman_predict_test.close()
  f_Raman_dfpt_test.close()
  

  # Calculate Raman spectrum if requested
  if (plot_Raman and not grid_search):
    print("\nNow computing Raman spectra\n")
    compute_Raman(name_raman_dfpt,name_raman_ML)

