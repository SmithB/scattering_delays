import numpy as np
from pathlib import Path

def ice_n_of_lambda(λ):

   path = str(Path(__file__).parent.absolute())

   WB08=np.loadtxt(path+'/IOP_2008_ASCIItable.dat')
   log_lambda_table=np.log(WB08[:,0]*1.e-6)
   n=np.exp( np.interp(np.log(λ), log_lambda_table, np.log(WB08[:,1])))\
       + 1j*np.exp(np.interp(np.log(λ), log_lambda_table, np.log(WB08[:,2])));
   return n

