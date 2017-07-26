# Taku Ito
# function run a single permutation on NM3 and save out file
# 04/25/2017

import sys
import numpy as np
import scipy.stats as stats
import multiprocessing as mp
import os
os.environ['OMP_NUM_THREADS'] = str(1)
import warnings
warnings.filterwarnings('ignore')
import time

import networkinformationtransfer as n2n

# Set basic parameters
datadir = '/work/ti61/projects/ModalityControl/data/networkITE_data/'
runLength = 4648
subjNums = ['032', '033', '037', '038', '039', '045', 
            '013', '014', '016', '017', '018', '021', 
            '023', '024', '025', '026', '027', '031', 
            '035', '046', '042', '028', '048', '053', 
            '040', '049', '057', '062', '050', '030', '047', '034']

# Load in network array
networkdef = np.loadtxt(datadir + 'network_array.csv', delimiter=',')


# Load in network keys (each network associated with a number in network array)
networkmappings = {'fpn':7, 'vis':1, 'smn':2, 'con':3, 'dmn':6, 'aud1':8, 'aud2':9, 'dan':11}
# Force aud2 key to be the same as aud1 (merging two auditory networks)
aud2_ind = np.where(networkdef==networkmappings['aud2'])[0]
networkdef[aud2_ind] = networkmappings['aud1']
# Redefine new network mappings with no aud1/aud2 distinction
networkmappings = {'fpn':7, 'vis':1, 'smn':2, 'con':3, 'dmn':6, 'aud':8, 'dan':11}

def runPermutations(cycle):
    """
    cycle - indicates which permutation cycle. Each cycle runs 100 permutations
    """
    # Specify output directory
    outdir = '/work/ti61/projects/ModalityControl/data/results/NetworkToNetworkITE/'

    ####
    # Load in resting-state FC data
    fcmat = {}
    for subj in subjNums:
        datafile = datadir + 'FC_Estimates/'
        fcmat[subj] = np.loadtxt(datafile + subj + '_multregconn_restfc.csv', delimiter=',')

    ####
    permutations = np.arange(100*cycle,100*cycle+100)

    # Make sure it's random!
    randomseed = np.random.randint(0,100000000,size=(len(permutations)*len(subjNums),1)) 
    seedcount = 0
    for nPermutation in permutations:
        # Begin multiprocessing scheme
        tstart = time.time()
        inputs = []
        for subj in subjNums:
            tmp_mat = fcmat[subj].copy()
            inputs.append((subj,tmp_mat,randomseed[seedcount]))
            seedcount += 1

        cpucount = mp.cpu_count()
        pool = mp.Pool(processes=cpucount-1)
        results = pool.map_async(informationTransferMappingWrapper,inputs).get()
        pool.close()
        pool.join()

        ####
        # Collect results from multiprocessing
        # Instantiate empty matrices
        ruledims = ['logic','sensory','motor']
        ite_matrix = {}
        for ruledim in ruledims:
            ite_matrix[ruledim] = np.zeros((len(networkmappings),len(networkmappings),len(subjNums)))

        # Fill out matrix with results
        scount = 0
        for result in results:
            for ruledim in ruledims:
                ite_matrix[ruledim][:,:,scount] = result[ruledim]
                
            scount += 1
            
        tfinish = time.time()
        print 'Completed permutation', nPermutation, 'in', tfinish-tstart, 'seconds'

        ####
        # Create dictionary that reflects network ordering for matrix rows and columns
        netkeys = {0:'vis',1:'smn',2:'con',3:'dmn',4:'fpn', 5:'aud', 6:'dan'}
        num_networks=len(netkeys)

        ####
        # Compute average rho's across subjects
        outdir = '/work/ti61/projects/ModalityControl/data/results/NetworkToNetworkITE/'
        avg_rho = {}
        for ruledim in ruledims:
            avg_rho[ruledim] = np.zeros((num_networks,num_networks))
            
            for net1 in netkeys:
                
                for net2 in netkeys:
                    # Skip if net1 and net2
                    if net1==net2: 
                        continue
                        
                    # Store results
                    avg_rho[ruledim][net1,net2] = np.mean(ite_matrix[ruledim][net1,net2,:])
                    

            np.savetxt(outdir + 'NetworkITE_NullModel_' + ruledim + '_Permutation' + str(nPermutation) + '.csv', avg_rho[ruledim], delimiter=',')


def informationTransferMappingWrapper((subj,fcmat,seedcount)):
    """
    A wrapper so we can use multiprocessing to run subjects in parallel
    """
    out = n2n.networkToNetworkInformationTransferMapping(subj,fcmat,null=True,seed=seedcount)
    return out 
