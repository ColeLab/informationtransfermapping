# Taku Ito
# 09/1/2016
# Script to transfer python's 'importRuleTimingsV3' function to a Matlab readable format

import scripts3_functions as func
import numpy as np

outdir = '/projects2/ModalityControl2/data/CPROTaskIdentifiers/'

subjNums = ['032', '033', '037', '038', '039', '045', 
            '013', '014', '016', '017', '018', '021', 
            '023', '024', '025', '026', '027', '031', 
            '035', '046', '042', '028', '048', '053', 
            '040', '049', '057', '062', '050', '030', '047', '034']
ruledims = ['logic','sensory','motor']

for ruledim in ruledims:
    for subj in subjNums:
        rules, rulesmb = func.importRuleTimingsV3(subj, ruledim, hrfconv=True)

        rulembarraymatlab = []
        for rule in rulesmb:
            tmp = np.asarray(rulesmb[rule].keys())
            # Add a 1 to it for matlab indices
            tmp += 1

            # Add to array
            rulembarraymatlab.append(tmp)
            outarray = np.asarray(rulembarraymatlab)
            # Tranpose to miniblocks x rule
            outarray = outarray.T

            # Save outfile
            outname = subj + '_' + ruledim + '_miniblocksPerRule_MatlabIndices.csv'
            np.savetxt(outdir + outname, outarray, delimiter=',')
