# Taku Ito
# Script that loads Glasser parcels into python
# Also loads KK's network partitions into python

import numpy as np
import nibabel as nib

def loadGlasserParcels():
    basedir = '/projects2/ModalityControl2/data/GlasserKKPartition/'
    filename = 'Q1-Q6_RelatedParcellation210.LR.CorticalAreas_dil_Colors.32k_fs_LR.csv'
    parcels = np.loadtxt(basedir + filename, delimiter=',')
#    parcels = nib.load(basedir + filename)
#    parcels = parcels.get_data()
#    parcels = np.squeeze(parcels)
    return parcels


def loadGlasserNetworks():
    basedir = '/projects2/ModalityControl2/data/GlasserKKPartition/'
    filename = 'Glasser_KK_Partitions.csv'
    #filename = 'Glasser_KK_Partitions_old.csv'

    networks = np.loadtxt(basedir + filename, delimiter=',')
    return networks

