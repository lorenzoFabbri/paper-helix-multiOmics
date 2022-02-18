#!/bin/bash

choice="$1"

# GGMs
if [[ "$choice" == "GGM" ]]
then
scp -r "/home/lorenzo/Documents/university/PhD/papers/paper1_helix_multiOmics/paper-helix-multiOmics/results/ggm/" lfabbri@172.20.10.115:/PROJECTES/HELIX/lorenzoF/paper-helix-multiOmics/results
fi

# Relevance Networks
if [[ "$choice" == "RN" ]]
then
scp -r "/home/lorenzo/Documents/university/PhD/papers/paper1_helix_multiOmics/paper-helix-multiOmics/results/networks/relevance/" lfabbri@172.20.10.115:/PROJECTES/HELIX/lorenzoF/paper-helix-multiOmics/results/networks
fi

# All output data
if [[ "$choice" == "data" ]]
then
scp -r "/home/lorenzo/Documents/university/PhD/papers/paper1_helix_multiOmics/paper-helix-multiOmics/data" lfabbri@172.20.10.115:/PROJECTES/HELIX/lorenzoF/paper-helix-multiOmics/
fi

# Metadata
if [[ "$choice" == "meta" ]]
then
scp -r /home/lorenzo/Documents/university/PhD/papers/paper1_helix_multiOmics/paper-helix-multiOmics/data/*.csv lfabbri@172.20.10.115:"/PROJECTES/HELIX/HELIX_analyses/Omics_variability/Final\ subsets/Final\ subsets\ v2/data_lorenzo/"
fi
