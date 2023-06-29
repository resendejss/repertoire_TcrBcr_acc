#!/bin/bash
PATH_SAMPLES="/media/jean/hd_130gb/Projetos/BigDataIPPPP/diversidadeTcrBcr/TCGA-ACC/4-extractionTcrBcr/outputTrust4"

while read SAMP \n
    do
    python3 trust-stats.py -r $PATH_SAMPLES/${SAMP} > output_trust-stats/${SAMP}
    done < samplesNames.txt


