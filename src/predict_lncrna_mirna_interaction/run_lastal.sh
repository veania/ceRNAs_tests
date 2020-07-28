#!/bin/bash

HOME_DIR=/home/aelizarova
PROJECT_DIR=${HOME_DIR}/ceRNAs/ceRNAs_tests
LASTDB=${PROJECT_DIR}/data/human_mir_db
#LNC_SEQ_PATH=${PROJECT_DIR}/data/lncRNAs_U.fasta
LNC_SEQ_PATH=${PROJECT_DIR}/"out/KD.gene.ID|ASO_id_U.fasta"
OUT_DIR=${PROJECT_DIR}/out/lnc_mir_prediction

[ -d ${OUT_DIR} ] || mkdir ${OUT_DIR}
lastal -v -f tab ${LASTDB} ${LNC_SEQ_PATH} -l 1 -m 4 -d 4 -e 10 > ${OUT_DIR}/myalns.tab

echo 'lastal finsished calculation'