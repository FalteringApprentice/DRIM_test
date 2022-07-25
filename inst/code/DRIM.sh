#!/bin/bash
echo "resolution ratio : $1"
echo "Column name of cell type : $2"
echo "Expansion threshold : $3"
python sc_st_gene_charge.py 
python spot_pre.py $1 
python GS_mapping_HVG_gene.py $1 $2
python newRegineGrowing_use_RG.py $1 $3 
python comb_mapping_spot.py $1 
python it_final_celltype.py $1 $2