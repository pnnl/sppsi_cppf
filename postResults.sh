#!/bin/bash

# rm -rf postProc
rhos='rho_sgl_edge_pos_mobile,rho_sgl_edge_neg_mobile,rho_sgl_screw_pos_mobile,rho_sgl_screw_neg_mobile,rho_sgl_edge_pos_immobile,rho_sgl_edge_neg_immobile,rho_sgl_screw_pos_immobile,rho_sgl_screw_neg_immobile,rho_dip_edge,rho_dip_screw'

# postResults --co $rhos --separation x,y,z --split 20grains16x16x16_shearXY.spectralOut
cd postProc
# for i in {1..12}
# do
#  addCalculation --label "$i_rho" --formula "#$i_rho_sgl_edge_pos_mobile# + #$i_rho_sgl_edge_neg_mobile# + #$i_rho_sgl_screw_pos_mobile# + #$i_rho_sgl_screw_neg_mobile# + #$i_rho_sgl_edge_pos_immobile# + #$i_rho_sgl_edge_neg_immobile# + #$i_rho_sgl_screw_pos_immobile# + #$i_rho_sgl_screw_neg_immobile# + #$i_rho_dip_edge# + #$i_rho_dip_screw#" 20grains16x16x16_shearXY_inc10.txt
# done

# commands to calculate ss curve
#### cd postProc
### addCauchy *.txt
### addStrainTensors --left --logarithmic *.txt
### addMises -s Cauchy *.txt
### addMises -e "ln(V)" *.txt


fname='20grains16x16x16_shearXY_inc10.txt'
vtk_rectilinearGrid $fname
vtk_addRectilinearGridData --data $rhos --vtk *.vtr $fname

# results for ss curve
mkdir postProc_ss
postResults --cr f,p -d postProc_ss *.spectralOut

# results for dislocation density
mkdir postProc_dislo
postResults --co edge_density,dipole_density -d postProc_dislo --separation x,y,z --split *.spectralOut


