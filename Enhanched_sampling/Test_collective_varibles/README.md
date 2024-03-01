#################################################################### Test code for enhanced sampling of aniline + anhydride ##########################################################


####TEST SYSTEM######
aniline 

#####PATHS##### - Reactants####
/dfs7/elizaml4/sebasatc/3_nano_gap_regan/4_Au_supercell_surface_with_all_chemicical_reactants/2_AIMD_runs/4_222_Au_super_cell/1_aniline/20_angstrom_vacuum

#####PATHS#####  - Products#### 
/dfs7/elizaml4/sebasatc/3_nano_gap_regan/5_Au_supercell_surface_with_all_chemicical_products/2_AIMD_runs/3_222_Au_cell_system/1_aniline/1_phenylpentanamide_on_top_gold_surface/3_both_surafce_with_aniline_and_carboxylic/2_Experimental_gap_test/1_20A_gap_space_with_relax_20A_gap_anhydride

#####NOTES####
The Collective variables is d (C*-X*)
where C* is the Crabon atom where the reaction happens
where X* is the atom that is being reacted, for aniline + anhydride X* = Nitrogen
where d is the distance in A from C* and X*

The second Collective variables is d (C*- O)
where C* is the carbon where the reactions happens 
where O is the oxygen atom of the carboxlic acid that is formed
where d is the dustance in A from C* to O


####More Notes####
Figure 1a shows the plot of C*-X* and C*-O from XDATCAR_1
where XDATCAR_1 = AIMD of 2x2x2 Au surafce with 20A gap of aniline with anhydride and of XDATCAR_2, where XDATCAR_2 = same system but products. 

Figure 1b is to confirm periodic boundaries

Figure 1c is the chemical mechanism of aniline + anhydride 


####### Final Notes#####

Contact Sebastian at sebasatc@uci.edu for questions



############################################################################### THANK YOU ############################################################################################
