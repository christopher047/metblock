

###rename tables for publication ###
## tables ###
file.copy("../tables/only_in_blocks.tab", "../tables/table4.tab", overwrite=T) 

## figures ###
file.copy("../figures/region_density.pdf", "../figures/figure2.pdf", overwrite=T) 


### sup2 is segments ###
file.copy("../results/segs_trim.tab", "../supplementary/supp_file3_segs.tab", overwrite=T) 
### sup3 is dmr_seq genesets ###
file.copy("../tables/dmrseq_genrich.tab", "../supplementary/supp_file4_dmrseq_gs.tab", overwrite=T)  
### sup4 is metilene genesets ###
file.copy("../tables/metiline_genrich.tab", "../supplementary/supp_file5_metilene_gs.tab", overwrite=T)
### sup5 is blocks genesets ###
file.copy("../tables/blocks_genrich.tab", "../supplementary/supp_file6_blocks_gs.tab", overwrite=T)
### sup6 blocks ###
file.copy("../results/raw_blocks_trim.tab", "../supplementary/supp_file8_blocks_stats.tab", overwrite=T)
### sup7 is index###
file.copy("../results/raw_index_trim.tab", "../supplementary/supp_file9_block_index.tab", overwrite=T) 
### sup8 is regions comparison ####
file.copy("../figures/test_01.pdf", "../supplementary/supp_figure_1.pdf", overwrite=T)  
### sup9 is 30 regions ###
### pdfunite ../figures/test_*.pdf ../supplementary/supp_file9_regions30.pdf 
### sup 10 is overview of CpGs found per sample ###
file.copy("../raw/sample_overview.tab", "../supplementary/supp_table_samples.tab", overwrite=T) 




