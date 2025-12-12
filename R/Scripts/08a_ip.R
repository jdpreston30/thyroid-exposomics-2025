#!!!!! Make version with largest fragment included as well



#!!!!! CHECK QUANT FRAG, if pass then level 1
#- 8.4.7: TEEP (CP3182)
#- 8.4.9 : Prosulfuron (CP2365)



#!!!! Consider making level 2 but move to smaller x window
#- 8.4.17: o-Cresol (CP3066)
# potential biomatrix effect
aps$CP3066
F2_S2_CP3066

#!!! Level 2 if quant fragment passes
#- 8.4.12: o-Anisidine (CP3021)
aps$CP3021
F1_S1_CP3021


#! TBD/check closely
#- 8.4.3: 2-Nitroaniline (CP2212) 
F6_S1_CP2212_revised <- remove_standard(aps$CP2212$F6_S1_CP2212, xl = 11.9, xu = 12.15)
#- 8.4.8: Menthone (CP3148)
F5_S1_CP3148_revised <- remove_standard(aps$CP3148$F5_S1_CP3148, xl = 8.1, xu = 8.3)

#!!! Multiple plots need to check
#- 8.4.13: Vernolate (CP2487)
for (plot_name in names(aps$CP2487)) {
  write_small(aps$CP2487[[plot_name]])
}

