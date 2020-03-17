# competency
Analysis of inoculation tests for determining competencies of hosts for P. ramorum.

The most common tree species in the Big Sur area were used to test the competency of P. ramorum, the causative agent of Sudden Oak Death. We inoculated leaves from 32 individuals for each species with inoculum for the treatment and water for the controls. We quantified spore production of sporangia, the spores important for transmission, and chlamydospores, the spores potentially important for survival in the soil. 

Models used for manuscript analysis located in rscripts/analysis2019 folder. I ran 4 seperate poisson glmms using brms:
M1: sporangia ~ species
M2: chlamydospores ~ species
M3: sporangia ~ lesion_size + (lesion_size|species)
m4: chlamydospores ~ lesion_size + (lesion_size|species)

Our results showed that once challenged by P. ramorum, almost all plant species were infected and produced spores to some extent. We found sporangia production was greatest in tanoak and bay laurel and the difference was insignificant, and even though other species produced much less, quantities were non-zero. Thus, additional species may play a previously unrecognized role in local transmission.
