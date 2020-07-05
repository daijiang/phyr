## code to prepare oldfield dataset

library(dplyr)

tree <- scan("http://www.plosone.org/article/fetchSingleRepresentation.action?uri=info:doi/10.1371/journal.pone.0007071.s002", 
             what=character(0)) %>% #download data
  ape::read.tree(text=.) ## read tree into ape

plot(tree, cex = 0.5)

comm_df <- readr::read_csv("http://www.plosone.org/article/fetchSingleRepresentation.action?uri=info:doi/10.1371/journal.pone.0007071.s001") 
comm_df

comm_long_df <- comm_df %>%
  tidyr::gather(sp, Abundance, -`Plot ID`, -`Habitat Type`) %>%
  dplyr::rename(site_orig = `Plot ID`) %>%
  dplyr::mutate(disturbance = ifelse(`Habitat Type` == "disturbed", 1, 0),
                site = paste0(site_orig, "_", `Habitat Type`),
                site_orig = as.integer(site_orig),
                Abundance = as.integer(Abundance),
                pres = ifelse(Abundance > 0, 1L, 0L)) %>% ## convert to presence / absence
  dplyr::rename(habitat_type = `Habitat Type`,
                abundance = Abundance)

comm_long_df

oldfield <- list(phy = tree,
                 data = comm_long_df)


usethis::use_data(oldfield, overwrite = TRUE)
