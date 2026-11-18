#cendr H^2 input
heritability_6reps <- total_regression_df %>%
  filter(! condition == "MO") %>%
  filter(! condition == "PV") %>%
  filter(! condition == "PA") %>%
  filter(! strain == "VC20019") %>%
  filter(! strain == "JU1400") %>%
  filter(! strain == "CB4856")

heritability_6reps_rep <- heritability_6reps %>%
  mutate(replicate_group = case_when(
    replicate %in% c("R01", "R02", "R05", "R07", "R09", "R11", "R13") ~ "1",
    replicate %in% c("R03", "R04", "R06", "R08", "R10", "R12", "R14") ~ "2",
    TRUE ~ NA_character_  # If there are other replicates not specified
  ))

AssayNumber <- c()
for (i in 1:nrow(heritability_6reps_rep)) {
  each_condition <- heritability_6reps_rep[i, 3]
  each_rep <- heritability_6reps_rep[i, 10]
  
  
  if (each_condition == "M114" & each_rep == "1") {
    AssayNumber <- append(AssayNumber, 1)
  } else if (each_condition == "M114" & each_rep == "2") {
    AssayNumber <- append(AssayNumber, 1)
  } else if (each_condition == "M117" & each_rep == "1") {
    AssayNumber <- append(AssayNumber, 2)
  } else if (each_condition == "M117" & each_rep == "2") {
    AssayNumber <- append(AssayNumber, 2)
  } else if (each_condition == "M395" & each_rep == "1") {
    AssayNumber <- append(AssayNumber, 3)
  } else if (each_condition == "M395" & each_rep == "2") {
    AssayNumber <- append(AssayNumber, 3)
  } else {
    print("error")
  }
}
heritability_6reps_rep$AssayNumber <- AssayNumber
heritability_6reps_input <- data.frame(AssayNumber= heritability_6reps_rep$AssayNumber,
                                       Strain = heritability_6reps_rep$strain,
                                       TraitName= rep("slope", nrow(heritability_6reps_rep)),
                                       Replicate = as.numeric(heritability_6reps_rep$replicate_group),
                                       Value = heritability_6reps_rep$z_score)
write_tsv(heritability_6reps_input, "regression_heritability_3assays_input.tsv")






