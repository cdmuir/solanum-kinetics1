write_rds(estimates_by_model2, "objects/estimates_by_model.rds")

ggplot(
  filter(estimates_by_model2, type == "b"),
  aes(
    estimate,
    model,
    xmin = `q2.5`,
    xmax = `q97.5`,
    color = overlap_zero
  )
) +
  facet_grid(explanatory ~ response, scales = "free_x") +
  geom_vline(xintercept = 0,
             color = "grey",
             linetype = "dashed") +
  geom_pointinterval()

ggplot(
  filter(estimates_by_model2, type == "cor_phy"),
  aes(
    estimate,
    model,
    xmin = `q2.5`,
    xmax = `q97.5`,
    color = overlap_zero
  )
) +
  facet_grid(explanatory ~ response, scales = "free_x") +
  geom_vline(xintercept = 0,
             color = "grey",
             linetype = "dashed") +
  geom_pointinterval()

selected_model = ?
  
  tbl_comparison |>
  mutate(selected = model == selected_model) |>
  write_rds("objects/tbl-comparison.rds")


write_rds(fits$fit[[as.numeric(str_remove(selected_model, "model"))]], "objects/best_model.rds")




# Write "best" model
# Note: I reran the models several times and the order of the top models
# changed, consistent with the difference in LOOIC being caused by sampling
# variability. After reviewing model estimates, I determined that model 6 is the
# clearest to interpret.

# Previous version (lowest LOOIC)
# write_rds(fits$fit[[as.numeric(str_extract(rownames(looic_table)[1], "\\d+"))]], "objects/best_model.rds")

# Current version (model 6)
write_rds(fits$fit[[as.numeric(str_remove(selected_model, "model"))]], "objects/best_model.rds")


