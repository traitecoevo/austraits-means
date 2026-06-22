# Assign each record to an Australian state/territory from its coordinates.
# A fast, approximate point-in-state rule (no spatial join): adequate for the
# coarse state grouping used in the species panels and summaries.
assign_state_fast <- function(data) {
  data |>
    dplyr::mutate(
      vic_nsw_border_lat = -34 - 0.25 * (decimal_longitude - 141),
      state = dplyr::case_when(
        is.na(decimal_longitude) | is.na(decimal_latitude) ~ NA_character_,
        decimal_latitude < -39 ~ "TAS",
        dplyr::between(decimal_longitude, 148.75, 149.40) & dplyr::between(decimal_latitude, -35.95, -35.10) ~ "ACT",
        decimal_longitude < 129 ~ "WA",
        decimal_longitude < 138 & decimal_latitude > -26 ~ "NT",
        decimal_longitude < 141 ~ "SA",
        decimal_latitude > -29 ~ "QLD",
        decimal_latitude < vic_nsw_border_lat ~ "VIC",
        TRUE ~ "NSW"
      )
    ) |>
    dplyr::select(-vic_nsw_border_lat) |>
    dplyr::relocate(state, .after = decimal_latitude)
}
