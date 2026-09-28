explorer_metadata <- read.delim(
  "data-raw/explorer.tsv",
  encoding = "latin1",
  na.strings = ""
) |>
  subset(select = c("Name", "Formula", "CAS", "HMDB", "LMP")) |>
  subset(!is.na(Formula)) |>
  tail(30)

usethis::use_data(explorer_metadata, overwrite = TRUE)
