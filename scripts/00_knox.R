source("R/common.R")
source("R/network_core.R")
source("R/knox_core.R")
cases <- validate_network_input(read_private_csv("cases_network.csv"))
xy <- as.matrix(cases[, c("utm_x", "utm_y")])
result <- knox_grid(xy, cases$onset_date)
safe_write_csv(result, "knox_all_cases.csv")
# The original analysis also considered local-only cases.
local <- cases$case_type == "local"
if (sum(local) >= 3L) safe_write_csv(knox_grid(xy[local, , drop = FALSE], cases$onset_date[local]),
                                   "knox_local_cases.csv")
message("Knox results written to results/private.")
