dcs_list_options = c("type", "kerns", "p_order", "drv", "var_est",
                     "IPI_options")

dcs_list_kernels = c("MW_200", "MW_210", "MW_220", "MW_320", "MW_420",
                     "MW_421", "MW_422")

usethis::use_data(dcs_list_kernels, dcs_list_options,
                  overwrite = TRUE, internal = TRUE)
