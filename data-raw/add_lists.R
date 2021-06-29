dcs_list_options = c("type", "drv", "var_est", "qarma_order", "order_max",
                     "kerns", "infl_exp", "infl_par", "delta", "const_window",
                     "p_order")

dcs_list_kernels = c("MW_200", "MW_210", "MW_220", "MW_320", "MW_420",
                     "MW_421", "MW_422")

usethis::use_data(dcs_list_kernels, dcs_list_options,
                  overwrite = TRUE, internal = TRUE)
