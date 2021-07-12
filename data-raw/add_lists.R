dcs_list_options = c("type", "kerns", "p_order", "drv", "var_est",
                     "IPI_options")

dcs_list_kernels = c("MW_200", "MW_210", "MW_220", "MW_320", "MW_420",
                     "MW_421", "MW_422")

dcs_list_IPI = c("infl_exp", "infl_par", "delta", "const_window")

dcs_list_var_est = c("iid", "qarma", "qarma_gpac", "qarma_bic", "sarma", "lm")

usethis::use_data(dcs_list_kernels, dcs_list_options, dcs_list_IPI,
                  dcs_list_var_est,
                  overwrite = TRUE, internal = TRUE)
