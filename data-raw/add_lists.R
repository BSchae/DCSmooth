dcs_list_options = c("type", "drv", "var_est", "qarma_order", "order_max",
                     "kerns", "infl_exp", "infl_par", "delta", "const_window",
                     "p_order")

dcs_list_kernels = c("MW200", "MW210", "MW220", "MW320", "MW420", "MW421", 
                     "MW422")

usethis::use_data(dcs_list_kernels, dcs_list_options,
                  overwrite = TRUE, internal = TRUE)
