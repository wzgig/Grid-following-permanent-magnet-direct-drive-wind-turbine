# Grid-following-permanent-magnet-direct-drive-wind-turbine
跟网型直驱风机的解析模型，采用稳态运行工作点进行建模仿真
copy中全为电动机惯例，wangce.m中机侧为电动机，网侧为发电机。但是经过更改，两种算出来的结果是一致的，包括电流正负号 Added Copy_of_wangce.m, Copy_of_wangceRef_and_Para.m, 和 Copy_of_wangce_Main.m as new files. Updated wangce.m, wangceRef_and_Para.m, 和 wangce_Main.m to use P_dc instead of P_g for DC power calculation and added P_dc to the output plots. This improves the accuracy of the DC link power computation and aligns the main and reference files.
