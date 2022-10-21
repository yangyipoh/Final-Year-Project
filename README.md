# FYP 2022 - Drowsiness Detection and Anti-Sleep System for Drivers

The following repository containes the work completed by
1. James Gazzard
2. Alexander Schmedje
3. Yang Yi Poh

## Folder structure
```
Final-Year-Project
│   README.md    
│
└─── blinking-rate-detection
│   │   face_mesh.py
|   |   record.py
│   
└─── seatbelt-unit
|   │   ecg_adjust.m
|   |   main.m
|   |
|   └─── Arduino
|   |   |   firmware.ino
|   |   |   SparkFun_KX13X_Arduino_Library-main.zip
|   |
|   └─── functions
|   |   |   calib_rotation.m
|   |   |   edr1.m
|   |   |   findpeaks1.m
|   |   |   import_ctx_data.m
|   |   |   kernel_matrix.m
|   |   |   KPCA_EDR.m
|   |   |   kpca.m
|   |   |   Nonlinear_filter.m
|   |   |   preimage_rbf.m
|   |   |   RBF_kernel.m
|   |   |   Rpeaks_EDR.m
|   |   |   RSA_resp.m
|   |
|   └─── shell
|       |   PI3MK3M_shell.gcode
|       |   shell.SLDPRT
|       |   shell.STL
|
└─── steering-visualisation
|   |   plotting.py
|
└─── temperature-sensor\temperature_firmware
    |   temperature_firmware.ino
```
