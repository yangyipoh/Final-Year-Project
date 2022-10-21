# FYP 2022 - Drowsiness Detection and Anti-Sleep System for Drivers

The following repository containes the following code
1. Blinking rate detection code
2. Seatbelt code with EDR, firmware, and case
3. Steering wheel output visualisation
4. Temperature sensor code
5. Modified CARLA code

## Folder structure
```
Final-Year-Project
│   README.md    
│
└─── blinking-rate-detection
│   │   blinking.py
|
└─── CARLA-modified
|   |   manual_control_steeringwheel.py
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
## Requirements
To install all the packages needed for ```blinking-rate-detection``` and ```steering-visualisation``` (generated using pip freeze in Python 3.9.12), run the following commmand:
```
pip install -r requirements.txt
```
Note that CARLA will have its set of requirements to run as well. I would suggest setting up two different environments using Anaconda.

## Contact
Should anyone have any questions about the code, feel free to contact me at yangyipoh.pyy@gmail.com