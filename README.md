# ODS-Simulator
This repo contains the MATLAB programs written in support of the Master's thesis, "Development and Validation of Guidance, Navigation, and Control Strategies for CubeSat Earth Observation Missions" written by Elijah Vautour. The orbit determination system (ODS) simulator is a realistic simulation environment to design and test ODS algorithms. The ODS simulator was developed to design and test onboard ODS strategies. Rather than simulating the space environment itself, the ODS simulator reproduces the navigation conditions which an onboard ODS would be expected to operate in.

NOTE: This repo does not include the necessary datafiles from the GRACE-FO mission needed to run the OPG.m drive file. The programs are provided merely for viewing purposes. Please reach out to elijah.vautour@dal.ca for more details on the necessary datafiles and/or other questions about the programs.

File Descriptions:

1) ElijahVautour2022.pdf: This PDF file contains the Master's thesis, "Development and Validation of Guidance, Navigation, and Control Strategies for CubeSat Earth Observation Missions" written by Elijah Vautour. It is suggested that the viewer read chapters 3 and 4 (at the very least) to understand the software contained in this repo.

2) OPG.m: The OPG.m program is the main drive file used to test the different navigation algorithms contained in the NAV class.

3) GPS.m: The GPS class creates a simulated GPS sensor object which in turn generates GPS observables and single point solution (SPS) state estimates.

4) EKF.m: The EKF class contains the functions used in the reduced dynamic extended Kalman filter (RDEKF) algorithm and generates a configurable EKF filter object which stores user-defined filter properties.

5) NAV.m: The NAV class contains the actual navigation algorithms which are tested in the OPG.m simulator.

6) OPD.m: The OPD class contains the functions used to predict flybys with an Earth based targets viewing cone.

7) DATA.m: The DATA class contains the functions used for communication and data handling in the OPG.m simulator.

8) PARAM.m: The PARAM class contains the list of constant parameters used in the OPG.m simulator including natural constants and mathematical conversions.

9) PARSE.m: The PARSE class is used to parse the raw GRACE-FO datasets into MATLAB data structures for testing.
