# e2a
Repository for everything e2a. In this repository, you will find:

    - skim_tree.cpp : this is the flagship code of the repository
                         it skims data for a selected reaction type, applying PID and 
			 fiducial cuts
    - calibration_data : a directory with data files for applying cuts and corrections
    - maps : a directory containing acceptance maps, and the code used to create them
    - write_tree : a directory for code used to convert BOS files into a root format
                      suitable for skimming
    - Rey_neutrons : a directory for Rey's e2a work
    - Axel_analysis : a directory for Axel's e2a work

For a minimal working example of how one might analyze skimmed data,
    see the program Axel_analysis/example_analysis.cpp

Features we would like to add in the future:
	 - Acceptance maps
	 - Resolution maps
	 - Background subtraction codes