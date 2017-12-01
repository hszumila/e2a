# Install the calibration data so that programs can read it.
INSTALL (FILES el_Epratio_mom_4461_2250.root 
	FCP_4461_2250.dat EMCP_4461_2250.dat protdeltat_mom_4461_2250.root
	PFID_4461_2250.dat
	run_table.dat
	prot_mom_corr_18338_18438.dat
	vz_He_4461.root
        DESTINATION "~/.e2a/")
