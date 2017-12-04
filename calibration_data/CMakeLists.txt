# Install the calibration data so that programs can read it.
INSTALL (FILES
	el_Epratio_mom_2261.root el_Epratio_mom_4461.root
	FCP_2261_2250.dat FCP_4461_2250.dat
	EMCP_4461_2250.dat
	protdeltat_mom_2261_2250.root protdeltat_mom_4461_2250.root
	PFID_4461_2250.dat
	run_table.dat
	prot_mom_corr_18338_18438.dat
	vz_4461_4He.root
	vz_2261_56Fe.root
	vz_2261_4He.root
	vz_2261_3He.root
	vz_2261_12C.root
        DESTINATION "~/.e2a/")
