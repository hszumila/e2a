# Install the calibration data so that programs can read it.
INSTALL (FILES
	el_Epratio_mom_2261.root el_Epratio_mom_4461.root
	el_Epratio_mom_4461_sim.root
	FCP_2261_2250.dat FCP_4461_2250.dat
	EMCP_2261_2250.dat EMCP_4461_2250.dat
	protdeltat_mom_2261_2250.root protdeltat_mom_4461_2250.root
	pipdeltat_mom_2261_2250.root pipdeltat_mom_4461_2250.root
	pimdeltat_mom_2261_2250.root pimdeltat_mom_4461_2250.root
	PFID_2261_2250.dat PFID_4461_2250.dat
	run_table.dat
	prot_mom_corr_18338_18438.dat
	vz_4461_4He.root vz_4461_12C.root
	vz_2261_4He.root vz_2261_12C.root vz_2261_3He.root vz_2261_56Fe.root 
	n_pathlength_corr_4461.dat
        DESTINATION "~/.e2a/")
