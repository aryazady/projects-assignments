proc phiPsi {id} {
	set protein [atomselect $id "protein and name CA" frame 0]
	set phi_psi [$protein get {phi psi}]
	
	set fp [open "phi_psi_angels_for_CA.txt" w]
	foreach pp $phi_psi {
		puts $fp $pp
	}
	close $fp
}