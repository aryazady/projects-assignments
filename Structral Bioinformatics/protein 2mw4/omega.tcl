proc findHydrogen {molID} {
    set protein [atomselect $molID "protein"]
    set angel 20
    set cutoff 3
    set hBonds [measure hbonds $cutoff $angel $protein]

    set data ""
    foreach accDonHyd $hBonds {
        foreach atom $accDonHyd {
            lappend data [[atomselect $molID "protein and index $atom"] get resid]
        }
    }

    set f [open "amino_acid_psi_phi_omega.txt" w]
    puts $f "Resid\tResname\tStructure\tPhi\t\t\t\t\tPsi\t\t\t\t\tOmega\t\t\t\tRadius of Gyration"
    foreach resid [lsort -unique $data] {
        set fourAtom [atomselect $molID "protein and backbone and resid $resid to [expr $resid + 1]"]
        set CA1 [lindex [$fourAtom get index] 1]
        set C [lindex [$fourAtom get index] 2]
        set N [lindex [$fourAtom get index] 0]
        set CA2 [lindex [$fourAtom get index] 5]

        set omega [measure dihed "$CA1 $C $N $CA2"]
        set rGyration [measure rgyr [atomselect $molID "resid $resid"]]
        set resname [[atomselect $molID "name CA and resid $resid"] get resname]
        set structure [[atomselect $molID "name CA and resid $resid"] get structure]
        set phi [[atomselect $molID "name CA and resid $resid"] get phi]
        set psi [[atomselect $molID "name CA and resid $resid"] get psi]

        puts $f "$resid\t\t$resname\t\t$structure\t\t\t$phi\t$psi\t$omega\t$rGyration"
    }
    close $f
    return
}