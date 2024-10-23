proc six {id} {
    set sulfids [atomselect $id "protein and resname CYS and element S"]
    set sulfidsID [$sulfids get index]
    set sulfidsCount [$sulfids num]

    set f [open "6_answer.txt" w]
    for {set atom1 0} {$atom1 < $sulfidsCount} {incr atom1} {
        for {set atom2 [expr $atom1 + 1]} {$atom2 < $sulfidsCount} {incr atom2} {
            set dist [measure bond [list [lindex $sulfidsID $atom1] [lindex $sulfidsID $atom2]]]
            if {$dist < 3} {
                set resid1 [atomselect $id "index [lindex $sulfidsID $atom1]"]
                set resid2 [atomselect $id "index [lindex $sulfidsID $atom2]"]
                puts $f "Disulfid Bonds between resid: [$resid1 get resid] & [$resid2 get resid]"
            }
        }
    }
    close $f
}