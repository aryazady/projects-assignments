proc rmsd {f} {
	set id [mol load pdb $f]
	set frames [molinfo $id get numframes]
	
	set fp [open "rmsd.txt" w]
	puts $fp "RMSD between Frames"
	for {set i 0} {$i < $frames} {incr i} {
		for {set j [expr $i + 1]} {$j < $frames} {incr j} {
			set value [measure rmsd [atomselect $id "protein" frame $j] [atomselect $id "protein" frame $i]]
			puts $fp "($i, $j) : $value"
		}
	}
	close $fp
}