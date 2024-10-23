proc secStruct {id} {
	set frames [molinfo $id get numframes]
	for {set f 0}  {$f < $frames} {incr f} {
		set protein [atomselect $id "protein and name CA" frame $f]
		set secondry [$protein get structure]
		set c 0.0
		set e 0.0
		set h 0.0
		set t 0.0
		set len [llength $secondry]
		foreach s $secondry {
			if {$s == "C"} {
				set c [expr $c + 1]
			} elseif {$s == "E"} {
				set e [expr $e + 1]
			} elseif {$s == "H"} {
				set h [expr $h + 1]
			} elseif {$s == "T"} {
				set t [expr $t + 1]
			}
		}
		
		set c [format "%.2f" [expr $c / $len * 100]]
		set e [format "%.2f" [expr $e / $len * 100]]
		set h [format "%.2f" [expr $h / $len * 100]]
		set t [format "%.2f" [expr $t / $len * 100]]
		
		set fp [open "secondry_structure_percentage.txt" a]
		puts $fp "Frame $f"
		puts $fp "secondry strcuture : $secondry"
		puts $fp "C: $c %\nE: $e %\nH: $h %\nT: $t %\n"
		close $fp
	}
}