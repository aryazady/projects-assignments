proc aminoAcid {id} {
	set protein [atomselect $id "protein and name CA"]
	set primary [$protein get resname]
	
	set counters {}
	foreach item $primary {
		dict incr counters $item
	}
	
	set fp [open "primary_structre.txt" w]
	puts $fp "Primary Structure: $primary"
	puts $fp "Each has:"
	dict for {item count} $counters {
		puts $fp "${item}: $count"
	}
	close $fp
}