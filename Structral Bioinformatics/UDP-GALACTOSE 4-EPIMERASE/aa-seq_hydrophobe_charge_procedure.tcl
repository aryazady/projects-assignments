proc five {id} {
	set protein [atomselect $id "protein and name CA"]
	set primary [$protein get resname]

	set hydrophob [[atomselect 0 "protein and name CA and hydrophobic"] num]
	set positive [[atomselect $id "protein and name CA and charged and acidic"] num]
	set negative [[atomselect $id "protein and name CA and charged and basic"] num]
	set polar_neutral [[atomselect $id "protein and name CA and not charged and polar"] num]

	set allAA [llength $primary]

	set hydrophobPercentage [expr $hydrophob * 100.0 / $allAA]
	set positivePercentage [expr $positive * 100.0 / $allAA]
	set negativePercentage [expr $negative * 100.0 / $allAA]
	set polar_neutralPercentage [expr $polar_neutral * 100.0 / $allAA]

	set fp [open "5_answer.txt" w]
	puts $fp "primary = $primary\n\nhydrophob = $hydrophobPercentage\npositive = $positivePercentage\nnegative = $negativePercentage\npolar_neutral = $polar_neutralPercentage"
	close $fp
}