package se.lth.immun

import se.lth.immun.chem.Peptide

object AveragineSample {

	val MASS_LOW = 200
	val MASS_HIGH = 8000
	
	val AVERAGINE_ISO_DISTS = 
		for (m <- MASS_LOW to MASS_HIGH) yield {
			val elemComp = Peptide.averagine(m)
			elemComp.getIsotopeDistribution
		}
	
	def atMass(m:Double) = 
		AVERAGINE_ISO_DISTS(math.min(math.max(MASS_LOW, math.round(m).toInt), MASS_HIGH) - MASS_LOW)
}