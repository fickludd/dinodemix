package se.lth.immun

import se.jt.Params
import java.io.File

class DinoMapperParams(val name:String, val version:String) extends Params {

	import Params._
	
	// USER EXPOSED PARAMS
	val verbose 		= false 	## "increase details in output"
	val matchPPM		= 10.0 		## "match threshhold in PPM"
	val matchPreRT		= 60.0		## "maximun allowed pre-emption of Ms2 compared to feature (sec)"
	val matchPostRT		= 60.0		## "maximun allowed delay of Ms2 compared to feature (sec)"
	val outDir			= ""		## "output directory (by default same as input mzML)"
	val outName			= ""		## "basename for output files (by default same as input mzML)"
	
	val dinoFeatures = ReqString("Csv-Merged dinosaur feature file")
	val searchResults = ReqString("search results from TPP (interact.pep.xls)")

	
	def outBase = {
		val mzMLFile = new File(dinoFeatures)
		val dir = 
			if (outDir.value != "") outDir.value
			else mzMLFile.getParent
		val name =
			if (outName.value != "") outName.value
			else stripExts(mzMLFile.getName)
		(dir, name) 
	}
	
	def stripExt(path:String, ext:String) =
		if (path.toLowerCase.endsWith(ext))
			path.dropRight(ext.length)
		else path
	
	def stripExts(path:String) =
		stripExt(stripExt(path, ".csv"), ".features")
	
			
	val outSep = "\t"
	val outQuote = "\""
}