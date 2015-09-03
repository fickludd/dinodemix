package se.lth.immun

import se.jt.Params
import java.io.File

class DemixerParams extends Params {

	import Params._
	
	val mzML 			= ReqString("MzML with spectra to deconvolute")
	
	val msFeatures = ""				## "MsFeature file with features from the mzML"
	
	val origPrecMzDiff = 0.1 		## "include original precursor unless it's within diff to dino feature"
	val maxNewSpectra = 3			## "The maximum new precursored spectra to introduce (most intense precursors will be used)"
	val precursorGuessPPM = 5.0		## "Allowed PPM difference to merge different complementary precursor guesses"
	
	val featureMinZ = 2				## "Minimum feature charge state to consider for Demixing"
	val featureMaxZ = 100			## "Maximum feature charge state to consider for Demixing"
	
	val readDebugFreq = 200			## "Output read status row every nth spectrum if verbose"
	val verbose = false				## "set to enable a lot of output"
	
	val gzipOutput 		= false		## "set to gzip out mzML"
	val outDir			= ""		## "output directory (by default same as input mzML)"
	val outName			= ""		## "basename for output files (by default same as input mzML)"
	
	def outBase = {
		val mzMLFile = new File(mzML)
		val dir = 
			if (outDir.value != "") outDir.value
			else mzMLFile.getParent
		val name =
			if (outName.value != "") outName.value
			else stripExts(mzMLFile.getName)
		(dir, name) 
	}
	
	def outFile = {
		val (dir, name) = outBase
		new File(dir, name)
	}
	
	def stripExt(path:String, ext:String) =
		if (path.toLowerCase.endsWith(ext))
			path.dropRight(ext.length)
		else path
	
	def stripExts(path:String) =
		stripExt(stripExt(path, ".gz"), ".mzML")
		
}