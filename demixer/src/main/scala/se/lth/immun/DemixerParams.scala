package se.lth.immun

import se.jt.Params
import java.io.File

class DemixerParams extends Params {

	import Params._
	
	val mzML 			= ReqString("Input MzML with spectra")
	
	val msFeatures = ""				## "MsFeature file with features from the mzML"
	
	val origPrecMzDiff = 0.1 		## "include original precursor unless it's within diff to dino feature"
	val maxNewSpectra = 3			## "The maximum new precursored spectra to introduce (most intense precursors will be used)"
	val precursorGuessPPM = 5.0		## "Allowed PPM difference to merge different complementary precursor guesses"
	
	val featureMinZ = 2				## "Minimum feature charge state to consider for Demixing"
	val featureMaxZ = 100			## "Maximum feature charge state to consider for Demixing"
	val origMinZ = 2				## "Minimum original precursor charge state to keep"
	val origMaxZ = 100				## "Maximum original precursor charge state to keep"
	
	val readDebugFreq = 200			## "Output read status row every nth spectrum if verbose"
	val verbose = false				## "set to enable a lot of output"
	val specIntenseHistograms = false	## "set to enable spectrum intensity histogram output"
	val maxOutput = false			## "debug output mode"
	val discardOutput = false		## "set to discard output (for debugging)"
	
	val massCalib = ""				## "Path to mass calibration file"
	
	val includeMs1 = false			## "set to include MS1 spectra in output"
	val fakeIsolationWindow = false	## "set to fake isolation window according to demixed precursor mass"
	val indexedMzML = false			## "write indexed mzML"
	
	val gzipOutput 		= false		## "set to gzip out mzML"
	val outDir			= ""		## "output directory (by default same as input mzML)"
	val outName			= ""		## "basename for output files (by default same as input mzML)"
	val pipe			= ""		## "Full path were output should be written"
	
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
	
	def outFile = 
		if (pipe.value != "") new File(pipe)
		else {
			val (dir, name) = outBase
			new File(dir, name + ".demix.mzML" + (if (gzipOutput) ".gz" else ""))
		}
	
	def stripExt(path:String, ext:String) =
		if (path.toLowerCase.endsWith(ext.toLowerCase))
			path.dropRight(ext.length)
		else path
	
	def stripExts(path:String) =
		stripExt(stripExt(path, ".gz"), ".mzML")
		
}