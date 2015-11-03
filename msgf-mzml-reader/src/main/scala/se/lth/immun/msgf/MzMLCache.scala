package se.lth.immun.msgf

import java.io.File
import java.io.FileReader
import java.io.FileInputStream
import java.io.BufferedReader
import java.io.InputStreamReader
import java.util.zip.GZIPInputStream
import se.lth.immun.xml.XmlReader
import se.lth.immun.mzml.MzML
import se.lth.immun.mzml.MzMLDataHandlers
import se.lth.immun.mzml.ghost.GhostSpectrum

import scala.collection.mutable.ArrayBuffer

import se.lth.immun.mzml.CvParam

object MzMLCache {
	
	val MS_LEVEL_ACC = "MS:1000511"
	val BASE_PEAK_INT_ACC = "MS:1000505"
	val BASE_PEAK_MZ_ACC = "MS:1000504"
	val TIC_ACC = "MS:1000285"
	val SCAN_START_TIME_ACC = "MS:1000016"
	val SCAN_WINDOW_MIN_ACC = "MS:1000501"
	val SCAN_WINDOW_MAX_ACC = "MS:1000500"
	val ISOLATION_WINDOW_TARGET = "MS:1000827"
	val ISOLATION_WINDOW_LOWER_OFF = "MS:1000828"
	val ISOLATION_WINDOW_UPPER_OFF = "MS:1000829"
	val SELECTED_ION_MZ_ACC = "MS:1000744"
	val CHARGE_STATE_ACC = "MS:1000041"
	val PEAK_INTENSITY_ACC = "MS:1000042"
	val CENTROID_ACC = "MS:1000127"
	val PROFILE_ACC = "MS:1000128"
	val THERMO_MONO_MZ_NAME = "[Thermo Trailer Extra]Monoisotopic M/Z:"
		
	val CID_ACC = "MS:1000133"
	val ETD_ACC = "MS:1000598"
	val HCD_ACC = "MS:1000422"
	val PQD_ACC = "MS:1000599"
	val KNOWN_ACTIVATION = Array(CID_ACC, ETD_ACC, HCD_ACC, PQD_ACC)
	
	case class SpecPrec(
			precMz:Double, 
			precZ:Int, 
			precInt:Double, 
			isolationWindowTargetMz:Double,
			activationAcc:Array[String]
		)
	
	case class SpecData(
			id:String,
			scanNum:Int,
			isCentroided:Boolean,
			msLevel:Int,
			precs:Array[SpecPrec],
			mzs:Array[Double],
			intensities:Array[Double]
		)
}

class MzMLCache(val specFile:File, msLevelMin:Int = 2, msLevelMax:Int = 2, intensityThreshold:Double = 0.0) {

	import MzMLCache._
	
	val spectra = new ArrayBuffer[SpecData]
	val mzML = read(specFile, spectra)
	
	
	def read(f:File, spectra:ArrayBuffer[SpecData]):MzML = {
		val r = getReader(specFile)
		val dh = new MzMLDataHandlers(
				n => {},
				s => {
					val gs = GhostSpectrum.fromSpectrum(s)
					if (gs.msLevel >= msLevelMin && gs.msLevel <= msLevelMax)
						spectra += toSpecData(gs)
				},
				n => {},
				c => {}
			)
		MzML.fromFile(r, dh, null)
	}
	
	
	def toSpecData(gs:GhostSpectrum):SpecData = {
		
		def cvDouble(cvParams:Seq[CvParam], acc:String) =
			cvParams.find(_.accession == acc).map(_.value.get.toDouble)
		def cvInt(cvParams:Seq[CvParam], acc:String) =
			cvParams.find(_.accession == acc).map(_.value.get.toInt)
		
		val isCentroided = gs.spectrum.cvParams.exists(_.accession == CENTROID_ACC)
		val thermoPrecMonoMz = 
			for {
				sl <- gs.spectrum.scanList
				scan <- sl.scans.headOption
				up <- scan.userParams.find(_.name == THERMO_MONO_MZ_NAME)
			} yield up.value.get.toDouble
			
		val specPrecs =
			(
				for (pc <- gs.spectrum.precursors) yield {
					val isoWinTarget = 
						for {
							iw <- pc.isolationWindow
							x <- cvDouble(iw.cvParams, ISOLATION_WINDOW_TARGET)
						} yield x
					
					val activations = 
						pc.activation.cvParams.map(_.accession).filter(KNOWN_ACTIVATION.contains)
						
					for {
						si <- pc.selectedIons
						precMz <- cvDouble(si.cvParams, SELECTED_ION_MZ_ACC)
					} yield {
						val precZ 			= cvInt(si.cvParams, CHARGE_STATE_ACC).getOrElse(2)
						val precIntensity 	= cvDouble(si.cvParams, PEAK_INTENSITY_ACC).getOrElse(1000.0)
						SpecPrec(precMz, precZ, precIntensity, isoWinTarget.getOrElse(precMz), activations.toArray)
					} 
				}).flatten
			
		
		val intensityFiltered = gs.mzs.zip(gs.intensities).filter(_._2 >= intensityThreshold)
		SpecData(
				gs.id,
				gs.spectrum.index,
				isCentroided,
				gs.msLevel,
				specPrecs.toArray,
				intensityFiltered.map(_._1).toArray,
				intensityFiltered.map(_._2).toArray
			)
			
	}
	
	
	def getReader(f:File):XmlReader = {
		val path = f.toString
		if (path.toLowerCase.endsWith(".mzml.gz"))
			new XmlReader(new BufferedReader(new InputStreamReader(
							new GZIPInputStream(new FileInputStream(f)))))
		else if (path.toLowerCase.endsWith(".mzml"))
			new XmlReader(new BufferedReader(new FileReader(f)))
		else
			throw new Exception("Unknown file format '%s'".format(path))
	}
	
	
	def spectrumIDFormatCvParam:CvParam =
		mzML.fileDescription.sourceFiles.flatMap(_.cvParams).find(cv => {
			val accNum = cv.accession.split(":").last.toLong
			(accNum >= 1000768 && accNum <= 1000777 
				|| accNum == 1000823 || accNum == 1000824 || accNum == 1000929
				|| accNum == 1001508 || accNum == 1001526 || accNum == 1001528
				|| accNum == 1001531 || accNum == 1001532
				|| accNum == 1001559 || accNum == 1001562
		)}).getOrElse(null)
	
}