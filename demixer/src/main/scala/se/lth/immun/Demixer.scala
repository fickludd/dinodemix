package se.lth.immun

import java.io.File
import java.io.FileReader
import java.io.FileWriter
import java.io.FileInputStream
import java.io.BufferedReader
import java.io.BufferedWriter
import java.io.InputStreamReader
import java.util.zip.GZIPInputStream
import java.util.Properties

import scala.collection.mutable.ArrayBuffer
import scala.collection.mutable.HashMap
import scala.io.Source

import se.jt.CLIApp

import se.lth.immun.xml.XmlReader
import se.lth.immun.xml.XmlWriter

import se.lth.immun.mzml._
import se.lth.immun.mzml.ghost.GhostSpectrum
import se.lth.immun.protocol.MsFeatureFile

import se.lth.immun.chem.Constants
import se.lth.immun.chem.IsotopeDistribution

import se.lth.immun.protocol.MsFeatures
import se.lth.immun.protocol.MsFeature

object Demixer extends CLIApp {
	
	val ISOLATION_WINDOW_TARGET = "MS:1000827"
	val ISOLATION_WINDOW_LOWER_OFF = "MS:1000828"
	val ISOLATION_WINDOW_UPPER_OFF = "MS:1000829"
	val SELECTED_ION_MZ_ACC = "MS:1000744"
	val CHARGE_STATE_ACC = "MS:1000041"
	val SCAN_START_TIME_ACC = "MS:1000016"
	val PEAK_INTENSITY_ACC = "MS:1000042"
	val MS_LEVEL_ACC = "MS:1000511"
				
	val C13C12_DIFF = 1.0033548378
	val ISOTOPE_PATTERN_DIFF = 1.00286864
		
	type PrecDef = Set[MzRange]
	case class MzRange(low:Double, high:Double, iw:IsolationWindow, a:Activation) {
		def holds(mz:Double) = mz >= low && mz < high
	}
	
	trait Ion {
		type Self <: Ion
		val self = this.asInstanceOf[Self]
		
		def z:Int
		def mz:Double
		def intensity:Double
		def within(ion:Ion, ppm:Double) = 
			2 * 1000000 * math.abs(mz - ion.mz) /	(mz + ion.mz) < ppm
		def withinOpt(ion:Ion, ppm:Double):Option[Self] =
			if (within(ion, ppm)) Some(self)
			else None
	}
	
	trait SpecType
	case object Ms1Spectrum extends SpecType
	case class OrigPrec(mz:Double, z:Int, intensity:Double, iw:IsolationWindow, a:Activation) extends SpecType with Ion {
		type Self = OrigPrec
	}
	case class IntFeature(intensity:Double, feature:MsFeature, iw:IsolationWindow, a:Activation) extends SpecType with Ion {
		type Self = IntFeature
		def mz = feature.mz
		def z = feature.z
	}
	
	case class PrecSuggestion(precMass:Double, precIntensity:Double, i1:Int, i2:Int)
	case class SpectrumSuggestion(precSuggestions:List[PrecSuggestion], iw:IsolationWindow, a:Activation) extends SpecType with Ion {
		type Self = SpectrumSuggestion
		def mass = precSuggestions.map(ps => ps.precMass * ps.precIntensity).sum / precSuggestions.map(_.precIntensity).sum
		def mz = mass / z + Constants.PROTON_WEIGHT
		def z = 2
		def intensity = precSuggestions.map(_.precIntensity).sum
	}
	
	case class SyncedPrec(mz:Double, intensity:Double, feat:Option[IntFeature], orig:Option[OrigPrec], comp:Option[SpectrumSuggestion]) {
		def specDef = feat.getOrElse(orig.getOrElse(comp.get))
	}
	
	case class OutSpectrumDef(precursor:SpecType, spectrum:Spectrum)
	
	case class DinoFeature(
		val mz:Double,
		val mostAbundantMz:Double,
		val z:Int,
		val rtStart:Double,
		val rtApex:Double,
		val rtEnd:Double,
		val fwhm:Double,
		val nIsotopes:Int,
		val nScans:Int,
		val averagineCorr:Double,
		val mass:Double,
		val massCalib:Double,
		val intensityApex:Double,
		val intensitySum:Double
	)
	
	val params = new DemixerParams
	
	var properties = new Properties
	properties.load(this.getClass.getResourceAsStream("/pom.properties"))
	val name 		= properties.getProperty("pom.artifactId")
	val version 	= properties.getProperty("pom.version")
	
	val t0 = System.currentTimeMillis
	val rt = Runtime.getRuntime
	
	def main(args:Array[String]):Unit = {
		
		failOnError(parseArgs(name, version, args, params, List("mzML"), None))
		val (outDir, outName) = params.outBase
		val outFile = new File(outDir, outName+".demix.mzML")
		
		
		println("  %s %s".format(name, version))
		println("      mzml: "+params.mzML.value)
		println("  features: "+params.msFeatures.value)
		println("  out file: "+outFile)
		
		//val features = readFeatures(new File(params.features)).sortBy(_.rtEnd)
		val features = 
			if (params.msFeatures.value != "")
				MsFeatureFile.read(new File(params.msFeatures.value), params.verbose)
			else MsFeatures(Nil, Nil)
		println("\n read %d features".format(features.features.length))
		
		
		val outSpectra = new ArrayBuffer[OutSpectrumDef]
		val r = getReader(params.mzML)
		
		if (params.verbose)
			println("\nindex   ms level   multiplex   time   heap(Mb)")
		
		val mzML = 
			MzML.fromFile(r, new MzMLDataHandlers(
				n => println(" reading %d spectra".format(n)),
				handleSpectrum(features, outSpectra),
				() => _,
				() => _
			))
			
		val w = new XmlWriter(new BufferedWriter(new FileWriter(outFile)))
	
		val nMs1 = outSpectra.map(_.precursor).collect { case Ms1Spectrum => 1 }.sum
		val nOrig = outSpectra.map(_.precursor).collect { case x:IntFeature => 1 }.sum
		val nDino = outSpectra.map(_.precursor).collect { case x:OrigPrec => 1 }.sum
		
		mzML.write(w, new MzMLDataWriters(
			outSpectra.length,
			writeSpectra(outSpectra),
			0,
			w => Unit
		))
	}
	
	
	
	var iStart = 0
	var iEnd = 0
	def handleSpectrum(
			msFeatures:MsFeatures, 
			spectra:ArrayBuffer[OutSpectrumDef]
	)(
			s:Spectrum
	):Unit = {
		
		val scanStartTimeOpt = s.scanList.flatMap(
				_.scans.head.cvParams.find(
					_.accession == SCAN_START_TIME_ACC
				).map(_.value.get.toDouble)
			)
			
		val msLevel = s.cvParams.find(_.accession == MS_LEVEL_ACC).map(_.value.get.toInt)
		
		val outSpectra = new ArrayBuffer[OutSpectrumDef]
			
		(msLevel, getPrecDef(s), scanStartTimeOpt) match {
			case (None, _, _) => 
				println("No Ms level information in spec num %d".format(s.index))
			
			case (Some(1), _, _) => {} // outSpectra += OutSpectrumDef(Ms1Spectrum, s)
			
			case (Some(2), None, _) =>
				println("Erroneous precursor definition at spec num %d".format(s.index))
			
			case (Some(2), _, None) =>
				println("No scanStartTime found for spec num %d".format(s.index))
			
			case (Some(2), Some(precDef), Some(scanStartTime)) =>
				iStart = msFeatures.features.indexWhere(_.rtStart < scanStartTime, iStart)
				iEnd = msFeatures.features.indexWhere(_.rtEnd < scanStartTime, iEnd)
				
				val origPcs = origPrecursors(s)
				
				val ms1Features = new ArrayBuffer[IntFeature]
				for {
					i <- iStart until iEnd
					mzRange <- precDef
				} {
					val f = msFeatures.features(i)
					val ionMass = featureMassOverlap(f, msFeatures.rtMap, scanStartTime, mzRange)
					if (ionMass > 0)
						ms1Features += IntFeature(ionMass, f, mzRange.iw, mzRange.a)
				}
				
				val t0 = System.currentTimeMillis
				val gs = GhostSpectrum.fromSpectrum(s)
				val complementaryFeatures = List[SpectrumSuggestion]() //spectralPrecursorGuesses(gs, precDef.head.iw, precDef.head.a) // pre-sorted
				
		
				if (params.verbose && msLevel.exists(_ == 2)) {
					println("\n ### SPECTRUM %6d n=%d".format(s.index, gs.mzs.length))
					println(gs.intensities.sorted.grouped(10).map(xs => "%.1e".format(xs.sum / xs.length)).mkString(" "))
				}
					//if (params.verbose)
				//	println(System.currentTimeMillis - t0 + " ms")
				
				val featurePcs = zip(ms1Features.sortBy(_.feature.mz), origPcs.sortBy(_.mz), complementaryFeatures, params.origPrecMzDiff)
				
				val (oldPcs, newPcs) = featurePcs.partition(_.orig.nonEmpty)
				
				for (syncedPrec <- oldPcs) {
					outSpectra += OutSpectrumDef(syncedPrec.specDef, s)
					if (params.verbose)
						println("%s  %s  %s  %.2e".format(
								if (syncedPrec.feat.nonEmpty) "FEAT" else "----",
								if (syncedPrec.orig.nonEmpty) "ORIG" else "----",
								if (syncedPrec.comp.nonEmpty) "COMP"+syncedPrec.comp.get.precSuggestions.length else "----",
								syncedPrec.intensity
							))
				}
					
				val descendingIntensity = newPcs.sortBy(_.intensity).reverse
					
				
					
				for (syncedPrec <- descendingIntensity.take(params.maxNewSpectra)) {
					outSpectra += OutSpectrumDef(syncedPrec.specDef, s)
					if (params.verbose)
						println("%s  ----  %s  new  %.2e".format(
								if (syncedPrec.feat.nonEmpty) "FEAT" else "----",
								if (syncedPrec.comp.nonEmpty) "COMP"+syncedPrec.comp.get.precSuggestions.length else "----",
								syncedPrec.intensity
							))
				}
				
					
				for (syncedPrec <- descendingIntensity.drop(params.maxNewSpectra).take(10)) {
					if (params.verbose)
						println("%s  ----  %s  dis  %.2e".format(
								if (syncedPrec.feat.nonEmpty) "FEAT" else "----",
								if (syncedPrec.comp.nonEmpty) "COMP"+syncedPrec.comp.get.precSuggestions.length else "----",
								syncedPrec.intensity
							))
				}
			
				
				
			case (Some(x), _, _) => 
				println("Cannot handle spectrum level %d in spec num %d".format(x, s.index))
		}
		
		if (params.verbose && s.index % params.readDebugFreq == 0 && false) 
			println("%10d %10d %10d %8d %8d Mb".format(
					s.index, 
					msLevel.getOrElse(0),
					outSpectra.length, 
					(System.currentTimeMillis - t0)/1000, 
					rt.totalMemory / 1000000
				))
					
		spectra ++= outSpectra
	}
	
	
	
	
	def featureMassOverlap(
			f:MsFeature, 
			rtMap:Seq[Double], 
			scanTime:Double, 
			mzRange:MzRange
	):Double = {
		val isotopeSizes = 
			for (h <- f.hills) yield {
				if (scanTime < rtMap(h.startIndex) || scanTime > rtMap(h.endIndex)) 0.0
				else {
					val i = (0 until h.intensity.length).minBy(i => math.abs(rtMap(h.startIndex + i) - scanTime))
					h.intensity(i)
				}
			}
		isotopeSizes.sum
	}
	
	
	
	def intensityOverlap(mzRange:MzRange, isoDist:IsotopeDistribution, feature:DinoFeature):Double = {
		var isoMass = 0.0
		for (i <- 0 until isoDist.intensities.length) {
			val isoMz = (feature.mass + i*ISOTOPE_PATTERN_DIFF) / feature.z + Constants.PROTON_WEIGHT
			if (mzRange.holds(isoMz))
				isoMass += isoDist.intensities(i)
		}
		isoMass * feature.intensitySum
	}
	
	
	
	def writeSpectra(spectra:Seq[OutSpectrumDef])(w:XmlWriter):Unit = 
		for (i <- 0 until spectra.length) { 
			val osd = spectra(i)
			osd.spectrum.index = i
			osd.precursor match {
				case Ms1Spectrum =>
					osd.spectrum.write(w)
				case df:IntFeature =>
					adjustByDino(df, osd.spectrum).write(w)
				case op:OrigPrec =>
					fixOrigPrec(op, osd.spectrum).write(w)
				case ss:SpectrumSuggestion =>
					
				
			}
		}
	
	
	def fixOrigPrec(op:OrigPrec, s:Spectrum):Spectrum = {
		val si = makeSelectedIon(op.mz, op.z, op.intensity)
		setPrecursor(si, op.iw, op.a, s)
	}
	
	
	
	def adjustByDino(df:IntFeature, s:Spectrum):Spectrum = {
		val si = makeSelectedIon(df.feature.mz, df.feature.z, df.intensity)
		setPrecursor(si, df.iw, df.a, s)
	}
	
	
	
	def adjustBySuggestion(ss:SpectrumSuggestion, s:Spectrum):Spectrum = {
		val si = makeSelectedIon(ss.mz, ss.z, ss.intensity)
		setPrecursor(si, ss.iw, ss.a, s)
	}
	
	
		
	def setPrecursor(
			si:SelectedIon, 
			iw:IsolationWindow,
			a:Activation,
			s:Spectrum
	):Spectrum = {
		
		val p = new Precursor
		p.isolationWindow = Some(iw)
		p.activation = a
		p.selectedIons += si
		
		val as = new Spectrum
		as.id = s.id
		as.defaultArrayLength = s.defaultArrayLength
		as.index = s.index
		as.dataProcessingRef = s.dataProcessingRef
		as.sourceFileRef = s.sourceFileRef
		as.spotID = s.spotID
		as.cvParams = s.cvParams
		as.userParams = s.userParams
		as.paramGroupRefs = s.paramGroupRefs
		as.scanList = s.scanList
		as.precursors += p
		as.products = s.products
		as.binaryDataArrays = s.binaryDataArrays
		
		as
	}
	
	
	
	def makeSelectedIon(mz:Double, z:Int, intensity:Double) = {
		val si = new SelectedIon 
		si.cvParams += new CvParam {
			accession = SELECTED_ION_MZ_ACC
			cvRef = "MS"
			name = "selected ion m/z"
			value = Some("%.8f".format(mz))
			unitCvRef = Some("MS")
			unitAccession = Some("MS:1000040")
			unitName = Some("m/z")
		}
		si.cvParams += new CvParam {
			accession = CHARGE_STATE_ACC
			cvRef = "MS"
			name = "charge state"
			value = Some(z.toString)
		}
		si.cvParams += new CvParam {
			accession = "MS:1000042"
			cvRef = "MS"
			name = "peak intensity"
			value = Some("%.2f".format(intensity))
			unitCvRef = Some("MS")
			unitAccession = Some("MS:1000131")
			unitName = Some("number of detector counts")
		}
		si
	}
	
	
	
	def zip(fs:Seq[IntFeature], origPcs:Seq[OrigPrec], complSpectra:Seq[SpectrumSuggestion], ppm:Double) = {
		var ifs = 0
		var iop = 0
		var icompl = 0
		
		def getOpt[X](xs:Seq[X], i:Int) =
			if (i < xs.length) Some(xs(i)) else None
		
		val res = new ArrayBuffer[SyncedPrec]
		while (ifs < fs.length || iop < origPcs.length || icompl < complSpectra.length) {
			val f 		= getOpt(fs, ifs)
			val op 		= getOpt(origPcs, iop)
			val compl 	= getOpt(complSpectra, icompl)
			
			val mz = Array(f, op, compl).flatten.minBy(_.mz)
			val fSyn = f.flatMap(_.withinOpt(mz, ppm))
			val opSyn = op.flatMap(_.withinOpt(mz, ppm))
			val complSyn = compl.flatMap(_.withinOpt(mz, ppm))
			
			if (fSyn.isDefined) ifs += 1
			if (opSyn.isDefined) iop += 1
			if (complSyn.isDefined) icompl += 1
			
			res += SyncedPrec(
					fSyn.getOrElse(opSyn.getOrElse(complSyn.get)).mz,
					fSyn.getOrElse(opSyn.getOrElse(complSyn.get)).intensity,
					fSyn,
					opSyn,
					complSyn
					)
		}
		res
	}
	
	
	
	def within(mz1:Double, mz2:Double, ppm:Double) = 
		2 * 1000000 * math.abs(mz1 - mz2) /	(mz1 + mz2) < ppm
	
	
		
	def getPrecDef(s:Spectrum):Option[PrecDef] = {
		val iws = 
			for {
				p <- s.precursors
				iw <- p.isolationWindow
			} yield isolationWindow2Range(iw, p.activation)
		if (iws.nonEmpty) Some(iws.toSet)
		else None
	}
	
		
		
	def isolationWindow2Range(iw:IsolationWindow, a:Activation):MzRange = {
		val iwTarget = iw.cvParams.find(_.accession == ISOLATION_WINDOW_TARGET)
		val iwLower = iw.cvParams.find(_.accession == ISOLATION_WINDOW_LOWER_OFF)
		val iwUpper = iw.cvParams.find(_.accession == ISOLATION_WINDOW_UPPER_OFF)
		
		(iwTarget, iwLower, iwUpper) match {
			case (Some(iwt), Some(iwl), Some(iwu)) =>
				val tmz = iwt.value.get.toDouble
				val lmz = iwl.value.get.toDouble
				val umz = iwu.value.get.toDouble
				MzRange(tmz-lmz, tmz+umz, iw, a)
			case _ =>
				throw new Exception("Erroneously defined isolation window: target=%s, lower=%s, upper=%s".format(iwTarget, iwLower, iwUpper))
				
		}
	}
	
	
	
	def origPrecursors(s:Spectrum):Seq[OrigPrec] = 
		for {
			p <- s.precursors
			si <- p.selectedIons
			iw <- p.isolationWindow
			(mz, z, int) <- selectedIon2Mz(si)
		} yield OrigPrec(mz, z, int, iw, p.activation)
	
	
	def selectedIon2Mz(si:SelectedIon):Option[(Double, Int, Double)] = 
		 for {
			 mzCV <- si.cvParams.find(_.accession == SELECTED_ION_MZ_ACC)
			 zCV <- si.cvParams.find(_.accession == CHARGE_STATE_ACC)
			 intCV <- si.cvParams.find(_.accession == PEAK_INTENSITY_ACC)
		 } yield (mzCV.value.get.toDouble, zCV.value.get.toInt, intCV.value.get.toDouble)
	
		
		
	
	val errs = new ArrayBuffer[String]
	def indexOfCol(header:Seq[String], col:String) = {
		val i = header.indexOf(col)
		if (i < 0)
			errs += "Could not find column '%s' in feature file".format(col)
		i
	}
		 
		 
	
	
	
	def getReader(path:String):XmlReader = {
		val f = new File(path)
		if (path.toLowerCase.endsWith(".mzml.gz"))
			new XmlReader(new BufferedReader(new InputStreamReader(
							new GZIPInputStream(new FileInputStream(f)))))
		else if (path.toLowerCase.endsWith(".mzml"))
			new XmlReader(new BufferedReader(new FileReader(f)))
		else
			throw new Exception("Unknown file format '%s'".format(path))
	}
	
	
	
	def readDinoFeatures(f:File):ArrayBuffer[DinoFeature] = {
		var headerParsed = false
		
		var iMZ = -1
		var iMZ_MOST_ABUNDANT = -1
		var iZ = -1
		var iRT_START = -1
		var iRT_APEX = -1
		var iRT_END = -1
		var iFWHM = -1
		var iN_ISO = -1
		var iN_SCAN = -1
		var iAVE_CORR = -1
		var iMASS = -1
		var iMASS_CAL = -1
		var iINT_APEX = -1
		var iINT_SUM = -1
		
		var res = new ArrayBuffer[DinoFeature]
		try {
			for (line <- Source.fromFile(f).getLines) {
				val cols = line.split("\t")map(_.trim)
				if (!headerParsed) {
					val header = cols
					iMZ = indexOfCol(header, "mz")
					iMZ_MOST_ABUNDANT = indexOfCol(header, "mostAbundantMz")
					iZ = indexOfCol(header, "charge")
					iRT_START = indexOfCol(header, "rtStart")
					iRT_APEX = indexOfCol(header, "rtApex")
					iRT_END = indexOfCol(header, "rtEnd")
					iFWHM = indexOfCol(header, "fwhm")
					iN_ISO = indexOfCol(header, "nIsotopes")
					iN_SCAN = indexOfCol(header, "nScans")
					iAVE_CORR = indexOfCol(header, "averagineCorr")
					iMASS = indexOfCol(header, "mass")
					iMASS_CAL = indexOfCol(header, "massCalib")
					iINT_APEX = indexOfCol(header, "intensityApex")
					iINT_SUM = indexOfCol(header, "intensitySum")
					headerParsed = true
				} else {
					
					res += DinoFeature(
						cols(iMZ).toDouble,
						cols(iMZ_MOST_ABUNDANT).toDouble,
						cols(iZ).toInt,
						cols(iRT_START).toDouble * 60,
						cols(iRT_APEX).toDouble * 60,
						cols(iRT_END).toDouble * 60,
						cols(iFWHM).toDouble,
						cols(iN_ISO).toInt,
						cols(iN_SCAN).toInt,
						cols(iAVE_CORR).toDouble,
						cols(iMASS).toDouble,
						cols(iMASS_CAL).toDouble,
						cols(iINT_APEX).toDouble,
						cols(iINT_SUM).toDouble
					)
				}
			}
		} catch {
			case e:Exception => {}
		}
		res
	}
	
	def spectralPrecursorGuesses(gs:GhostSpectrum, iw:IsolationWindow, a:Activation) = {
		
		def join(curr:SpectrumSuggestion, precs:List[PrecSuggestion], res:List[SpectrumSuggestion]):List[SpectrumSuggestion] = 
			precs match {
				case p::ps =>
					if (within(curr.mass, p.precMass, params.precursorGuessPPM))
						join(SpectrumSuggestion(p::curr.precSuggestions, iw, a), ps, res)
					else
						join(SpectrumSuggestion(List(p), iw, a), ps, curr::res)
				case Nil =>
					res.reverse
			}
		
		val precs = 
			for {
				i1 <- 0 until gs.mzs.length-1
				i2 <- i1+1 until gs.mzs.length
			} yield PrecSuggestion(
					gs.mzs(i1) + gs.mzs(i2) - 2*Constants.PROTON_WEIGHT, 
					gs.intensities(i1) + gs.intensities(i2), 
					i1, 
					i2
				)
		val sorted = precs.sortBy(_.precMass).toList
		
		join(SpectrumSuggestion(List(sorted.head), iw, a), sorted.tail, Nil).filter(_.precSuggestions.length > 2)
	}
}