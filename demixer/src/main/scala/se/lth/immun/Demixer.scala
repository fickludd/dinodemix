package se.lth.immun

import java.io.File
import java.io.FileReader
import java.io.FileWriter
import java.io.FileInputStream
import java.io.FileOutputStream
import java.io.BufferedReader
import java.io.BufferedWriter
import java.io.InputStreamReader
import java.io.OutputStreamWriter
import java.util.zip.GZIPInputStream
import java.util.zip.GZIPOutputStream
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
	
	val SECOND_ACC = "UO:0000010"
	val MINUTE_ACC = "UO:0000031"
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
	
	trait SpecType { def order:Double }
	case object Ms1Spectrum extends SpecType { val order = Double.MinValue }
	case class OrigPrec(mz:Double, z:Int, intensity:Double, iw:IsolationWindow, a:Activation) extends SpecType with Ion {
		type Self = OrigPrec
		def order = -intensity
	}
	case class IntFeature(intensity:Double, feature:MsFeature, iw:IsolationWindow, a:Activation) extends SpecType with Ion {
		type Self = IntFeature
		def mz = feature.mz
		def z = feature.z
		def order = -intensity
	}
	
	case class PrecSuggestion(precMass:Double, precIntensity:Double, i1:Int, i2:Int)
	case class SpectrumSuggestion(precSuggestions:List[PrecSuggestion], iw:IsolationWindow, a:Activation) extends SpecType with Ion {
		type Self = SpectrumSuggestion
		def mass = precSuggestions.map(ps => ps.precMass * ps.precIntensity).sum / precSuggestions.map(_.precIntensity).sum
		def mz = mass / z + Constants.PROTON_WEIGHT
		def z = 2
		def intensity = precSuggestions.map(_.precIntensity).sum
		def order = -intensity
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
		val outFile = params.outFile
		
		
		println("  %s %s".format(name, version))
		println("      mzml: "+params.mzML.value)
		println("  features: "+params.msFeatures.value)
		println("  out file: "+outFile)
		if (params.massCalib.value != "")
			println("mass calib: "+params.massCalib.value)
			
		
		//val features = readFeatures(new File(params.features)).sortBy(_.rtEnd)
		val msFeatures = 
			if (params.msFeatures.value != "")
				MsFeatureFile.read(new File(params.msFeatures.value), false) // params.verbose) <- too verbose
			else MsFeatures(Nil, Nil)
		
		val validFeatures = msFeatures.features.filter(f => 
				f.z >= params.featureMinZ && f.z <= params.featureMaxZ
			).sortBy(_.rtStart)
		
		val massCalib = 
			if (params.massCalib.value != "")
				parseMassCalib(params.massCalib.value)
			else
				0.0
				
		println(" using mass calib of %f ppm".format(massCalib))
			
		println("\n read %d features of which %d have z=%d-%d".format(
				msFeatures.features.length, 
				validFeatures.length, 
				params.featureMinZ.value, 
				params.featureMaxZ.value
			))
		
		
		val outSpectra = new ArrayBuffer[OutSpectrumDef]
		val r = getReader(params.mzML)
		
		val mzML = 
			MzML.fromFile(r, new MzMLDataHandlers(
				n => {
					println(" reading %d spectra".format(n))
					if (params.verbose)
						println("\nindex   ms level   multiplex   time     n datapoints  heap(Mb)%s"
								.format(if (params.specIntenseHistograms) "   spec.intens.histogram" else ""))
				},
				handleSpectrum(validFeatures, msFeatures.rtMap, outSpectra),
				() => _,
				() => _
			))
	
		val nMs1 = outSpectra.map(_.precursor).collect { case Ms1Spectrum => 1 }.sum
		val nOrig = outSpectra.map(_.precursor).collect { case x:OrigPrec => 1 }.sum
		val nFeature = outSpectra.map(_.precursor).collect { case x:IntFeature => 1 }.sum
		
		
		println("                n ms1: "+nMs1)
		println("    n orig precursors: "+nOrig)
		println(" n feature precursors: "+nFeature)
	
		val w = getWriter(outFile)
		
		mzML.write(w, new MzMLDataWriters(
			outSpectra.length,
			writeSpectra(outSpectra, massCalib),
			0,
			w => Nil
		), params.indexedMzML)
	}
	
	
	
	var iStart = 0
	var iEnd = 0
	def handleSpectrum(
			features:Seq[MsFeature], 
			rtMap:Seq[Double],
			spectra:ArrayBuffer[OutSpectrumDef]
	)(
			s:Spectrum
	):Unit = {
		
		val scanStartInSecondsOpt = getScanTime(s)
			
		val msLevel = s.cvParams.find(_.accession == MS_LEVEL_ACC).map(_.value.get.toInt)
		
		val outSpectra = new ArrayBuffer[OutSpectrumDef]
				
		fixIsolationWindowResolution(s)
			
		(msLevel, getPrecDef(s), scanStartInSecondsOpt) match {
			case (None, _, _) => 
				println("No Ms level information in spec num %d".format(s.index))
			
			case (Some(1), _, _) => 
				if (params.includeMs1)
					outSpectra += OutSpectrumDef(Ms1Spectrum, s)
			
			case (Some(2), None, _) =>
				println("Erroneous precursor definition at spec num %d".format(s.index))
			
			case (Some(2), _, None) =>
				println("No scanStartTime found for spec num %d".format(s.index))
			
			case (Some(2), Some(precDef), Some(scanTime)) =>
				
				def startLater(f:MsFeature) =
					f.rtEnd < scanTime
				
				while (iStart < features.length && startLater(features(iStart)))
					iStart += 1
				
				iEnd = math.max(iStart, iEnd)
				while (iEnd < features.length && features(iEnd).rtStart < scanTime)
					iEnd += 1
				
				val origPcs = origPrecursors(s).filter(op => op.z >= params.origMinZ && op.z <= params.origMaxZ)
				
				val ms1Features = new ArrayBuffer[IntFeature]
				for {
					i <- iStart until iEnd
					mzRange <- precDef
				} {
					val f = features(i)
					val ionMass = featureMassOverlap(f, rtMap, scanTime, mzRange)
					if (ionMass > 0)
						ms1Features += IntFeature(ionMass, f, 
							if (params.fakeIsolationWindow)
								iwFromFeature(f)
							else mzRange.iw, 
							mzRange.a)
				}
				
				val t0 = System.currentTimeMillis
				//val gs = GhostSpectrum.fromSpectrum(s) <- WARNING: REMOVES BINARY STRING FROM Spectrum
				val complementaryFeatures = List[SpectrumSuggestion]() //spectralPrecursorGuesses(gs, precDef.head.iw, precDef.head.a) // pre-sorted
				//if (params.verbose)
				//	println(System.currentTimeMillis - t0 + " ms")
				
				val featurePcs = zip(ms1Features.sortBy(_.feature.mz), origPcs.sortBy(_.mz), complementaryFeatures, params.origPrecMzDiff)
				
				val (oldPcs, newPcs) = featurePcs.partition(_.orig.nonEmpty)
				
				/*
				if (params.maxOutput && newPcs.nonEmpty) {
					println("\n ### SPECTRUM %6d n=%d".format(s.index, gs.mzs.length))
					if (complementaryFeatures.nonEmpty)
						println(gs.intensities.sorted.grouped(10).map(xs => "%.1e".format(xs.sum / xs.length)).mkString(" "))
				}
				*/
				
				for (syncedPrec <- oldPcs) {
					outSpectra += OutSpectrumDef(syncedPrec.specDef, s)
					if (params.maxOutput)
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
					if (params.maxOutput)
						println("%s  ----  %s  new  %.2e".format(
								if (syncedPrec.feat.nonEmpty) "FEAT" else "----",
								if (syncedPrec.comp.nonEmpty) "COMP"+syncedPrec.comp.get.precSuggestions.length else "----",
								syncedPrec.intensity
							))
				}
				
					
				for (syncedPrec <- descendingIntensity.drop(params.maxNewSpectra).take(10)) {
					if (params.maxOutput)
						println("%s  ----  %s  dis  %.2e".format(
								if (syncedPrec.feat.nonEmpty) "FEAT" else "----",
								if (syncedPrec.comp.nonEmpty) "COMP"+syncedPrec.comp.get.precSuggestions.length else "----",
								syncedPrec.intensity
							))
				}
			
				
				
			case (Some(x), _, _) => 
				println("Cannot handle spectrum level %d in spec num %d".format(x, s.index))
		}
		
		if (params.verbose && s.index % params.readDebugFreq == 0) {
			val binaryStrings = s.binaryDataArrays.map(_.binary)
			val gs = GhostSpectrum.fromSpectrum(s)
			println("%10d %10d %10d %8d %6d %8d Mb   %s".format(
					s.index, 
					msLevel.getOrElse(0),
					outSpectra.length, 
					(System.currentTimeMillis - t0)/1000,
					gs.mzs.length,
					rt.totalMemory / 1000000,
					if (params.specIntenseHistograms)
						gs.intensities.sorted.grouped(10).map(xs => "%.1e".format(xs.sum / xs.length)).mkString(" ")
					else ""
				))
			for ((bda, s) <- s.binaryDataArrays.zip(binaryStrings))
				bda.binary = s
		}
					
		if (!params.discardOutput)
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
				if (
					scanTime < rtMap(h.startIndex) || scanTime > rtMap(h.endIndex) ||
					h.mz < mzRange.low - h.mzErr || h.mz > mzRange.high + h.mzErr
						) 0.0
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
	
	
	
	def writeSpectra(
			spectra:Seq[OutSpectrumDef], 
			massCalib:Double
	)(
			w:XmlWriter
	):Seq[MzML.OffsetRef] = {
		val t0 = System.currentTimeMillis
		if (params.verbose)
			println(" index      spec type  time")
		
		var i = 0
		val byOrigAndOrder = 
			spectra.groupBy(_.spectrum.index).values
				.map(_.sortBy(_.precursor.order)).toArray
		for {
			specGroup <- byOrigAndOrder
			j <- 0 until specGroup.length
		} yield {
			val osd = specGroup(j)
			osd.spectrum.index = i
			i += 1
			
			if (params.verbose && i % params.readDebugFreq == 0)
				println("%10d %10s %10ds".format(
						i,
						specTypeToString(osd.precursor),
						(System.currentTimeMillis - t0)/1000			
					))
			
			val s = 
				osd.precursor match {
					case Ms1Spectrum =>
						osd.spectrum
					case df:IntFeature =>
						adjustByDino(df, massCalib, osd.spectrum)
					case op:OrigPrec =>
						fixOrigPrec(op, massCalib, osd.spectrum)
					case ss:SpectrumSuggestion =>
						adjustBySuggestion(ss, massCalib, osd.spectrum)
				}
			
			val scanNum = parseScanNum(osd.spectrum.id)
			s.id = s.id + " scan=%d precRank=%d".format(scanNum, j)
			
			s.write(w, null, None, None)
		}
	}
	
	
	def specTypeToString(specType:SpecType) = 
		specType match {
			case Ms1Spectrum 	=> "MS1"
			case df:IntFeature 	=> "FEAT"
			case op:OrigPrec 	=> "ORIG"
			case ss:SpectrumSuggestion => "COMPL"
		}
	
	
	
	def fixOrigPrec(op:OrigPrec, massCalib:Double, s:Spectrum):Spectrum = {
		val si = makeSelectedIon(op.mz - massCalib * op.mz / 1e6, op.z, op.intensity)
		val ms2rt = getScanTime(s).get
		val id = "ORIG precInt=%.2e rtInSec=%.1f".format(op.intensity, ms2rt)
		setPrecursor(si, op.iw, op.a, s, id)
	}
	
	
	
	def adjustByDino(df:IntFeature, massCalib:Double, s:Spectrum):Spectrum = {
		val si = makeSelectedIon(df.feature.mz - massCalib * df.feature.mz / 1e6, df.feature.z, df.intensity)
		val id = "FEAT precInt=%.2e rtInSec=%.1f featApexInt=%.1e".format(df.intensity, df.feature.rtApex, df.feature.intensityApex)
		setPrecursor(si, df.iw, df.a, s, id)
	}
	
	
	
	def adjustBySuggestion(ss:SpectrumSuggestion, massCalib:Double, s:Spectrum):Spectrum = {
		val si = makeSelectedIon(ss.mz - massCalib * ss.mz / 1e6, ss.z, ss.intensity)
		val id = "COMPL_FRAG"
		setPrecursor(si, ss.iw, ss.a, s, id)
	}
	
	
	
	val scanNumRE = """scan=(\d+)""".r.unanchored
	def parseScanNum(id:String) = {
		id match {
			case scanNumRE(num) => num.toInt
			case _ => 
				throw new Exception("Couldn't parse scan number from spectrum title '%s'".format(id))
		}
	}
	
	
		
	def setPrecursor(
			si:SelectedIon, 
			iw:IsolationWindow,
			a:Activation,
			s:Spectrum,
			id:String
	):Spectrum = {
		
		val p = new Precursor
		p.isolationWindow = Some(iw)
		p.activation = a
		p.selectedIons += si
		
		val as = new Spectrum
		as.id = id
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
	
	
	
	def iwFromFeature(f:MsFeature):IsolationWindow = {
		val iw = new IsolationWindow
		iw.cvParams += new CvParam {
			accession = ISOLATION_WINDOW_TARGET
			cvRef = "MS"
			name = "isolation window target m/z"
			value = Some("%.8f".format(f.mz))
			unitCvRef = Some("MS")
			unitAccession = Some("MS:1000040")
			unitName = Some("m/z")
		}
		/*
		iw.cvParams += new CvParam {
			accession = ISOLATION_WINDOW_LOWER_OFF
			cvRef = "MS"
			name = "isolation window lower offset"
			value = Some("2.0")
			unitCvRef = Some("MS")
			unitAccession = Some("MS:1000040")
			unitName = Some("m/z")
		}
		iw.cvParams += new CvParam {
			accession = ISOLATION_WINDOW_UPPER_OFF
			cvRef = "MS"
			name = "isolation window upper offset"
			value = Some("2.0")
			unitCvRef = Some("MS")
			unitAccession = Some("MS:1000040")
			unitName = Some("m/z")
		}
		*/
		iw
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
	
	
	def getScanTime(s:Spectrum):Option[Double] = 
		s.scanList.flatMap(
				_.scans.head.cvParams.find(
					_.accession == SCAN_START_TIME_ACC
				).map(cv => {
					cv.unitAccession match {
						case Some(SECOND_ACC) =>
							cv.value.get.toDouble
						case Some(MINUTE_ACC) =>
							cv.value.get.toDouble * 60.0
						case _ =>
							error(s, "Unknown unit of scan time!")
					}
				})
			)
	
	
			
	def error(s:Spectrum, msg:String) =
		throw new Exception("[MZML ERR SCAN=%d] %s".format(s.index, msg))
			
			
	
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
	
	
	
	def fixIsolationWindowResolution(s:Spectrum) = {
		for {
			p 			<- s.precursors
			si 			<- p.selectedIons
			iw 			<- p.isolationWindow
			mzCV 		<- si.cvParams.find(_.accession == SELECTED_ION_MZ_ACC)
			iwTargetCV 	<- iw.cvParams.find(_.accession == ISOLATION_WINDOW_TARGET)
		} {
			if (within(mzCV.value.get.toDouble, iwTargetCV.value.get.toDouble, 1000))
				iwTargetCV.value = Some(mzCV.value.get)
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
	
	
	
	def getWriter(f:File):XmlWriter = {
		val name = f.getName
		if (name.toLowerCase.endsWith(".mzml.gz"))
			XmlWriter(f, true)
		else if (name.toLowerCase.endsWith(".mzml"))
			XmlWriter(f, false)
		else
			throw new Exception("Unknown file format '%s'".format(name))
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
	
	
	def parseMassCalib(path:String):Double = {
		val r = new BufferedReader(new FileReader(new File(path)))
		val massCalibHeader = r.readLine
		r.close
		
		val meanRE = """mean mz diff=([-\d.]+)ppm""".r.unanchored
		massCalibHeader match {
			case meanRE(diffPPM) =>
				diffPPM.toDouble
			case _ => 0.0
		}
	}
	
	
	// REMEMBER to fake IsolationWindow depending on params!
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