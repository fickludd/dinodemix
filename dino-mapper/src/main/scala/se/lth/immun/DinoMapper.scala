package se.lth.immun

import se.jt.CLIApp
import java.util.Properties
import java.io.File
import java.io.FileWriter
import java.io.BufferedWriter

import scala.io.Source
import scala.collection.mutable.ArrayBuffer
import scala.collection.mutable.HashMap

object DinoMapper extends CLIApp {

	
	case class DinoFeature(
		val mz:Double,
		val z:Int,
		val rtStart:Double,
		val rtApex:Double,
		val rtEnd:Double,
		val nIsotopes:Int,
		val nScans:Int,
		val averagineCorr:Double,
		val mass:Double,
		val massCalib:Double,
		val intensityApex:Double,
		val intensitySum:Double
	)
	
	case class TppPeptide(
		val probability:Double,
		val scanIndex:Int,
		val peptide:Int,
		val protein:Int,
		val theoMass:Double,
		val empMz:Double,
		val empMass:Double,
		val z:Int,
		val rt:Double
	)
	
	
	var params:DinoMapperParams = _
	val errs = new ArrayBuffer[String]
	
	def indexOfCol(header:Seq[String], col:String) = {
		val i = header.indexOf(col)
		if (i < 0)
			errs += "Could not find column '%s' in feature file".format(col)
		i
	}
	
	
	def main(args:Array[String]):Unit = {
		
		var properties = new Properties
    	properties.load(this.getClass.getResourceAsStream("/pom.properties"))
    	val name 		= properties.getProperty("pom.artifactId")
    	val version 	= properties.getProperty("pom.version")
    	
    	params = new DinoMapperParams(name, version)
    	
		val t0 	= System.currentTimeMillis
		
		failOnError(parseArgs(name, version, args, params, List("dinoFeatures", "searchResults"), None))
		
		val (outDir, outName) = params.outBase
		
		println(name + " "+version)
    	println("   dino feature file: " + params.dinoFeatures.value)
    	println("  search result file: " + params.searchResults.value)
    	println("             out dir: " + outDir)
    	println("            out name: " + outName)
		println()
		
		println("reading input files...")
		
		val features = readFeatures(new File(params.dinoFeatures))
		println("features from %d files read".format(features.size))
		println("%25s %8s".format("file", "n feats"))
		for ((fileName, feats) <- features)
			println("%25s %8d".format(fileName, feats.length))
		println
		
		val (searchResults, peptides, proteins) = readSearchResults(new File(params.searchResults))
		println("search results from %d files read".format(searchResults.size))
		println("%25s %8s".format("file", "n peps"))
		for ((fileName, peps) <- searchResults)
			println("%25s %8d".format(fileName, peps.length))
		println
		
		if (errs.nonEmpty) {
			for (err <- errs)
				println(err)
			System.exit(1)
		}
		
		val t1 	= System.currentTimeMillis
		println("matching features to peptides...")
		
		val matches = 
			for (fileName <- features.keys.toSeq) 
				yield (fileName, matchFeatures(features(fileName), searchResults(fileName)))
				
		println("found %d matches".format(matches.map(_._2.map(_._2.count(_.isDefined)).sum).sum))
		val t2 	= System.currentTimeMillis
		println("writing matches...")
		
		val outFile = new File(outDir, outName+".matches.csv")
		writeMatches(matches, outFile, peptides, proteins)
		
		val t3 	= System.currentTimeMillis
		println("total time: "+niceTiming(t3-t0))
	}
	
	
	
	
	def writeMatches(
			matches:Seq[(String, Map[Int, IndexedSeq[Option[(TppPeptide, DinoFeature)]]])], 
			f:File,
			peptides:Array[String],
			proteins:Array[String]
	) = {
		val w = new BufferedWriter(new FileWriter(f))
		
		def writeRow(qoute:Boolean)(a:Any*) = 
			w.write(a.map(_ match {
				case s:String => 
					if (qoute) params.outQuote + s + params.outQuote
					else s
				case x => x.toString
			}).mkString(params.outSep) + "\n")
		
			
			
		writeRow(false)(
				"file",
				"rtStart",
				"rtApex",
				"rtEnd",
				"rtMs2",
				"nIsotopes",
				"nScans",
				"averagineCorr",
				"mzFeat",
				"mzPrec",
				"massFeat",
				"massFeatCalib",
				"massTheoretical",
				"massPrecursor",
				"z",
				"intensityApex",
				"intensitySum",
				"peptide",
				"protein"
			)
		
		for {
			(fileName, zMap) <- matches
			(z, pepFeats) <- zMap
			pepFeat <- pepFeats
			(p, f) <- pepFeat
		} {
			writeRow(true)(
				fileName,
				f.rtStart,
				f.rtApex,
				f.rtEnd,
				p.rt,
				f.nIsotopes,
				f.nScans,
				f.averagineCorr,
				f.mz,
				p.empMz,
				f.mass,
				f.massCalib,
				p.theoMass,
				p.empMass,
				f.z,
				f.intensityApex,
				f.intensitySum,
				peptides(p.peptide),
				proteins(p.protein)
			)
		}
		
		w.close()
	}
	
	
	
	def matchFeatures(
			feats:ArrayBuffer[DinoFeature], 
			peps:ArrayBuffer[TppPeptide]
	) = {
		for ((z, zPeps) <- peps.groupBy(_.z)) yield
			(z, matchZFeatures(feats.filter(_.z == z), zPeps))
	}
	
	
	def matchZFeatures(
			feats:ArrayBuffer[DinoFeature], 
			peps:ArrayBuffer[TppPeptide]
	) = {
		val fSort = feats.sortBy(_.rtStart)
		val pepSort = peps.sortBy(_.rt)
		
		var iFeatMin = 0
		var iFeatMax = 0
		for (i <- 0 until pepSort.length) yield {
			val pep = pepSort(i)
			while (iFeatMin < fSort.length && (fSort(iFeatMin).rtEnd + params.matchPostRT) < pep.rt)
				iFeatMin += 1
			iFeatMax = math.max(iFeatMax, iFeatMin)
			
			while (iFeatMax < fSort.length && (fSort(iFeatMax).rtStart - params.matchPreRT) < pep.rt)
				iFeatMax += 1
				
			
			val matchingFeats = new ArrayBuffer[DinoFeature]	
			for (j <- iFeatMin until iFeatMax) {
				val feat = fSort(j)
				if (pep.rt >= (feat.rtStart-params.matchPreRT) && pep.rt <= (feat.rtEnd+params.matchPostRT))
					if (within(pep.empMz, feat.mz, params.matchPPM))
						matchingFeats += feat
			}
			if (matchingFeats.isEmpty) None
			else Some(pep -> matchingFeats.minBy(f => math.abs(f.mz - pep.empMz)))
		}
	}
	
	
	def within(mz1:Double, mz2:Double, ppm:Double) = 
		2 * 1000000 * math.abs(mz1 - mz2) /	(mz1 + mz2) < ppm
	
	
	
	def readFeatures(f:File):HashMap[String, ArrayBuffer[DinoFeature]] = {
		var headerParsed = false
		
		var iMZ = -1
		var iZ = -1
		var iRT_START = -1
		var iRT_APEX = -1
		var iRT_END = -1
		var iN_ISO = -1
		var iN_SCAN = -1
		var iAVE_CORR = -1
		var iMASS = -1
		var iMASS_CAL = -1
		var iINT_APEX = -1
		var iINT_SUM = -1
		var iFILE = -1
		
		var res = new HashMap[String, ArrayBuffer[DinoFeature]]
		try {
			for (line <- Source.fromFile(f).getLines) {
				val cols = line.split("\t")map(_.trim)
				if (!headerParsed) {
					val header = cols
					iMZ = indexOfCol(header, "monoisoMz")
					iZ = indexOfCol(header, "charge")
					iRT_START = indexOfCol(header, "rtStart")
					iRT_APEX = indexOfCol(header, "rtApex")
					iRT_END = indexOfCol(header, "rtEnd")
					iN_ISO = indexOfCol(header, "nIsotopes")
					iN_SCAN = indexOfCol(header, "nScans")
					iAVE_CORR = indexOfCol(header, "averagineCorr")
					iMASS = indexOfCol(header, "mass")
					iMASS_CAL = indexOfCol(header, "massCalib")
					iINT_APEX = indexOfCol(header, "intensityApex")
					iINT_SUM = indexOfCol(header, "intensitySum")
					iFILE = indexOfCol(header, "file")
					headerParsed = true
				} else {
					val fileName = cols(iFILE).trim
					if (!res.contains(fileName))
						res(fileName) = new ArrayBuffer
					
					res(fileName) += DinoFeature(
						cols(iMZ).toDouble,
						cols(iZ).toInt,
						cols(iRT_START).toDouble * 60,
						cols(iRT_APEX).toDouble * 60,
						cols(iRT_END).toDouble * 60,
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
	
	
	def readSearchResults(
			f:File
	):(HashMap[String, ArrayBuffer[TppPeptide]], Array[String], Array[String]) = {
		var headerParsed = false
		
		var iPROB = -1
		var iSCAN = -1
		var iPEP = -1
		var iPROT = -1
		var iMASS_THEO = -1
		var iMZ_EMP = -1
		var iMASS_EMP = -1
		var iZ = -1
		var iRT = -1
		var iFILE = -1
		
		val peptides = new HashMap[String, Int]
		val proteins = new HashMap[String, Int]
		
		def pepToIndex(pep:String) = {
			if (!peptides.contains(pep))
				peptides(pep) = peptides.size
			peptides(pep)
		}
		def protToIndex(prot:String) = {
			if (!proteins.contains(prot))
				proteins(prot) = proteins.size
			proteins(prot)
		}
		
		var res = new HashMap[String, ArrayBuffer[TppPeptide]]
		try {
			for (line <- Source.fromFile(f).getLines) {
				val cols = line.split("\t")map(_.trim)
				if (!headerParsed) {
					val header = cols
					iPROB = indexOfCol(header, "probability")
					iSCAN = indexOfCol(header, "start_scan")
					iPEP = indexOfCol(header, "peptide")
					iPROT = indexOfCol(header, "protein")
					iMASS_THEO = indexOfCol(header, "calc_neutral_pep_mass")
					iMZ_EMP = indexOfCol(header, "MZratio")
					iMASS_EMP = indexOfCol(header, "precursor_neutral_mass")
					iZ = indexOfCol(header, "assumed_charge")
					iRT = indexOfCol(header, "retention_time_sec")
					iFILE = indexOfCol(header, "file")
					headerParsed = true
				} else {
					val fileName = cols(iFILE).trim
					if (!res.contains(fileName))
						res(fileName) = new ArrayBuffer
					
					res(fileName) += TppPeptide(
						cols(iPROB).toDouble,
						cols(iSCAN).toInt,
						pepToIndex(cols(iPEP)),
						protToIndex(cols(iPROT)),
						cols(iMASS_THEO).toDouble,
						cols(iMZ_EMP).toDouble,
						cols(iMASS_EMP).toDouble,
						cols(iZ).toInt,
						cols(iRT).toDouble
					)
				}
			}
		} catch {
			case e:Exception => {}
		}
		
		
		val pepArr = new Array[String](peptides.size)
		for ((pep, i) <- peptides) pepArr(i) = pep
		
		val protArr = new Array[String](proteins.size)
		for ((prot, i) <- proteins) protArr(i) = prot
		
		(res, pepArr, protArr)
	}
}