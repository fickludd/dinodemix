package se.lth.immun.protocol

import java.io.File
import java.io.BufferedInputStream
import java.io.FileInputStream
import java.io.BufferedOutputStream
import java.io.FileOutputStream

import se.lth.immun.protocol.MSFeatureProtocol.MsgSize
import se.lth.immun.protocol.MSFeatureProtocol.RtMap
import se.lth.immun.protocol.MSFeatureProtocol.RtUnit
import se.lth.immun.protocol.MSFeatureProtocol.Feature
import se.lth.immun.protocol.MSFeatureProtocol.Hill

import com.google.protobuf.CodedInputStream

import scala.collection.JavaConversions._
import scala.collection.mutable.ArrayBuffer



case class MsFeatures(rtMap:Seq[Double], features:Seq[MsFeature])
case class MsFeature(
		val mz:Double,
		val z:Int,
		val mass:Double,
		val rtApex:Double,
		val rtStart:Double,
		val rtEnd:Double,
		val intensityApex:Double,
		val intensitySum:Double,
		val averagineCorr:Option[Double],
		val hills:Seq[MsHill]
)

case class MsHill(
		val startIndex:Int,
		val endIndex:Int,
		val mz:Double,
		val mzErr:Double,
		val fwhm:Option[Double],
		val rtApex:Double,
		val intensityApex:Double,
		val intensity:Seq[Double]
)



object MsFeatureFile {

	def read(f:File, verbose:Boolean) = {
		val r = new ProtoBufferedFileInputStream(2048*2048, f)
		
		if (!parseAndCodeSize(r, verbose))
			throw new Exception("Error reading MsgSize!")
		
		val rtMap = RtMap.parseFrom(r.cis)
		val mult = rtMap.getUnit match {
			case RtUnit.MINUTE => 1.0
			case RtUnit.SECOND => 1 / 60.0
		}
		val rts = rtMap.getRtList.map(_ * mult)
		
		val features = new ArrayBuffer[MsFeature]
		while (parseAndCodeSize(r, verbose)) 
			features += readMsFeature(r, rts)
		
		MsFeatures(rts, features)
	}
	
	def readMsFeature(r:ProtoBufferedFileInputStream, rtMap:Seq[Double]):MsFeature = {
		val f = Feature.parseFrom(r.cis)
		val hills = f.getHillList.map(toMsHill(_, rtMap))
		val rtApex = 
			if (f.hasRtApex) f.getRtApex 
			else hills.maxBy(_.intensityApex).rtApex
		val intensityApex =
			if (f.hasIntensityApex) f.getIntensityApex
			else hills.map(_.intensityApex).max
		val intensitySum =
			if (f.hasIntensitySum) f.getIntensitySum
			else hills.map(_.intensity.sum).sum
		val averagineCorr =
			if (f.hasAveragineCorr) Some(f.getAveragineCorr.toDouble)
			else None
		MsFeature(
			f.getMz, f.getZ, f.getMass, 
			rtApex, 
			hills.map(h => rtMap(h.startIndex)).min,
			hills.map(h => rtMap(h.endIndex)).max,
			intensityApex, intensitySum, 
			averagineCorr, hills)
	}
	
	def toMsHill(hill:Hill, rtMap:Seq[Double]):MsHill = {
		val intensity = hill.getIntensityList.map(_.toDouble)
		val maxLocalIndex = (0 until intensity.length).maxBy(intensity)
		val rtApex =
			if (hill.hasRtApex) hill.getRtApex 
			else rtMap(hill.getStartIndex + maxLocalIndex)
		val intensityApex =
			if (hill.hasIntensityApex) hill.getIntensityApex 
			else intensity(maxLocalIndex).toDouble
		MsHill(
				hill.getStartIndex,
				hill.getEndIndex,
				hill.getMz,
				hill.getMzErr,
				if (hill.hasFwhm) Some(hill.getFwhm.toDouble) else None,
				rtApex,
				intensityApex,
				intensity)
	}
	
	def parseAndCodeSize(r:ProtoBufferedFileInputStream, verbose:Boolean) = {
		if (r.ensure(5)) {
			val n = MsgSize.parseFrom(r.cis).getSize
			if (verbose)
				println("reading msg of n=%d bytes".format(n))
			r.ensure(n)
		} else false
	}
	
	def write(f:File, features:MsFeatures, rtUnit:RtUnit = RtUnit.MINUTE, verbose:Boolean) = {
		val w = new BufferedOutputStream(new FileOutputStream(f))
		val rtMap = RtMap.newBuilder.setUnit(rtUnit)
		for (rt <- features.rtMap) rtMap.addRt(rt)
		
		val rtMapMsg = rtMap.build
		writeMsg(w, rtMapMsg.toByteArray, verbose)
		if (verbose)
			println(rtMapMsg.toString)
		for (f <- features.features) {
			val featureMsg = toFeature(f)
			writeMsg(w, featureMsg.toByteArray, verbose)
			if (verbose)
				println(featureMsg.toString)
		}
		w.close
	}
	
	def toFeature(f:MsFeature):Feature = {
		val b = Feature.newBuilder.setMz(f.mz)
			.setZ(f.z)
			.setMass(f.mass)
			.setRtApex(f.rtApex.toFloat)
			.setIntensityApex(f.intensityApex.toFloat)
			.setIntensitySum(f.intensitySum.toFloat)
		f.averagineCorr.foreach(ac => b.setAveragineCorr(ac.toFloat))
		for (h <- f.hills) b.addHill(toHill(h))
		b.build
	}
	
	def toHill(h:MsHill):Hill = {
		val b = Hill.newBuilder.setStartIndex(h.startIndex)
			.setEndIndex(h.endIndex)
			.setMz(h.mz)
			.setMzErr(h.mzErr)
			.setRtApex(h.rtApex.toFloat)
			.setIntensityApex(h.intensityApex.toFloat)
		h.fwhm.foreach(x => b.setFwhm(x.toFloat))
		for (x <- h.intensity) b.addIntensity(x.toFloat)
		b.build
	}
	
	def writeMsg(w:BufferedOutputStream, bytes:Array[Byte], verbose:Boolean = false) = {
		val n = bytes.length
		val sizeMsgBytes = MsgSize.newBuilder().setSize(n).build.toByteArray
		if (verbose)
			println("writing msg, n(bytes): %d+%d".format(sizeMsgBytes.length, n))
		w.write(sizeMsgBytes)
		w.write(bytes)
	}
}