package se.lth.immun.protocol

import java.io.InputStream
import java.io.BufferedInputStream
import java.io.File
import java.io.FileInputStream

import com.google.protobuf.CodedInputStream

class ProtoBufferedFileInputStream(size:Int, file:File) {

	val buffer = new Array[Byte](size)
	val bis = new BufferedInputStream(new FileInputStream(file))
	var cis = CodedInputStream.newInstance(buffer)
	var atEOF = false
	
	def ensure(n:Int):Boolean = {
		if (atEOF) return false
		var bytesRead = 0
		var readStatus = 0
		while (bytesRead < n && readStatus >= 0) {
			readStatus = bis.read(buffer, bytesRead, n-bytesRead)
			if (readStatus >= 0)
				bytesRead += readStatus
			else atEOF = true
		}
		cis = CodedInputStream.newInstance(buffer)
		cis.pushLimit(bytesRead)
		return n == bytesRead
	}
	
	def close = {
		bis.close
	}
}