/**
 * 
 */
package org.rcsb.codec;

import java.io.BufferedInputStream;
import java.io.ByteArrayInputStream;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.zip.GZIPInputStream;
import static org.rcsb.codec.CodecConstants.*;

/**
 * StructureInflator decodes a HESC file. The HESC file is a gzipped byte array with the following structure:
 * 
 * Header:
 * Magic number        : 4 bytes
 * Major version number: 1 byte
 * Minor version number: 1 byte
 * Compression method  : 1 byte
 * 
 * Data records:
 *        byte        1                  1          n bytes
 *            +---------------------+---------------+-----
 *            |lower case record id | record length | data ..
 *            +---------------------+---------------+-----
 *            
 *        byte        1                  4          n bytes
 *            +---------------------+---------------+-----
 *            |upper case record id | record length | data ..
 *            +---------------------+---------------+-----
 * @author Peter
 *
 */
public class StructureInflator {	
	private DataInputStream inStream;
	private StructureInflatorInterface inflator;
	
	private byte majorVersion;
	private byte minorVersion;
	private int compressionMethod;

	private long fileSize = 0;
	private long fileSizeCompressed = 0;
	private long readTime = 0;

	public StructureInflator(StructureInflatorInterface inflator) {
		this.inflator = inflator;
	}
	
	/**
	 * @return the fileSize
	 */
	public long getFileSize() {
		return fileSize;
	}

	/**
	 * @return the fileSizeCompressed
	 */
	public long getFileSizeCompressed() {
		return fileSizeCompressed;
	}

	/**
	 * @return the readWriteTime
	 */
	public long getReadTime() {
		return readTime;
	}
	
	public void read(String fileName) throws Exception {
		fileSize = 0;
		fileSizeCompressed = 0;
		readTime = 0;
		
		long start = System.nanoTime();
		openFile(fileName);		
		readHeader();
		readData();
		inStream.close();
		
		readTime = System.nanoTime() - start;

	}
	
	/**
	 * Inflates a structure from an unbuffered input stream
	 * @param inputStream an unbuffered input stream
	 * @throws Exception
	 */
	
	public void read(InputStream inputStream) throws Exception {
		if (inputStream instanceof BufferedInputStream) {
			throw new IOException("StructureInflator: the input stream should not be buffered.");	
		}
		fileSize = 0;
		fileSizeCompressed = 0;
		readTime = 0;
		
		long start = System.nanoTime();
		
		inStream = new DataInputStream(new BufferedInputStream(new GZIPInputStream(inputStream)));
		readHeader();
		readData();
		inStream.close();
		
		readTime = System.nanoTime() - start;
	}
	
	public void read (byte[] data) throws Exception {
		if (data == null) {
			throw new IOException();
		}
		fileSize = 0;
		fileSizeCompressed = 0;
		readTime = 0;
		
		long start = System.nanoTime();
		
		openByteArray(data);		
		readHeader();
		readData();
		inStream.close();
		
		readTime = System.nanoTime() - start;
	}
	
	public void readNext () throws IOException {	
		long start = System.nanoTime();
		
		readHeader();	
		readData();

		long time = System.nanoTime() - start;
		readTime += time;
	}

	private void openFile(String fileName) throws Exception {
		if (! fileName.endsWith(CODEC_FILE_EXTENSION)) {
			throw new IOException("StructureInflator: File name has invalid extension: " + fileName + ". File extension .hesc required");
		}

		File file = new File(fileName);
		fileSizeCompressed = file.length();
        // TODO should check size
		byte[] data = new byte[(int) fileSizeCompressed];
		DataInputStream dis = new DataInputStream(new FileInputStream(file));
		dis.readFully(data);
		dis.close();
		openByteArray(data);
	}
	
	private void openByteArray(byte[] data) throws Exception {
		fileSizeCompressed = data.length;
		inStream = new DataInputStream(new BufferedInputStream(new GZIPInputStream(new ByteArrayInputStream(data), 8192)));
	}
	
	private void readHeader() throws IOException {
		String magicNumber = readMagicNumber();
		if (! magicNumber.equals(MAGIC_NUMBER)) {
			throw new IOException("Invalid file format: magic number is: " + magicNumber +" Expected: " + MAGIC_NUMBER);
		}
		
	    majorVersion = inStream.readByte();
		minorVersion = inStream.readByte();

		boolean validVersion = false;
		
		for (int i = 0; i < MAJOR_VERSIONS.length; i++) {
	
		   if (majorVersion == MAJOR_VERSIONS[i] && minorVersion == MINOR_VERSIONS[i]) {
			   validVersion = true;
			   break;
		   }
			
		}

		if (! validVersion) {
			throw new IOException("Invalid file format: version: " + majorVersion + "." + minorVersion);
		}

		compressionMethod = inStream.readByte();
	}

	private String readMagicNumber() throws IOException {
		byte[] bytes = new byte[MAGIC_NUMBER.length()];
		inStream.read(bytes);
		return new String(bytes);
	}
	
	private void readData() throws IOException {
		StructureDecoder decoder = StructureDecoder.getDecoder(minorVersion, majorVersion, compressionMethod, inStream, inflator);
		decoder.decode();
	}
	
	public void close() throws IOException {
		if (inStream != null) {
			inStream.close();
		}
	}
}
