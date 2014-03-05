/**
 * 
 */
package org.rcsb.codec;

import static org.rcsb.codec.CodecConstants.AMINO_ACID;
import static org.rcsb.codec.CodecConstants.BFACTOR;
import static org.rcsb.codec.CodecConstants.BYTE2_INTEGER_MARKER;
import static org.rcsb.codec.CodecConstants.BYTE2_ENCODED_MARKER;
import static org.rcsb.codec.CodecConstants.BYTE4_SHORT_MARKER;
import static org.rcsb.codec.CodecConstants.BYTE4_INTEGER_MARKER;
import static org.rcsb.codec.CodecConstants.BYTE4_ENCODED_MARKER;
import static org.rcsb.codec.CodecConstants.SHORT_COORDINATE_TYPE;
import static org.rcsb.codec.CodecConstants.INTEGER_COORDINATE_TYPE;
import static org.rcsb.codec.CodecConstants.ENCODED_COORDINATE_TYPE;
import static org.rcsb.codec.CodecConstants.BO_SCALE;
import static org.rcsb.codec.CodecConstants.BO_PRECISION;
import static org.rcsb.codec.CodecConstants.CHAIN;
import static org.rcsb.codec.CodecConstants.COORD;
import static org.rcsb.codec.CodecConstants.END;
import static org.rcsb.codec.CodecConstants.GINFO;
import static org.rcsb.codec.CodecConstants.GROUP;
import static org.rcsb.codec.CodecConstants.HEAD;
import static org.rcsb.codec.CodecConstants.MODEL;
import static org.rcsb.codec.CodecConstants.NUCLEOTIDE;
import static org.rcsb.codec.CodecConstants.NUCLEOTIDE_BOND_LENGTH;
import static org.rcsb.codec.CodecConstants.NUCLEOTIDE_TAIL_ATOM_NAME;
import static org.rcsb.codec.CodecConstants.OCCUPANCY;
import static org.rcsb.codec.CodecConstants.PEPTIDE_BOND_LENGTH;
import static org.rcsb.codec.CodecConstants.PEPTIDE_TAIL_ATOM_NAME;
import static org.rcsb.codec.CodecConstants.SEQUENCE;
import static org.rcsb.codec.CodecConstants.STRUCTURE;
import static org.rcsb.codec.CodecConstants.TAIL;
import static org.rcsb.codec.CodecConstants.XYZ_PRECISION;

import java.io.DataInputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * StructureDecoderImpl1 decodes the data section of a High-Efficiency Structure Codec (HESC) encoded byte array (compression method 1). 
 * The StructureInflatorInterface must be implemented to use this class.
 * 
 * The data section of HESC is a set of records, terminated by the END record.
 * The following records have the general format: record id, record length in number of bytes, data ..
 * 
 *                         +----------+---------------+---+---+---+---+---+---+---+--
 *                         |record id | record length | data ...
 *                         +----------+---------------+---+---+---+---+---+---+---+--
 *  STRUCTURE              |        s |          byte | model count (int), homogeneousModel (boolean)
 *  MODEL                  |        m |          byte | chain count (int)
 *  CHAIN                  |        c |          byte | sequenceIndex (int), chain id (4 bytes), groupCount (int)
 *  GROUP                  |        g |          byte | groupIndex (int), [groupNumber (int)]
 *  GINFO                  |        I |           int | atomCount (short), flag (1 byte), group information (len - 4*atomCount - 3 bytes), bond information (4*atomCount shorts)
 *  SEQUENCE               |        Q |           int | sequence string using 1-letter codes (note, this record is not used currently)
 *  COORD                  |        X |           int | atom coordinates encoded as list of integers and shorts
 *  BFACTOR                |        T |           int | b factors encoded as list of integer and shorts
 *  OCCUPANCY              |        O |           int | occupancies (len/2 shorts)
 *  END                    |        e |          none | none
 *  
 *  Note: 
 *     Lower case record ids use 1 byte for the record length 
 *     Upper case record ids use 4 bytes (int) for the record length
 *     Data items in [ ] are optional, the record length will indicate if these data are present
 *                                  
 * @author Peter Rose
 *
 */
public class StructureDecoderImpl1 extends StructureDecoder {
	private DataInputStream inStream = null;

	// byte arrays for temporary data
	private byte[] b1 = new byte[1];
	private byte[] b2 = new byte[2];
	private byte[] b3 = new byte[3];
	private byte[] b4 = new byte[4];
	private int[] buffer = new int[4];

	private int modelCount = 0;
	private int sequenceIndex = 0;
	private int chainCount = 0;
	private int groupNumber = 0;
	private int intType = 4;
	private int byteOffset = 0;
	
	private boolean homogeneousModel = false;
	private List<String> sequences = new ArrayList<String>();
	private List<Integer> sequenceIndices = new ArrayList<Integer>();
	private List<String> chainIds = null;
	private List<Integer> chainCounts = new ArrayList<Integer>();
	private List<Integer> groupCounts = new ArrayList<Integer>();
	private List<Integer> groupIndices = new ArrayList<Integer>();
	private List<Integer> groupNumbers = new ArrayList<Integer>();
	private List<String[]> groupInfo = new ArrayList<String[]>();
	private List<int[]> bondInfo = new ArrayList<int[]>();
	private List<Byte> flagInfo = new ArrayList<Byte>();
	private List<Integer> bFactors = new ArrayList<Integer>();
	private int[] occupancy = null;
	
	private boolean useBfactor = false;
	private boolean useOccupancy = false;

	private StructureInflatorInterface inflator;

	private static final int MAX_GROUP_SIZE = 1024;

	public StructureDecoderImpl1(DataInputStream dataInputStream, StructureInflatorInterface inflator) {
		this.inStream = dataInputStream;
		this.inflator = inflator;
	}
	
	public void decode() throws IOException {	
		modelCount = 0;
		chainCount = 0;
		homogeneousModel = false;
		
		byte recordId = 0;
		
		while ((recordId = (byte) inStream.readByte()) != END) {

			switch (recordId) {
			case STRUCTURE:
				readStructureRecord();
				break;
			case MODEL:
				readModelRecord();
				break;
			case SEQUENCE:
				readSequenceRecord();
				break;
			case CHAIN:
				readChainRecord();
				break;
			case GROUP:
				readGroupRecord();
				break;
			case GINFO:
				readGInfoRecord();
				break;
			case COORD:
				readCoordRecord();
				break;
			case BFACTOR:
				readBFactorRecord();
				break;
			case OCCUPANCY:
				readOccupancyRecord();
				break;
			default:
				throw new IOException("StructureDecoder: Invalid record: " + recordId);
			}
		}
	}
	
	private void readSequenceRecord() throws IOException {
		int len = inStream.readInt();
		String sequence = readFixedLengthString(len);
//		System.out.println("readSequenceRecord: " + sequence);
		sequences.add(sequence);
	}

	private void readStructureRecord() throws IOException {
	//	System.out.println("readStructureRecord");
		inStream.skipBytes(1);
		modelCount = inStream.readInt();	
		homogeneousModel = inStream.readBoolean();
	//	System.out.println("modelCount: " + modelCount);
	//	System.out.println("homogeneous: " + homogeneousModel);
		chainIds = new ArrayList<String>(chainCount);
		groupCounts = new ArrayList<Integer>(chainCount);
	}
	
	private void readModelRecord() throws IOException {
		inStream.skipBytes(1);
		chainCounts.add(inStream.readInt());
	}
	
	private void readChainRecord() throws IOException {
//		System.out.println("readChainRecord");
		inStream.skipBytes(1);
		sequenceIndex = inStream.readInt();
		sequenceIndices.add(sequenceIndex);
		String chainId = readFixedLengthString(4);
		chainIds.add(chainId);
		int groupCount = inStream.readInt();	
		groupCounts.add(groupCount);
		groupNumber = 0;	
	}
	
	private void readGroupRecord() throws IOException {
		int len = inStream.readByte();
		int groupIndex = inStream.readInt();
		groupIndices.add(groupIndex);

		if (len == 8) {
			groupNumber = inStream.readInt();
			groupNumbers.add(groupNumber);
		} else {
			// if no explicit group number is given, then the group number are sequential
			groupNumbers.add(++groupNumber);
		}
	}
	
	private void readGInfoRecord() throws IOException {
//		System.out.println("readGInfoRecord");
		int len = inStream.readInt();
		int aCount = inStream.readShort();
		len -= 4*aCount + 3;
		byte flags = inStream.readByte();	
		String[] info = readGroupInfo(len);
		int[] bondList = new int[aCount*2];
		for (int bl = 0; bl < aCount*2; bl++) {
			bondList[bl] = (int)inStream.readShort();
		}
		groupInfo.add(info);
		bondInfo.add(bondList);
		flagInfo.add(flags);
	}
	
	private void readBFactorRecord() throws IOException {
//		System.out.println("readBFactorRecord");
		int len = inStream.readInt();
		byteOffset = 0;
		intType = 4;
		while (byteOffset < len) {
			int val = readNextInt();
			bFactors.add(val);
		}
		useBfactor = true;
	}
	
	private void readOccupancyRecord() throws IOException {
	//	System.out.println("readOccupancyRecord");
		int len = inStream.readInt();	
		int n = len/2;
		occupancy = new int[n];
		for (int i = 0; i < n; i++) {
			occupancy[i] = inStream.readShort();
		}
		useOccupancy = true;
	}

	private void readCoordRecord() throws IOException {
	//			System.out.println("readCoordRecord");
		inStream.skipBytes(4);
		int[] x = new int[MAX_GROUP_SIZE];
		int[] y = new int[MAX_GROUP_SIZE];
		int[] z = new int[MAX_GROUP_SIZE];
		int[] b = new int[MAX_GROUP_SIZE];

		inflator.setModelCount(modelCount);
		
		byteOffset = 0;
		intType = 4;

		int chainIndex = 0;
		int groupIndex = 0;
		int atomIndex = 0;

		for (int m = 0; m < modelCount; m++) {
			if (homogeneousModel) {
				chainIndex = 0;
				groupIndex = 0;
				chainCount = chainCounts.get(0);
			} else {
			    chainCount = chainCounts.get(m);
			}
			
			inflator.setModelInfo(m, chainCount);

			for (int i = 0; i < chainCount; i++) {	
				// integer atom coordinates and b-factor
				int xOffset = 0;
				int yOffset = 0;
				int zOffset = 0;
				int bOffset = 0;
				
				// integer atom coordinates and b-factor for the 
				// polymer tail atom from the previous group (residue)
				int xTail = 0;
				int yTail = 0;
				int zTail = 0;
				int bTail = 0;
				boolean hasTail = false;
						
				String chainId = chainIds.get(chainIndex);
				int seqIndex = sequenceIndices.get(chainIndex); // not used currently
				int groupCount = groupCounts.get(chainIndex);
				
				inflator.setChainInfo(chainId, groupCount);
				chainIndex++;

				for (int j = 0; j < groupCount; j++) {
					groupNumber = groupNumbers.get(groupIndex);
					int gIndex = groupIndices.get(groupIndex);
					groupIndex++;

					String[] info = groupInfo.get(gIndex);
					int[] bondList = bondInfo.get(gIndex);
					int atomCount = bondList.length/2;
	//				System.out.println("Getting goup info: atomCount: " + atomCount);
					byte flags = flagInfo.get(gIndex);

					boolean isAminoAcid = (flags & AMINO_ACID) != 0;
					boolean isNucleotide = (flags & NUCLEOTIDE) != 0;
				
					boolean hasHead =  (flags & HEAD) != 0;
					if (! hasTail) {
						xTail = 0;
						yTail = 0;
						zTail = 0;
						bTail = 0;
					}

					int index = 0;
					String groupName = info[index++];
					char insertionCode = info[index++].charAt(0);

					int polymerType = 0;
					if (isAminoAcid) {
						polymerType = 1;
					} else if (isNucleotide) {
						polymerType = 2;
					} 
					
					inflator.setGroupInfo(groupName, groupNumber, insertionCode, polymerType, atomCount);

					// no group is larger than MAX_GROUP_SIZE, but
					// just in case that ever changes, accommodate large sizes
					if (atomCount > MAX_GROUP_SIZE) {
						x = new int[atomCount];
						y = new int[atomCount];
						z = new int[atomCount];
						b = new int[atomCount];
					}

					for (int k = 0; k < atomCount; k++) {
						String atomName = info[index++];
						String atomNameTrimmed = atomName.trim();
						String element = info[index++];
						char altLoc = info[index++].charAt(0);

						int bondLength = 0;
						if (bondList[k] >=0) {
							bondLength = bondList[atomCount+k];
						} else if (k == 0 && isAminoAcid && hasHead && hasTail) {
							bondLength = PEPTIDE_BOND_LENGTH;
						} else if (k == 0 && isNucleotide && hasHead && hasTail) {
							bondLength = NUCLEOTIDE_BOND_LENGTH;
						}

						int[] xyz = decodeCoords(bondLength);
						
						if (bondList[k] >=0) {
							// use coordinates from a previous atom in this group
							int reference = bondList[k];
							xOffset = x[reference];
							yOffset = y[reference];
							zOffset = z[reference];
							bOffset = b[reference];
						} else if (k == 0 && hasTail && hasHead) {
							xOffset = xTail;
							yOffset = yTail;
							zOffset = zTail;
							bOffset = bTail;
						}

						xOffset += xyz[0];
						yOffset += xyz[1];
						zOffset += xyz[2];
						if (useBfactor) {
							bOffset += bFactors.get(atomIndex);
						}

						int occ = BO_SCALE;
						if (useOccupancy) {
							occ = occupancy[atomIndex];
						}
						atomIndex++;

						x[k] = xOffset;
						y[k] = yOffset;
						z[k] = zOffset;
						b[k] = bOffset;

						inflator.setAtomInfo(atomName, 0, altLoc, xOffset*XYZ_PRECISION, yOffset*XYZ_PRECISION, zOffset*XYZ_PRECISION, occ*BO_PRECISION, bOffset*BO_PRECISION, element);

						if ((isAminoAcid && atomNameTrimmed.equals(PEPTIDE_TAIL_ATOM_NAME) || (isNucleotide && atomNameTrimmed.equals(NUCLEOTIDE_TAIL_ATOM_NAME)))) {
							xTail = xOffset;
							yTail = yOffset;
							zTail = zOffset;
							bTail = bOffset;
						}			
					}
					
					hasTail =  (flags & TAIL) != 0;
				}
			}
		}
	}
	
	private int readNextInt() throws IOException {
		int v = 0;
		switch (intType) {
		case SHORT_COORDINATE_TYPE:
			v = inStream.readShort();
			byteOffset += 2;
			
			switch (v) {
			case BYTE2_INTEGER_MARKER: 
				intType = INTEGER_COORDINATE_TYPE;
				return readNextInt();
			case BYTE2_ENCODED_MARKER: 
				intType = ENCODED_COORDINATE_TYPE;
				return readNextInt();
			}
			break;

		case INTEGER_COORDINATE_TYPE:
			v = inStream.readInt();
			byteOffset += 4;

			switch (v) {
			case BYTE4_SHORT_MARKER: 
				intType = SHORT_COORDINATE_TYPE;
				return readNextInt();
			case BYTE4_ENCODED_MARKER: 
				intType = ENCODED_COORDINATE_TYPE;
				return readNextInt();
			}
			break;
			
		case ENCODED_COORDINATE_TYPE:
			v = inStream.readInt();
			byteOffset += 4;

			switch (v) {
			case BYTE4_SHORT_MARKER: 
				intType = SHORT_COORDINATE_TYPE;
				return readNextInt();
			case BYTE4_INTEGER_MARKER: 
				intType = INTEGER_COORDINATE_TYPE;
				return readNextInt();
			}
			break;
		}

		return v;
	}
	
	private int[] decodeCoords(int bondLength) throws IOException {
		int v = readNextInt();
		if (intType == ENCODED_COORDINATE_TYPE) {
			// decodes the deltaX (buffer[0]), deltaY (buffer[1]), and deltaZ (buffer[2]) coordinates
			// from a single 32-bit integer value, using the standard bond length as a parameter
			BitEncoder.fromInt(v, bondLength, b4, buffer);
		} else {
			// decodes the deltaX (buffer[0]), deltaY (buffer[1]), and deltaZ (buffer[2]) coordinates
			// from three integer values
			buffer[0] = v;
			buffer[1] = readNextInt();
			buffer[2] = readNextInt();
		}
		return buffer;
	}

	private String[] readGroupInfo(int length) throws IOException {
		int len = 2 + 3 * (length - 4)/7;
		String[] info = new String[len];

		inStream.read(b3); 
		info[0] = new String(b3); // group name
		inStream.read(b1); 
		info[1] = new String(b1); // insertion code
		for (int i = 2; i < len; i += 3) {
			inStream.read(b4);
			info[i] = new String(b4); // atom name including spaces
			inStream.read(b2);
			info[i+1] = new String(b2).trim(); // element name
			inStream.read(b1);
			info[i+2] = new String(b1); // alternative location
		}
		return info;
	}

	private String readFixedLengthString(int length) throws IOException {
		byte[] bytes = new byte[length];
		inStream.read(bytes);
		return new String(bytes);
	}
}

