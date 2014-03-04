package org.rcsb.codec;

public final class CodecConstants {
	public static final String MAGIC_NUMBER = "HESC"; // High Efficiency Structure Codec
	public static final String CODEC_FILE_EXTENSION = ".hesc";
	
	// The following two arrays list the supported versions major/minor version numbers
	// Both arrays must of the same size
	public static final byte[] MAJOR_VERSIONS = {0};
	public static final byte[] MINOR_VERSIONS = {0};
	
	public static final byte AMINO_ACID = 1;
	public static final byte NUCLEOTIDE = (1 << 1);
	public static final byte HETATOM = (1 << 2);
	public static final byte HEAD = (1 << 3); // flag indicates that monomer has a head atom (N or P)
	public static final byte TAIL = (1 << 4); // flag indicates that monomer has a tail atom (C or O3')

	public static int XYZ_PRECISION = 1000;
	public static float XYZ_SCALE = 0.001f;
	public static int B_PRECISION = 100;
	public static float B_SCALE = 0.01f;

	public static final int PEPTIDE_BOND_LENGTH = (int)Math.round(1.325 * XYZ_PRECISION);
	public static final int NUCLEOTIDE_BOND_LENGTH = (int)Math.round(1.6 * XYZ_PRECISION);
	
	public static final String PEPTIDE_HEAD_ATOM_NAME = "N";
	public static final String PEPTIDE_TAIL_ATOM_NAME = "C";
	public static final String NUCLEOTIDE_HEAD_ATOM_NAME = "P";
	public static final String NUCLEOTIDE_TAIL_ATOM_NAME = "O3'";
	
	public static final int BYTE2_MIN_VALUE = Short.MIN_VALUE;
	public static final int BYTE2_MAX_VALUE = Short.MAX_VALUE-5;
	public static final int BYTE4_MAX_VALUE = Integer.MAX_VALUE-5;
	
	public static final int SHORT_COORDINATE = 2;
	public static final int INTEGER_COORDINATE = 4;
	public static final int ENCODED_COORDINATE = 5;

	public static final int BYTE2_SHORT_MARKER = BYTE2_MAX_VALUE + SHORT_COORDINATE;
	public static final int BYTE2_INTEGER_MARKER = BYTE2_MAX_VALUE + INTEGER_COORDINATE;
	public static final int BYTE2_ENCODED_MAKRKER = BYTE2_MAX_VALUE + ENCODED_COORDINATE;

	public static final int BYTE4_SHORT_MARKER = BYTE4_MAX_VALUE + SHORT_COORDINATE;
	public static final int BYTE4_INTEGER_MARKER = BYTE4_MAX_VALUE + INTEGER_COORDINATE;
	public static final int BYTE4_ENCODED_MAKRKER = BYTE4_MAX_VALUE + ENCODED_COORDINATE;

	public static final byte BLANK = 32;
	
	// Structure record identifiers. Record identifiers in lower case have a record length
	// specified by a 1-byte values, whereas record identifiers in uppper case have
	// the record length specified as an integer (4 bytes)
	public static final byte STRUCTURE = 's';
	public static final byte MODEL = 'm';
	public static final byte CHAIN = 'c';
	public static final byte SEQUENCE = 'Q';
	public static final byte GINFO = 'I';
	public static final byte GROUP = 'g';
	public static final byte ATOM = 'a';
	public static final byte BOND = 'B';
	public static final byte COORD = 'X';
	public static final byte BFACTOR = 'T';
	public static final byte OCCUPANCY = 'O';
	public static final byte END = 'e';
}
