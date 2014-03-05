package org.rcsb.codec;

public final class CodecConstants {
	/**
	 * The magic number must appear as the first 4 bytes of the 
	 * High Efficiency Structure Codec (HESC) byte array. Note, .hesc files
	 * have been compressed with standard compression algorithms (e.g., gzip). The HESC
	 * magic number will be the first 4 bytes after decompression by the standard 
	 * compression tool (gzip). The compression tool may vary depending on the
	 * codec version or compression method.
	 */
	public static final String MAGIC_NUMBER = "HESC";
	
	/**
	 * File extension for an HESC files.
	 */
	public static final String CODEC_FILE_EXTENSION = ".hesc";

	/**
	 * Array of supported major version numbers.
	 */
	public static final byte[] MAJOR_VERSIONS = {0};
	
	/**
	 * Array of supported minor version numbers. Note, the major and minor version
	 * number must be specified in pairs in both array. The length of both arrays
	 * must be identical.
	 */
	public static final byte[] MINOR_VERSIONS = {0};
	
	/**
	 * A flag that indicates an amino acid group (residue).
	 */
	public static final byte AMINO_ACID = 1; 
	
	/** 
	 * A flag that indicates a nucleotide group (residue).
	 */
	public static final byte NUCLEOTIDE = (1 << 1);
	
	/**
	 * A flag that indicates a non-polymer group (residue).
	 */
	public static final byte NON_POLYMER = (1 << 2);
	
	/**
	 * A flag that indicates the presence of a polymer head atom in a group (residue)
	 * (N for amino acid, P for nucleotide).
	 */
	public static final byte HEAD = (1 << 3);
	
	/**
	 * A flag that indicates the presence of a polymer tail atom in a group (residue)
	 * (C for amino acid, O3' for nucleotide).
	 */
	public static final byte TAIL = (1 << 4);

	/**
	 * A factor that converts the x, y, and z coordinates of an atom to an integer without loss of precision.
	 */
	public static final int XYZ_SCALE = 1000;
	
	/**
	 * The precision of the x, y, and z coordinates.
	 */
	public static final float XYZ_PRECISION = 0.001f;
	
	/** 
	 * A factor that converts the b factor (temperature factor) and occupancy value to an integer without loss of precision.
	 */
	public static final int BO_SCALE = 100;
	
	/**
	 * The precision of a b factor (temperature factor) and occupancy value. Note, new mmCIF files support a precision of 0.001
	 * TODO need to check precision dynamically. It can't be a constant going forward.
	 */
	public static final float BO_PRECISION = 0.01f;

	/**
	 * Standard peptide bond length (N - C) expressed as an integer.
	 */
	public static final int PEPTIDE_BOND_LENGTH = (int)Math.round(1.325 * XYZ_SCALE);
	
	/**
	 * Standard nucleotide bond length (P - O3') expressed as an integer.
	 */
	public static final int NUCLEOTIDE_BOND_LENGTH = (int)Math.round(1.6 * XYZ_SCALE);
	
	/**
	 * Atom name of the polymer head atom of an amino acid.
	 */
	public static final String PEPTIDE_HEAD_ATOM_NAME = "N";
	
	/**
	 * Atom name of the polymer tail atom of an amino acid.
	 */
	public static final String PEPTIDE_TAIL_ATOM_NAME = "C";
	
	/**
	 * Atom name of the polymer head atom of a nucleotide.
	 */
	public static final String NUCLEOTIDE_HEAD_ATOM_NAME = "P";
	
	/**
	 * Atom name of the polymer tail atom of a nucleotide.
	 */
	public static final String NUCLEOTIDE_TAIL_ATOM_NAME = "O3'";
	
	/**
	 * Minimum value of a short (2 byte signed integer) used in the encoding of integer values.
	 */
	public static final int BYTE2_MIN_VALUE = Short.MIN_VALUE;
	
	/** 
	 * Maximum value of a short (2 byte signed integer) used in the encoding of integer values.
	 */
	public static final int BYTE2_MAX_VALUE = Short.MAX_VALUE-5;
	
	/** 
	 * Maximum value of an 4 byte signed integer used in the encoding of integer values.
	 */
	public static final int BYTE4_MAX_VALUE = Integer.MAX_VALUE-5;
	
	/**
	 * A type of integer, encoded as a short (2 byte signed) integer.
	 */
	
	public static final int SHORT_COORDINATE_TYPE = 2;
	
	/**
	 * A type of integer, encoded in a 4 byte signed integer
	 */
	public static final int INTEGER_COORDINATE_TYPE = 4;
	
	/**
	 * A type of integer, that encodes x, y, z coordinate values
	 */
	public static final int ENCODED_COORDINATE_TYPE = 5;

	/**
	 * A 2-byte marker indicating that the next value is of type INTEGER_COORDINATE_TYPE.
	 */
	public static final int BYTE2_INTEGER_MARKER = BYTE2_MAX_VALUE + INTEGER_COORDINATE_TYPE;
	
	/**
	 * A 2-byte marker indicating that the next value is of type ENCODED_COORDINATE_TYPE.
	 */
	public static final int BYTE2_ENCODED_MARKER = BYTE2_MAX_VALUE + ENCODED_COORDINATE_TYPE;

	/**
	 * A 4-byte marker indicating that the next value is of type HORT_COORDINATE_TYPE.
	 */
	public static final int BYTE4_SHORT_MARKER = BYTE4_MAX_VALUE + SHORT_COORDINATE_TYPE;
	
	/**
	 * A 4-byte marker indication that the next value is of type INTEGER_COORDINATE_TYPE
	 */
	public static final int BYTE4_INTEGER_MARKER = BYTE4_MAX_VALUE + INTEGER_COORDINATE_TYPE;
	
	/**
	 * A 4-byte maker indicating that the next values is of type ENCODED_COORDINATE_TYPE
	 */
	public static final int BYTE4_ENCODED_MARKER = BYTE4_MAX_VALUE + ENCODED_COORDINATE_TYPE;

	public static final byte BLANK = 32;
	
	// Structure record identifiers. Record identifiers in lower case have a record length
	// specified by 1-byte values, whereas record identifiers in upper case have
	// the record length specified as integers (4 bytes)
	/**
	 * An identifier for a structure record.
	 */
	public static final byte STRUCTURE = 's';
	
	/**
	 * An identifier for a model record.
	 */
	public static final byte MODEL = 'm';
	
	/** 
	 * An identifier for a chain record.
	 */
	public static final byte CHAIN = 'c';
	
	/**
	 * An identifier for a sequence record.
	 */
	public static final byte SEQUENCE = 'Q';
	
	/**
	 * An identifier for a group (residue) information record.
	 */
	public static final byte GINFO = 'I';
	
	/**
	 * An identifier for a group (residue) record
	 */
	public static final byte GROUP = 'g';

	/**
	 * An identifier of a bond record (reserved for the future).
	 */
	public static final byte BOND = 'B';
	
	/**
	 * An identifier for a coordinate record.
	 */
	public static final byte COORD = 'X';
	
	/**
	 * An identifier for a b factor (temperature factor) record.
	 */
	public static final byte BFACTOR = 'T';
	
	/**
	 * An identifier for a occupancy record.
	 */
	public static final byte OCCUPANCY = 'O';
	
	/**
	 * An identifier for an end record. This record identifier
	 * is the last byte in a .hesc file.
	 */
	public static final byte END = 'e';
}
