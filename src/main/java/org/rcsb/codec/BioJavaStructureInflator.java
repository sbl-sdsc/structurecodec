package org.rcsb.codec;

import java.util.ArrayList;

import org.biojava.nbio.structure.AminoAcidImpl;
import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.AtomImpl;
import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.ChainImpl;
import org.biojava.nbio.structure.Element;
import org.biojava.nbio.structure.Group;
import org.biojava.nbio.structure.HetatomImpl;
import org.biojava.nbio.structure.NucleotideImpl;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureImpl;
import org.biojava.nbio.structure.io.PDBParseException;
import org.rcsb.codec.StructureInflatorInterface;

public class BioJavaStructureInflator implements StructureInflatorInterface {
	private Structure structure;
	private int modelCount = 0;
	private int modelNumber = 0;
	private int atomCount = 1;
	private Chain chain;
	private Group group;

	public BioJavaStructureInflator() {
		this.structure = new StructureImpl();
	}
	
	public Structure getStructure() {
		return this.structure;
	}
	
	@Override
	public void setModelCount(int modelCount) {
//		System.out.println("model count: " + modelCount);
		this.atomCount = 1;
		this.modelCount = modelCount;
	}

	@Override
	public void setModelInfo(int modelNumber, int chainCount) {
//		System.out.println("modelNumber: " + modelNumber + " chainCount: " + chainCount);
		this.modelNumber = modelNumber;
		this.atomCount = 1;
		this.structure.addModel(new ArrayList<Chain>(chainCount));
	}

	@Override
	public void setChainInfo(String chainId, int groupCount) {
//		System.out.println("chainId: " + chainId + " groupCount: " + groupCount);
		chain = new ChainImpl();
		chain.setChainID(chainId.trim());
		structure.addChain(chain, modelNumber);
	}

	@Override
	public void setGroupInfo(String groupName, int groupNumber,
			char insertionCode, int polymerType, int atomCount) {
//		System.out.println("groupName: " + groupName + " groupNumber: " + groupNumber + " insertionCode: " + insertionCode + " polymerType: " + polymerType + " atomCount: " + atomCount);
		switch (polymerType) {
		case 1: 
			group = new AminoAcidImpl();
			break;
		case 2:
			group = new NucleotideImpl();
			break;
		default:
			group = new HetatomImpl();
		}
				

			group.setPDBName(groupName);

		group.setResidueNumber(chain.getChainID().trim(), groupNumber, insertionCode);
		group.setAtoms(new ArrayList<Atom>(atomCount));
	//	System.out.println("getting group: " + group);
		if (polymerType != 0) {
			chain.getSeqResGroups().add(group);
		}
		if (atomCount > 0) {
		
			chain.addGroup(group);	
		}
	}

	@Override
	public void setAtomInfo(String atomName, int serialNumber, char alternativeLocationId, float x,
			float y, float z, float occupancy, float temperatureFactor,
			String element) {
		Atom atom = new AtomImpl();	
		atom.setPDBserial(atomCount++);
//		System.out.println("BioJavaStructureInflator: " + atomCount);
		atom.setName(atomName.trim());
		atom.setElement(Element.valueOfIgnoreCase(element));
		atom.setAltLoc(alternativeLocationId);
		atom.setX(x);
		atom.setY(y);
		atom.setZ(z);
		atom.setOccupancy(occupancy);
		atom.setTempFactor(temperatureFactor);
//		System.out.println(atom);
		group.addAtom(atom);
	}
}
