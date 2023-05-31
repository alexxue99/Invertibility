package invertibility;

import java.util.Arrays;

/**
 * Class that contains information about the sequences at the
 * leaves of a tree.
 */
public class TreeLeaves {
	private String[] seqLeaves;

	public TreeLeaves(String[] seqLeaves) {
		this.seqLeaves = seqLeaves;
	}

	public int[] getSeqLeavesLengths() {
		return Arrays.stream(seqLeaves).mapToInt(String::length).toArray();	
	}

	private int numOnes(String s) {
		int num = 0;
		for (char c: s.toCharArray())
			if (c == '1')
				num++;
		return num;
	}

	public int[] getSeqNumOnes() {
		return Arrays.stream(seqLeaves).mapToInt(s -> numOnes(s)).toArray();
	}

	public int[] getSeqNumZeros() {
		return Arrays.stream(seqLeaves).mapToInt(s -> s.length() - numOnes(s)).toArray();
	}

	public String[] getSeq() {
		return seqLeaves;
	}

	public int getN() {
		return seqLeaves.length;
	}
}
