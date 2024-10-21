package invertibility;

import java.util.Arrays;

/**
 * Class that contains information about the sampled sequences.
 */
public class LeafSamples {
	private String[] samples;

	/**
	 * Constructor to initialize a LeafSamples object.
	 * 
	 * @param samples the sampled sequences
	 */
	public LeafSamples(String[] samples) {
		this.samples = samples;
	}

	/**
	 * @return the lengths of the sequences as an int[]
	 */

	public int[] getSamplesLengths() {
		return Arrays.stream(samples).mapToInt(String::length).toArray();
	}

	private int numOnes(String s) {
		int num = 0;
		for (char c : s.toCharArray())
			if (c == '1')
				num++;
		return num;
	}

	/**
	 * @return the number of ones at the sequences as an int[]
	 */
	public int[] getSeqNumOnes() {
		return Arrays.stream(samples).mapToInt(s -> numOnes(s)).toArray();
	}

	/**
	 * @return the number of zeros at the sequences as an int[]
	 */
	public int[] getSeqNumZeros() {
		return Arrays.stream(samples).mapToInt(s -> s.length() - numOnes(s)).toArray();
	}

	public String[] getSeq() {
		return samples;
	}

	public int getN() {
		return samples.length;
	}
}
