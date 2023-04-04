package invertibility;

// import static org.junit.Assert.assertEquals;

import java.util.List;

/**
 * Class that contains information about the sequences at the
 * leaves of a tree.
 */
public class TreeLeaves {
	private List<String> seqLeaves;

	public TreeLeaves(List<String> seqLeaves) {
		this.seqLeaves = seqLeaves;
	}

	public List<Integer> getSeqLeavesLengths() {
		return seqLeavesLengths;
	}
}
