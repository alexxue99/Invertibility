package invertibility;

// import static org.junit.Assert.assertEquals;

import java.util.List;

/**
 * Class that contains information about the lengths of the sequences at the
 * leaves of a tree.
 */
public class TreeLeaves {
	private List<Integer> seqLeavesLengths;

	// Prereq: at least one leaf length must be given.
	public TreeLeaves(List<Integer> seqLeavesLengths) {
		// assertEquals(seqLeavesLengths.size() > 0, true);
		this.seqLeavesLengths = seqLeavesLengths;
	}

	public List<Integer> getSeqLeavesLengths() {
		return seqLeavesLengths;
	}
}
