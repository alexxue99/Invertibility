package invertibility;

import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.Random;
import java.util.stream.Collectors;

/** Class used to simulate the length process, given a tree. */
public class TreeSimul {
	private double lambda;
	private double mu;
	private double nu;
	private double pi0;
	private String root;
	private int rootVertex;
	private int M;
	private HashMap<Integer, LinkedList<TimeToVertex>> adj; // (v2, t) is in the list for key v1 if (v1, v2) is a
															// directed
	// edge with time t
	private HashMap<Integer, String> seq; // maps vertices v to their strings
	private LinkedList<Leaf> seqLeaves; // contains information about the sequences at the leaves

	/**
	 * Class that represents a vertex and the time corresponding to a directed edge ending at the
	 * vertex.
	 */
	private static class TimeToVertex {
		private int v2;
		private double time;

		public TimeToVertex(int v2, double time) {
			this.v2 = v2;
			this.time = time;
		}

		public int getVertex() {
			return v2;
		}

		public double getTime() {
			return time;
		}
	}

	/** Class that represents a directed edge. */
	public static class Edge {
		private int v1;
		private int v2;

		public Edge(int v1, int v2) {
			this.v1 = v1;
			this.v2 = v2;
		}

		public int getVertex1() {
			return v1;
		}

		public int getVertex2() {
			return v2;
		}
	}

	/**
	 * Class that represents a leaf.
	 * 
	 * @param v
	 *            is the vertex.
	 * @param seq
	 *            is the sequence associated with the vertex.
	 */
	public static class Leaf {
		private int v;
		private String seq;

		public Leaf(int v, String seq) {
			this.v = v;
			this.seq = seq;
		}

		public int getVertex() {
			return v;
		}

		public String getSequence() {
			return seq;
		}
	}

	/**
	 * Constructor to initialize tree.
	 * 
	 * @param lambda
	 *            is insertion rate.
	 * @param mu
	 *            is deletion rate.
	 * @param nu
	 *            is substitution rate.
	 * @param pi0
	 *            is probability that a character is 0.
	 * @param rootVertex
	 *            is the integer corresponding to the root vertex.
	 * @param root
	 *            is the string at the root.
	 * @param edges
	 *            is the list of directed edges.
	 * @param t
	 *            contains the times, where t[i] corresponds to the time for the ith
	 *            edge.
	 */
	public TreeSimul(double lambda, double mu, double nu, double pi0, int rootVertex, String root,
			LinkedList<Edge> edges, double[] t) {
		assert lambda >= 0 && mu >= 0 && nu >= 0 && pi0 >= 0 && pi0 <= 1;
		this.lambda = lambda;
		this.mu = mu;
		this.nu = nu;
		this.pi0 = pi0;
		this.root = root;
		this.rootVertex = rootVertex;
		M = root.length();

		// set up the adjacency map
		adj = new HashMap<>();
		Iterator<Edge> edgeIter = edges.iterator();
		while (edgeIter.hasNext()) {
			int v1 = edgeIter.next().getVertex1();
			adj.put(v1, new LinkedList<TimeToVertex>());
		}

		edgeIter = edges.iterator();
		int i = 0;
		while (edgeIter.hasNext()) {
			Edge curEdge = edgeIter.next();
			int v1 = curEdge.getVertex1();
			int v2 = curEdge.getVertex2();
			adj.get(v1).add(new TimeToVertex(v2, t[i++]));
		}

		seq = new HashMap<Integer, String>();
		seq.put(rootVertex, root);

		seqLeaves = new LinkedList<Leaf>();

		runLengthProcess();
	}

	/** Returns a variable with exponential distribution with parameter var. */
	private double getExp(double var) {
		Random rand = new Random();
		return Math.log(1 - rand.nextDouble()) / (-var);
	}

	/** Returns a 0 with probability pi0 and a 1 with probability 1 - pi0. */
	private char getChar() {
		Random rand = new Random();
		return (rand.nextDouble() <= pi0) ? '0' : '1';
	}

	/**
	 * Evolves string s along an edge for time t, using the length process
	 * parameters.
	 */
	private String evolve(String s, double t) {
		while (t > 0) {
			double[] timings = new double[s.length() * 3];

			// get exponential variables under parameters lambda, mu, and nu, and then find
			// the smallest variable
			for (int i = 0; i < s.length(); i++) {
				timings[3 * i] = getExp(lambda);
				timings[3 * i + 1] = getExp(mu);
				timings[3 * i + 2] = getExp(nu);
			}

			double min = Double.MAX_VALUE;
			int minIndex = 0;

			for (int i = 0; i < timings.length; i++) {
				if (timings[i] < min) {
					minIndex = i;
					min = timings[i];
				}
			}

			t -= min;
			if (t < 0)
				break;

			int location = minIndex / 3;
			switch (minIndex % 3) {
			case 0: // insertion
				s = s.substring(0, location + 1) + getChar() + s.substring(location + 1);
				break;
			case 1: // deletion
				s = s.substring(0, location) + s.substring(location + 1);
				break;
			case 2: // substitution
				s = s.substring(0, location) + getChar() + s.substring(location + 1);
				break;
			}
		}

		return s;
	}

	/** Simulates the length process starting from vertex curVertex. */
	private void runFromVertex(int curVertex) {
		Iterator<TimeToVertex> iter = adj.getOrDefault(curVertex, new LinkedList<TimeToVertex>()).iterator();

		if (!iter.hasNext()) // curVertex is a leaf 
			seqLeaves.add(new Leaf(curVertex, seq.get(curVertex)));

		while (iter.hasNext()) {
			TimeToVertex next = iter.next();
			seq.put(next.getVertex(), evolve(seq.get(curVertex), next.getTime()));
			runFromVertex(next.getVertex());
		}
	}

	/** Simulates the length process on the whole tree. */
	public void runLengthProcess() {
		runFromVertex(rootVertex);
	}

	/**
	 * Transforms a TreeSimul object into a TreeLeaves object.
	 */
	public TreeLeaves toTreeLeaves() {
		return new TreeLeaves(seqLeaves.stream().map(leaf -> leaf.getSequence().length()).collect(Collectors.toList()));
	}

	public double getLambda() {
		return lambda;
	}

	public double getMu() {
		return mu;
	}

	public double getNu() {
		return nu;
	}

	public double getPi0() {
		return pi0;
	}

	public String getRoot() {
		return root;
	}

	public int getM() {
		return M;
	}

	public LinkedList<Leaf> getSeqLeaves() {
		return seqLeaves;
	}
}
