package bdmmprime.flow;

import bdmmprime.flow.BirthDeathMigrationDistribution;
import bdmmprime.parameterization.*;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import beast.base.inference.parameter.RealParameter;
import org.apache.commons.math.special.Gamma;
import org.junit.Test;

import static junit.framework.Assert.assertEquals;


/**
 * Created by Jeremie Scire (jscire) on 26.06.17.
 */

public class BirthDeathMigrationLikelihoodTest {

	double runtime;

	/**
     * The original tests were developed assuming BDSKY/BDMM-like behaviour, i.e. return an oriented
	 * tree probability unless r!=1 in which case return an un-oriented and unlabeled tree probability.
	 * In contrast, BDMM-Prime always returns a labeled tree probability.
     *
	 * This method exists to convert BDSKY/BDMM test probabilities to be labeled tree probabilities,
	 * allowing comparison with BDMM-Prime.
	 *
	 * @param density BDMM-prime probability density object
	 * @return conversion factor
	 */
	private double labeledTreeConversionFactor(bdmmprime.flow.BirthDeathMigrationDistribution density) {
		Tree tree = (Tree)density.treeInput.get();
		boolean SAmodel = density.parameterizationInput.get().getRemovalProbs()[0][0] != 1.0;
		double factor = - Gamma.logGamma(tree.getLeafNodeCount() +1);

		if (!SAmodel)
			factor += Math.log(2) * (tree.getLeafNodeCount() - tree.getDirectAncestorNodeCount() - 1);

		return factor;
	}

	/**
	 * Basic test for migration rate change 
	 * Reference from BDMM itself
     * Canonical parameterization
	 */
	@Test 
	public void testLikelihoodMigRateChangeBasicCanonical() {

		// Test for uncoloured tree

        String newick = "(t1[&state=0] : 1.5, t2[&state=1] : 0.5);";

		RealParameter originParam = new RealParameter("2.5");

		Parameterization parameterization = new CanonicalParameterization();
		parameterization.initByName(
		        "typeSet", new TypeSet(2),
                "processLength", originParam,
                "birthRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("2.0"), 2),
                "deathRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("1.0"), 2),
                "birthRateAmongDemes", new SkylineMatrixParameter(
                        null,
                        new RealParameter("0.0"), 2),
                "migrationRate", new SkylineMatrixParameter(
                        new RealParameter("1.0"),
                        new RealParameter("0.1 0.2"), 2),
                "samplingRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("0.5"), 2),
                "removalProb", new SkylineVectorParameter(
                        null,
                        new RealParameter("1.0"), 2));


		bdmmprime.flow.BirthDeathMigrationDistribution density = new bdmmprime.flow.BirthDeathMigrationDistribution();

		density.initByName(
		        "parameterization", parameterization,
                "frequencies", new RealParameter("0.5 0.5"),
                "tree", new TreeParser(newick,
                        false, false,
                        true,0),
                "typeLabel", "state"
		);

		double logL = density.calculateLogP();

		System.out.println("Birth-death result: " + logL + "\t- Test LikelihoodMigRateChange 1");

		// Reference BDMM (version 0.2.0) 22/06/2017
		assertEquals(-6.7022069383966025 - labeledTreeConversionFactor(density),
				logL, 1e-5);
	}
}