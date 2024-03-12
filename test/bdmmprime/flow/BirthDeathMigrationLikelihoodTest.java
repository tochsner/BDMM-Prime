package bdmmprime.flow;

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
     * <p>
     * This method exists to convert BDSKY/BDMM test probabilities to be labeled tree probabilities,
     * allowing comparison with BDMM-Prime.
     *
     * @param density BDMM-prime probability density object
     * @return conversion factor
     */
    private double labeledTreeConversionFactor(bdmmprime.flow.BirthDeathMigrationDistribution density) {
        Tree tree = (Tree) density.treeInput.get();
        boolean SAmodel = density.parameterizationInput.get().getRemovalProbs()[0][0] != 1.0;
        double factor = -Gamma.logGamma(tree.getLeafNodeCount() + 1);

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
                        true, 0),
                "typeLabel", "state"
        );

        double logL = density.calculateLogP();

        System.out.println("Birth-death result: " + logL + "\t- Test LikelihoodMigRateChange 1");

        // Reference BDMM (version 0.2.0) 22/06/2017
        assertEquals(-6.7022069383966025 - labeledTreeConversionFactor(density),
                logL, 1e-5);
    }

    /**
     * Basic test for migration rate change
     * Reference from BDMM itself
     * Epi parameterization
     */
    @Test
    public void testLikelihoodMigRateChangeBasicEpi() {

        // Test for uncoloured tree

        String newick = "(t1[&state=0] : 1.5, t2[&state=1] : 0.5);";

        RealParameter originParam = new RealParameter("2.5");

        Parameterization parameterization = new EpiParameterization();
        parameterization.initByName(
                "typeSet", new TypeSet(2),
                "processLength", originParam,
                "R0", new SkylineVectorParameter(
                        null,
                        new RealParameter(4.0 / 3.0 + " " + 4.0 / 3.0)),
                "becomeUninfectiousRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("1.5 1.5")),
                "R0AmongDemes", new SkylineMatrixParameter(
                        null,
                        new RealParameter("0.0 0.0")),
                "migrationRate", new SkylineMatrixParameter(
                        new RealParameter("1.0"),
                        new RealParameter("0.1 0.1 0.2 0.2")),
                "samplingProportion", new SkylineVectorParameter(
                        null,
                        new RealParameter(1.0 / 3.0 + " " + 1.0 / 3.0)),
                "removalProb", new SkylineVectorParameter(
                        null,
                        new RealParameter("1.0 1.0")));


        bdmmprime.flow.BirthDeathMigrationDistribution density = new bdmmprime.flow.BirthDeathMigrationDistribution();

        density.initByName(
                "parameterization", parameterization,
                "frequencies", new RealParameter("0.5 0.5"),
                "tree", new TreeParser(newick,
                        false, false,
                        true, 0),
                "typeLabel", "state"
        );

        double logL = density.calculateLogP();

        System.out.println("Birth-death result: " + logL + "\t- Test LikelihoodMigRateChange 1");

        // Reference BDMM (version 0.2.0) 22/06/2017
        assertEquals(-6.7022069383966025, logL, 1e-5);
    }

    /**
     * Basic test for removal-probability rate change
     * Reference from BDMM itself
     *
     * @throws Exception
     */
    @Test
    public void testLikelihoodRemovalProbChangeBasic() {

        String newick = "((1[&type=0]: 1.5, 2[&type=0]: 0)3[&type=0]: 3.5, 4[&type=0]: 4) ;";

        Parameterization parameterization = new EpiParameterization();
        parameterization.initByName(
                "typeSet", new TypeSet(1),
                "processLength", new RealParameter("6.0"),
                "R0", new SkylineVectorParameter(
                        null,
                        new RealParameter(String.valueOf(4.0 / 3.0))),
                "becomeUninfectiousRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("1.5")),
                "R0AmongDemes", new SkylineMatrixParameter(null, null),
                "migrationRate", new SkylineMatrixParameter(null, null),
                "samplingProportion", new SkylineVectorParameter(
                        null,
                        new RealParameter(String.valueOf(1.0 / 3.0))),
                "removalProb", new SkylineVectorParameter(
                        new RealParameter("1.0"),
                        new RealParameter("0.3 0.7")));

        bdmmprime.flow.BirthDeathMigrationDistribution density = new bdmmprime.flow.BirthDeathMigrationDistribution();
        density.initByName(
                "parameterization", parameterization,
                "frequencies", new RealParameter("1.0"),
                "tree", new TreeParser(newick, false, false, true, 0),
                "typeLabel", "type"
        );

        double logL = density.calculateLogP();

        // Reference BDMM (version 0.2.0) 22/06/2017
        assertEquals(-21.25413884159791 + labeledTreeConversionFactor(density),
                logL, 1e-5);

        bdmmprime.flow.BirthDeathMigrationDistribution densityExact = new bdmmprime.flow.BirthDeathMigrationDistribution();
        densityExact.initByName(
                "parameterization", parameterization,
                "frequencies", new RealParameter("1.0"),
                "tree", new TreeParser(newick, false, false, true, 0),
                "typeLabel", "type"
        );

        double logLExact = densityExact.calculateLogP();

        // Reference BDMM (version 0.2.0) 22/06/2017
        assertEquals(-21.25413884159791 + labeledTreeConversionFactor(density),
                logLExact, 1e-5);
    }

    /**
     * Direct comparison between numerical and analytical solutions for a tiny example with no rate changes.
     */
    @Test
    public void tinyAnalyticalTest() {
        String newick = "(1[&type=0]: 1.0, 2[&type=0]: 1.0): 1.0;";

        Parameterization parameterization = new CanonicalParameterization();
        parameterization.initByName(
                "typeSet", new TypeSet(1),
                "processLength", new RealParameter("2.0"),
                "birthRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("2.0")),
                "deathRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("0.5")),
                "birthRateAmongDemes", new SkylineMatrixParameter(null, null),
                "migrationRate", new SkylineMatrixParameter(null, null),
                "samplingRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("0.5")),
                "rhoSampling", new TimedParameter(
                        new RealParameter("2.0"),
                        new RealParameter("0.5")),
                "removalProb", new SkylineVectorParameter(
                        null,
                        new RealParameter("0.5")));

        bdmmprime.flow.BirthDeathMigrationDistribution density = new bdmmprime.flow.BirthDeathMigrationDistribution();
        density.initByName(
                "parameterization", parameterization,
                "frequencies", new RealParameter("1.0"),
                "tree", new TreeParser(newick, false, false, true, 0),
                "typeLabel", "type"
                );

        double logLnumerical = density.calculateLogP();

        bdmmprime.flow.BirthDeathMigrationDistribution densityExact = new bdmmprime.flow.BirthDeathMigrationDistribution();
        densityExact.initByName(
                "parameterization", parameterization,
                "frequencies", new RealParameter("1.0"),
                "tree", new TreeParser(newick, false, false, true, 0),
                "typeLabel", "type"
                );

        double logLanalytical = densityExact.calculateLogP();

        assertEquals(logLnumerical, logLanalytical, 1e-5);
    }

    /**
     * Two-state test for removal-probability rate change
     * Reference from BDMM itself
     */
    @Test
    public void testLikelihoodRemovalProbChangeTwoState() {

        String newick = "((1[&type=0]: 1.5, 2[&type=1]: 0)3[&type=0]: 3.5, 4[&type=1]: 4) ;";

        Parameterization parameterization = new EpiParameterization();
        parameterization.initByName(
                "typeSet", new TypeSet(2),
                "processLength", new RealParameter("6.0"),
                "R0", new SkylineVectorParameter(
                        new RealParameter("1.0"),
                        new RealParameter((4.0 / 3.0) + " 1.1"),
                        2),
                "becomeUninfectiousRate", new SkylineVectorParameter(
                        new RealParameter("1.0"),
                        new RealParameter("1.5 1.4"),
                        2),
                "R0AmongDemes", new SkylineMatrixParameter(
                        null,
                        new RealParameter("0.0"),
                        2),
                "migrationRate", new SkylineMatrixParameter(
                        new RealParameter("1.0"),
                        new RealParameter("0.2 0.3"),
                        2),
                "samplingProportion", new SkylineVectorParameter(
                        null,
                        new RealParameter("0.33"),
                        2),
                "removalProb", new SkylineVectorParameter(
                        new RealParameter("1.0"),
                        new RealParameter("0.3 0.4 0.7 0.6")));

        bdmmprime.flow.BirthDeathMigrationDistribution density = new bdmmprime.flow.BirthDeathMigrationDistribution();
        density.initByName(
                "parameterization", parameterization,
                "frequencies", new RealParameter("0.5 0.5"),
                "tree", new TreeParser(newick, false, false, true, 0),
                "typeLabel", "type");

        double logL = density.calculateLogP();

        // Reference BDMM (version 0.2.0) 29/03/2018
        assertEquals(-21.185194919464568 + labeledTreeConversionFactor(density), logL, 1e-5);
    }

    /**
     * Basic 1-dim test
     * No rate change, 1 state, no rho-sampling
     * Reference from BDSKY
     */
    @Test
    public void testLikelihood1dim() {

        Tree tree = new TreeParser( "((3[&state=0] : 1.5, 4[&state=0] : 0.5)[&state=0] : 1 , (1[&state=0] : 2, 2[&state=0] : 1)[&state=0] : 3)[&state=0];",
                false);

        Parameterization parameterization = new EpiParameterization();
        parameterization.initByName(
                "processLength", new RealParameter("6.0"),
                "typeSet", new TypeSet(1),
                "R0", new SkylineVectorParameter(
                        null,
                        new RealParameter("1.3333333334")),
                "becomeUninfectiousRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("1.5")),
                "samplingProportion", new SkylineVectorParameter(
                        null,
                        new RealParameter("0.33333333333")),
                "removalProb", new SkylineVectorParameter(
                        null,
                        new RealParameter("1.0")));

        bdmmprime.flow.BirthDeathMigrationDistribution density = new bdmmprime.flow.BirthDeathMigrationDistribution();
        density.initByName("parameterization", parameterization,
                "frequencies", new RealParameter("1.0"),
                "tree", tree,
                "typeLabel", "state"
                );

        assertEquals(-19.019796073623493 + labeledTreeConversionFactor(density), density.calculateLogP(), 1e-5);   // Reference BDSKY (version 1.3.3)
    }

    /**
     * 1-dim and 1 rate-change test
     * reference from BDSKY
     */
    @Test
    public void testLikelihoodRateChange1dim() {

        Tree tree = new TreeParser("((3[&state=0] : 1.5, 4[&state=0] : 0.5)[&state=0] : 1 , (1[&state=0] : 2, 2[&state=0] : 1)[&state=0] : 3)[&state=0];",
                false);

        Parameterization parameterization = new EpiParameterization();
        parameterization.initByName(
                "processLength", new RealParameter("6.0"),
                "typeSet", new TypeSet(1),
                "R0", new SkylineVectorParameter(
                        new RealParameter("3.0"),
                        new RealParameter("0.6666666667 1.3333333334")),
                "becomeUninfectiousRate", new SkylineVectorParameter(
                        new RealParameter("3.0"),
                        new RealParameter("4.5 1.5")),
                "samplingProportion", new SkylineVectorParameter(
                        new RealParameter("3.0"),
                        new RealParameter("0.4444444444 0.33333333333")),
                "removalProb", new SkylineVectorParameter(
                        null,
                        new RealParameter("1.0")));

        bdmmprime.flow.BirthDeathMigrationDistribution density = new bdmmprime.flow.BirthDeathMigrationDistribution();
        density.initByName("parameterization", parameterization,
                "frequencies", new RealParameter("1.0"),
                "tree", tree,
                "typeLabel", "state");

        assertEquals(-33.7573 + labeledTreeConversionFactor(density), density.calculateLogP(), 1e-4); // Reference BDSKY
    }

    /**
     * Basic tests on 2 types situations with migration or birth among demes
     * reference from R
     * @throws Exception
     */
    @Test
    public void testLikelihoodCalculationMigTiny() throws Exception {

        // migration and no birth among demes

        Tree tree = new TreeParser("(1[&state=0] : 1.5, 2[&state=1] : 0.5)[&state=0];", false);
        Parameterization parameterization = new EpiParameterization();
        parameterization.initByName(
                "processLength", new RealParameter("2.5"),
                "typeSet", new TypeSet(2),
                "R0", new SkylineVectorParameter(
                        null,
                        new RealParameter(Double.toString(4.0/3.0)), 2),
                "becomeUninfectiousRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("1.5"), 2),
                "samplingProportion", new SkylineVectorParameter(
                        null,
                        new RealParameter(Double.toString(1.0/3.0)), 2),
                "migrationRate", new SkylineMatrixParameter(
                        null,
                        new RealParameter("0.1"), 2),

                "removalProb", new SkylineVectorParameter(
                        null,
                        new RealParameter("1.0"), 2));

        bdmmprime.flow.BirthDeathMigrationDistribution density = new bdmmprime.flow.BirthDeathMigrationDistribution();
        density.initByName("parameterization", parameterization,
                "frequencies", new RealParameter("0.5 0.5"),
                "tree", tree,
                "typeLabel", "state"
                );

        assertEquals(-7.215222 + labeledTreeConversionFactor(density), density.calculateLogP(), 1e-6); // result from R

        // no migration, symmetric birth among demes

        parameterization.setInputValue("migrationRate", null);
        parameterization.setInputValue("R0AmongDemes", new SkylineMatrixParameter(
                null,
                new RealParameter("0.0666667"), 2));
        parameterization.initAndValidate();
        density.initAndValidate();

        assertEquals(-7.404888 + labeledTreeConversionFactor(density), density.calculateLogP(), 1e-6); // result from R

        // no migration, asymmetric birth among demes

        parameterization.setInputValue("R0AmongDemes", new SkylineMatrixParameter(
                null,
                new RealParameter("0.0666667 0.1"), 2));
        parameterization.initAndValidate();
        density.initAndValidate();

        assertEquals(-7.18723 + labeledTreeConversionFactor(density), density.calculateLogP(), 1e-6); // result from R


        // no migration, asymmetric R0, asymmetric birth among demes

        parameterization.setInputValue("R0", new SkylineVectorParameter(
                null,
                new RealParameter("2 1.3333333")));
        parameterization.initAndValidate();
        density.initAndValidate();

        assertEquals(-7.350649 + labeledTreeConversionFactor(density), density.calculateLogP(), 1e-6); // result from R

        // no migration, asymmetric R0, birth among demes, BU rate, samp proportion

        parameterization.setInputValue("R0", new SkylineVectorParameter(
                null,
                new RealParameter("2.0 1.5")));
        parameterization.setInputValue("becomeUninfectiousRate", new SkylineVectorParameter(
                null,
                new RealParameter("2.0 1.0")));
        parameterization.setInputValue("samplingProportion", new SkylineVectorParameter(
                null,
                new RealParameter("0.5 0.3")));
        parameterization.setInputValue("R0AmongDemes", new SkylineMatrixParameter(
                null,
                new RealParameter("0.1 0.5")));
        parameterization.initAndValidate();
        density.initAndValidate();

        assertEquals(-6.504139 + labeledTreeConversionFactor(density), density.calculateLogP(), 1e-6); // result from R

        // Same params as last test, swapped leaf states

        tree = new TreeParser("(1[&state=1] : 1.5, 2[&state=0] : 0.5)[&state=0];", false);
        density.setInputValue("tree", tree);
        density.initAndValidate();

        assertEquals(-7.700916 + labeledTreeConversionFactor(density), density.calculateLogP(), 1e-6); // result from R
    }
}