package bdmmprime.flow;

import bdmmprime.parameterization.*;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import beast.base.inference.parameter.RealParameter;
import org.apache.commons.math.special.Gamma;
import org.junit.Test;

import static junit.framework.Assert.assertEquals;

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
                "conditionOnSurvival", false,
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
                "conditionOnSurvival", false,
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
                        new RealParameter(String.valueOf(4.0/3.0))),
                "becomeUninfectiousRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("1.5")),
                "R0AmongDemes", new SkylineMatrixParameter(null, null),
                "migrationRate", new SkylineMatrixParameter(null, null),
                "samplingProportion", new SkylineVectorParameter(
                        null,
                        new RealParameter(String.valueOf(1.0/3.0))),
                "removalProb", new SkylineVectorParameter(
                        new RealParameter("1.0"),
                        new RealParameter("0.3 0.7")));

        bdmmprime.flow.BirthDeathMigrationDistribution density = new bdmmprime.flow.BirthDeathMigrationDistribution();
        density.initByName(
                "parameterization", parameterization,
                "frequencies", new RealParameter("1.0"),
                "tree", new TreeParser(newick, false, false, true, 0),
                "conditionOnSurvival", false,
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
                "conditionOnSurvival", false,
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
                "conditionOnSurvival", false,
                "typeLabel", "type"
                );

        double logLnumerical = density.calculateLogP();

        bdmmprime.flow.BirthDeathMigrationDistribution densityExact = new bdmmprime.flow.BirthDeathMigrationDistribution();
        densityExact.initByName(
                "parameterization", parameterization,
                "frequencies", new RealParameter("1.0"),
                "tree", new TreeParser(newick, false, false, true, 0),
                "conditionOnSurvival", false,
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
                        new RealParameter((4.0/3.0) + " 1.1"),
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
                "conditionOnSurvival", false,
                "typeLabel", "type");

        double logL = density.calculateLogP();

        // Reference BDMM (version 0.2.0) 29/03/2018
        assertEquals(-22.82747259570373, logL, 1e-5);
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
                "conditionOnSurvival", false,
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
                "conditionOnSurvival", false,
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
                "conditionOnSurvival", false,
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

    /**
     * Test migration
     * 2 types, migration, no birth among demes
     * Adapted from BDSKY
     */
    @Test
    public void testLikelihoodCalculationMig() {

        // uncoloured tree, asymmetric types
        Tree tree = new TreeParser(
                "((3[&type=0] : 1.5, 4[&type=1] : 0.5) : 1 , (1[&type=1] : 2, 2[&type=0] : 1) : 3);",
                false);

        Parameterization parameterization = new EpiParameterization();
        parameterization.initByName(
                "processLength", new RealParameter("6.0"),
                "typeSet", new TypeSet(2),
                "R0", new SkylineVectorParameter(
                        null,
                        new RealParameter((4.0 / 3.0) + " " + 5.0)),
                "becomeUninfectiousRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("1.5 1.25")),
                "samplingProportion", new SkylineVectorParameter(
                        null,
                        new RealParameter((1.0 / 3.0) + " " + (1.0/2.0))),
                "migrationRate", new SkylineMatrixParameter(
                        null,
                        new RealParameter("0.2 0.1")),

                "removalProb", new SkylineVectorParameter(
                        null,
                        new RealParameter("1.0"), 2));

        bdmmprime.flow.BirthDeathMigrationDistribution density = new bdmmprime.flow.BirthDeathMigrationDistribution();
        density.initByName("parameterization", parameterization,
                "frequencies", new RealParameter("0.5 0.5"),
                "tree", tree,
                "conditionOnSurvival", false,
                "typeLabel", "type"
        );

        assertEquals(-26.53293 + labeledTreeConversionFactor(density), density.calculateLogP(), 1e-5);
    }

    /**
     * Test of migration and infection among demes with rate changes
     * 2 types, no SA
     * Uncoloured tree
     * Reference from BDMM itself (version 0.2.0 28/06/2017)
     * @throws Exception
     */
    @Test
    public void testAmongRateChange() throws Exception {

        Tree tree = new TreeParser("((3[&type=0]:1.5,4[&type=1]:0.5):1,(1[&type=1]:1,2[&type=0]:1):3);",
                false);

        Parameterization parameterization = new EpiParameterization();
        parameterization.initByName(
                "processLength", new RealParameter("4.1"),
                "typeSet", new TypeSet(2),
                "R0", new SkylineVectorParameter(
                        new RealParameter("1.0"),
                        new RealParameter("6 5 2 2.5"), 2),
                "becomeUninfectiousRate", new SkylineVectorParameter(
                        new RealParameter("1.0"),
                        new RealParameter("0.5 0.55 0.45 0.6"), 2),
                "samplingProportion", new SkylineVectorParameter(
                        new RealParameter("1.0"),
                        new RealParameter("0.5 0.45 0.333333 0.35"), 2),
                "R0AmongDemes", new SkylineMatrixParameter(
                        new RealParameter("1.0"),
                        new RealParameter("1.1 1.3 1.2 1.15"), 2),
                "migrationRate", new SkylineMatrixParameter(
                        new RealParameter("1.0"),
                        new RealParameter("0.1 0.15 0.2 0.25")),
                "removalProb", new SkylineVectorParameter(
                        null,
                        new RealParameter("1.0"), 2));

        bdmmprime.flow.BirthDeathMigrationDistribution density = new bdmmprime.flow.BirthDeathMigrationDistribution();
        density.initByName("parameterization", parameterization,
                "frequencies", new RealParameter("0.5 0.5"),
                "tree", tree,
                "conditionOnSurvival", false,
                "typeLabel", "type"
        );

        assertEquals(-16.466832439520886 + labeledTreeConversionFactor(density), density.calculateLogP(), 1e-4); // result from BDMM, 28/06/2017
    }

    /**
     * Test of migration and infection among demes with rate changes
     * 2 types, no SA
     * Uncoloured tree
     * Reference from BDMM itself (version 0.2.0 28/06/2017)
     * @throws Exception
     */
    @Test
    public void testAmongNoRateChange() throws Exception {

        Tree tree = new TreeParser("((3[&type=1]:1.5,4[&type=1]:0.5):1,(1[&type=1]:2,2[&type=1]:1):3);",
                false);

        Parameterization parameterization = new EpiParameterization();
        parameterization.initByName(
                "processLength", new RealParameter("6.0"),
                "typeSet", new TypeSet(2),
                "R0", new SkylineVectorParameter(
                        null,
                        new RealParameter("0 0"), 2),
                "becomeUninfectiousRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("0 0.75"), 2),
                "samplingProportion", new SkylineVectorParameter(
                        null,
                        new RealParameter("0 0.7"), 2),
                "R0AmongDemes", new SkylineMatrixParameter(
                        null,
                        new RealParameter("0 2"), 2),
                "migrationRate", new SkylineMatrixParameter(
                        null,
                        new RealParameter("0.5 0")),
                "removalProb", new SkylineVectorParameter(
                        null,
                        new RealParameter("1.0"), 2));

        bdmmprime.flow.BirthDeathMigrationDistribution density = new bdmmprime.flow.BirthDeathMigrationDistribution();
        density.initByName("parameterization", parameterization,
                "frequencies", new RealParameter("1.0 0.0"),
                "tree", tree,
                "conditionOnSurvival", false,
                "typeLabel", "type"
        );

        assertEquals(-12.1441 + labeledTreeConversionFactor(density), density.calculateLogP(), 1e-4); // tanja's result from R
    }

    /**
     * Test of migration with 3 types
     * No rate change, no SA
     * Reference from BDMM itself
     * @throws Exception
     */
    @Test
    public void testMig3types() throws Exception {

        Tree tree = new TreeParser("((3[&type=2]:1.5,4[&type=1]:0.5):1,(1[&type=1]:1,2[&type=0]:1):3);",
                false);

        Parameterization parameterization = new EpiParameterization();
        parameterization.initByName(
                "processLength", new RealParameter("4.1"),
                "typeSet", new TypeSet(3),
                "R0", new SkylineVectorParameter(
                        null,
                        new RealParameter("6 2 5")),
                "becomeUninfectiousRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("0.5 0.45 0.55")),
                "samplingProportion", new SkylineVectorParameter(
                        null,
                        new RealParameter("0.5 0.333333 0.45")),
                "migrationRate", new SkylineMatrixParameter(
                        null,
                        new RealParameter("0.1 0.2 0.15 0.12 0.12 0.15")),
                "removalProb", new SkylineVectorParameter(
                        null,
                        new RealParameter("1.0"), 3));

        bdmmprime.flow.BirthDeathMigrationDistribution density = new bdmmprime.flow.BirthDeathMigrationDistribution();
        density.initByName("parameterization", parameterization,
                "frequencies", new RealParameter((1.0/3.0) + " " + (1.0/3.0) + " " + (1.0/3.0)),
                "tree", tree,
                "conditionOnSurvival", false,
                "typeLabel", "type"
        );

        assertEquals(-16.88601100061662 + labeledTreeConversionFactor(density), density.calculateLogP(), 1e-4); // result from BDMM, version 0.2.0, 06/07/2017
    }

    /**
     * Likelihood test from the Sasha's SA package.
     */
    @Test
    public void testSALikelihoodMini3() {
        String newick = "((1:1.0,2:0.0):1.0,3:0):0.0";

        Parameterization parameterization = new CanonicalParameterization();
        parameterization.initByName(
                "processLength", new RealParameter("10.0"),
                "birthRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("2.0")),
                "deathRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("0.99")),
                "samplingRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("0.5")),
                "removalProb", new SkylineVectorParameter(
                        null,
                        new RealParameter("0.9")));

        bdmmprime.flow.BirthDeathMigrationDistribution density = new bdmmprime.flow.BirthDeathMigrationDistribution();
        density.initByName(
                "parameterization", parameterization,
                "tree", new TreeParser(newick, false, false, true,0),
                "conditionOnSurvival", false
        );

        // this value was calculated by Sasha with Mathematica
        assertEquals(-25.3707 + labeledTreeConversionFactor(density),
                density.calculateLogP(), 1e-5); // likelihood conditioning on at least one sampled individual
    }

    /**
     * 1-dim and 1 rate-change test
     * reference from BDSKY
     */
    @Test
    public void testLikelihoodRateChangeCondOnSampling1dim() {

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
                        new RealParameter("2.4"),
                        new RealParameter("0.0 0.33333333333")),
                "removalProb", new SkylineVectorParameter(
                        null,
                        new RealParameter("1.0")));

        bdmmprime.flow.BirthDeathMigrationDistribution density = new bdmmprime.flow.BirthDeathMigrationDistribution();
        density.initByName("parameterization", parameterization,
                "frequencies", new RealParameter("1.0"),
                "conditionOnSurvival", true,
                "tree", tree,
                "typeLabel", "state"
        );

        double logPnumeric = density.calculateLogP();
        System.out.println("Numerical solution: " + logPnumeric);

        double logPanalytic = density.calculateLogP();
        System.out.println("Analytical solution: " + logPnumeric);

        assertEquals(logPnumeric, logPanalytic, 1e-4);
    }

    /**
     * Test infection among demes
     * No rate changes
     * Symmetric configuration
     * reference from R
     */
    @Test
    public void testLikelihoodCalculationInfAmongDemesSymmetric() {

        // uncoloured, symmetric tree

        Tree tree = new TreeParser("((t3[&type=1]:0.004214277605,t4[&type=1]:0.02157681391):0.229186993,(t2[&type=0]:0.624713651,t1[&type=1]:1.347400211):0.06231047755);",
                false);

        Parameterization parameterization = new CanonicalParameterization();
        parameterization.initByName(
                "processLength", new RealParameter(Double.toString(tree.getRoot().getHeight() + 0.02686563367)),
                "typeSet", new TypeSet(2),
                "birthRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("2.0"), 2),
                "deathRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("0.5"), 2),
                "samplingRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("0.5"), 2),
                "birthRateAmongDemes", new SkylineMatrixParameter(
                        null,
                        new RealParameter("1.0"), 2),
                "removalProb", new SkylineVectorParameter(
                        null,
                        new RealParameter("1.0"), 2));

        bdmmprime.flow.BirthDeathMigrationDistribution density = new bdmmprime.flow.BirthDeathMigrationDistribution();
        density.initByName("parameterization", parameterization,
                "frequencies", new RealParameter("0.5 0.5"),
                "conditionOnSurvival", true,
                "tree", tree,
                "typeLabel", "type"
        );

        //System.out.println("Log-likelihood " + logL + " - testLikelihoodCalculationInfAmongDemes \t");
        assertEquals(-5.1966118470881 + labeledTreeConversionFactor(density), density.calculateLogP(), 1e-3);

    }

    /**
     * Test infection among demes
     * No rate changes
     * Asymmetric configuration
     * reference from R
     */
    @Test
    public void testLikelihoodCalculationInfAmongDemesAsymmetric() {

        Tree tree = new TreeParser("((3[&type=1]:1.5,4[&type=0]:0.5):1,(1[&type=0]:2,2[&type=1]:1):3);",
                false);

        Parameterization parameterization = new CanonicalParameterization();
        parameterization.initByName(
                "processLength", new RealParameter("6.0"),
                "typeSet", new TypeSet(2),
                "birthRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("2.0 6.25"), 2),
                "deathRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("1.2 0.625"), 2),
                "samplingRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("0.3 0.625"), 2),
                "birthRateAmongDemes", new SkylineMatrixParameter(
                        null,
                        new RealParameter("0.2 0.1"), 2),
                "removalProb", new SkylineVectorParameter(
                        null,
                        new RealParameter("1.0"), 2));

        bdmmprime.flow.BirthDeathMigrationDistribution density = new bdmmprime.flow.BirthDeathMigrationDistribution();
        density.initByName("parameterization", parameterization,
                "frequencies", new RealParameter("0.5 0.5"),
                "conditionOnSurvival", true,
                "tree", tree,
                "typeLabel", "type"
        );

        assertEquals(-26.7939 + labeledTreeConversionFactor(density), density.calculateLogP(), 1e-5);  //result from R
    }

    /**
     * Basic test on sampled-ancestors lik. calculation.
     * 2 leaves, 1 SA. 1 type, no rho-sampling, no rate-change
     * Reference value from BDSKY (23/03/2017)
     */
    @Test
    public void testSALikelihoodMini() {

        Tree tree = new TreeParser("((3[&type=0]: 1.5, 6[&type=0]: 0)5[&type=0]: 3.5, 4[&type=0]: 4) ;",
                false);

        Parameterization parameterization = new EpiParameterization();
        parameterization.initByName(
                "processLength", new RealParameter("6.0"),
                "typeSet", new TypeSet(1),
                "R0", new SkylineVectorParameter(
                        null,
                        new RealParameter("1.5")),
                "becomeUninfectiousRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("1.5")),
                "samplingProportion", new SkylineVectorParameter(
                        null,
                        new RealParameter("0.2")),
                "removalProb", new SkylineVectorParameter(
                        null,
                        new RealParameter("0.9")));

        bdmmprime.flow.BirthDeathMigrationDistribution density = new bdmmprime.flow.BirthDeathMigrationDistribution();
        density.initByName("parameterization", parameterization,
                "frequencies", new RealParameter("1.0"),
                "conditionOnSurvival", true,
                "tree", tree,
                "typeLabel", "type"
        );

        assertEquals(-18.854438107814335 + labeledTreeConversionFactor(density), density.calculateLogP(), 1e-4); //Reference value from BDSKY (23/03/2017)
    }

    /**
     * Likelihood test from the Sasha's SA package.
     */
    @Test
    public void testSALikelihoodMini2() {
        String newick = "((1:1.5,2:0.5):0.5,3:0.0)4:0.0;";

        Parameterization parameterization = new CanonicalParameterization();
        parameterization.initByName(
                "processLength", new RealParameter("10.0"),
                "birthRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("2.0")),
                "deathRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("0.99")),
                "samplingRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("0.5")),
                "removalProb", new SkylineVectorParameter(
                        null,
                        new RealParameter("0.9")));

        bdmmprime.flow.BirthDeathMigrationDistribution density = new bdmmprime.flow.BirthDeathMigrationDistribution();
        density.initByName(
                "parameterization", parameterization,
                "conditionOnSurvival", true,
                "tree", new TreeParser(newick, false, false, true,0)
        );

        // this value was calculated by Sasha with Mathematica
        assertEquals(-22.08332 + labeledTreeConversionFactor(density),
                density.calculateLogP(), 1e-5); // likelihood conditioning on at least one sampled individual
    }

    /**
     * Test on sampled-ancestors lik. calculation with no sampled ancestor
     * No rate-change, one state, 4 tips
     * This state is just there in case something is broken with sampled ancestors,
     * helps for debugging if combined with testSALikelihoodMini for instance
     */
    @Test
    public void testSALikelihoodCalculationWithoutAncestors() {

        Tree tree = new TreeParser("((3[&type=0] : 1.5, 4[&type=0] : 0.5) : 1 , (1[&type=0] : 2, 2[&type=0] : 1) : 3);",
                false);

        Parameterization parameterization = new EpiParameterization();
        parameterization.initByName(
                "processLength", tree,
                "typeSet", new TypeSet(1),
                "R0", new SkylineVectorParameter(
                        null,
                        new RealParameter("1.5")),
                "becomeUninfectiousRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("1.5")),
                "samplingProportion", new SkylineVectorParameter(
                        null,
                        new RealParameter("0.3")),
                "removalProb", new SkylineVectorParameter(
                        null,
                        new RealParameter("0.9")));

        bdmmprime.flow.BirthDeathMigrationDistribution density = new bdmmprime.flow.BirthDeathMigrationDistribution();
        density.initByName("parameterization", parameterization,
                "frequencies", new RealParameter("1.0"),
                "conditionOnSurvival", true,
                "conditionOnRoot", true,
                "tree", tree,
                "typeLabel", "type"
        );

        // Conditioned on root:

        assertEquals(-15.545323363405362 + labeledTreeConversionFactor(density), density.calculateLogP(), 1e-4);

        // Conditioned on origin:

        parameterization.setInputValue("processLength", new RealParameter("10.0"));
        parameterization.initAndValidate();
        density.setInputValue("conditionOnRoot", false);
        density.initAndValidate();

        assertEquals(-25.991511346557598 + labeledTreeConversionFactor(density), density.calculateLogP(), 1e-4);
    }

    /**
     * Tests the case where we have direct ancestors (SA nodes).
     * Tests if two identical trees but with different newick representations have the same likelihood.
     */
    @Test
    public void testDirectAncestor() {

        // two identical trees up to rotation (the two root children are rotated)
        String newick1 = "((1[&type=0]: 1.5, 2[&type=1]: 0.0)3[&type=0]: 3.5, (4[&type=0]: 1.5, 5[&type=1]: 1.5)6[&type=0]: 3.5) ;";
        String newick2 = "((1[&type=0]: 1.5, 2[&type=1]: 1.5)3[&type=0]: 3.5, (4[&type=0]: 1.5, 5[&type=1]: 0.0)6[&type=0]: 3.5) ;";

        Parameterization parameterization = new EpiParameterization();
        parameterization.initByName(
                "typeSet", new TypeSet(2),
                "processLength", new RealParameter("6.0"),
                "R0", new SkylineVectorParameter(
                        null,
                        new RealParameter((4.0/3.0) + " 1.1"),
                        2),
                "becomeUninfectiousRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("1.5 1.4"),
                        2),
                "R0AmongDemes", new SkylineMatrixParameter(
                        null,
                        new RealParameter("0.0"),
                        2),
                "migrationRate", new SkylineMatrixParameter(
                        null,
                        new RealParameter("0.2 0.3"),
                        2),
                "samplingProportion", new SkylineVectorParameter(
                        null,
                        new RealParameter("0.33"),
                        2),
                "removalProb", new SkylineVectorParameter(
                        null,
                        new RealParameter("0.3 0.4")));

        bdmmprime.flow.BirthDeathMigrationDistribution density = new bdmmprime.flow.BirthDeathMigrationDistribution();
        density.initByName(
                "parameterization", parameterization,
                "tree", new TreeParser(newick1, false, false, true,0),
                "frequencies", new RealParameter("0.5 0.5"),
                "conditionOnSurvival", false,
                "typeLabel", "type"
        );

        density.setInputValue("tree", new TreeParser(newick1, false, false, true,0));
        density.initAndValidate();
        double logL1 = density.calculateLogP();

        density.setInputValue("tree", new TreeParser(newick2, false, false, true,0));
        density.initAndValidate();
        double logL2 = density.calculateLogP();

        assertEquals(logL1, logL2, 1e-5);
    }
}