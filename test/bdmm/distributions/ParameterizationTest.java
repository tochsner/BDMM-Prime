package bdmm.distributions;

import beast.core.parameter.RealParameter;
import org.junit.Assert;
import org.junit.Test;

public class ParameterizationTest {

    public double TOLERANCE = 1e-20;

    @Test
    public void basicTest() {

		RealParameter originParam = new RealParameter("2.0");

		Parameterization parameterization = new CanonicalParameterization();
		parameterization.initByName(
		        "nTypes", 2,
                "origin", originParam,
                "birthRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("4.0 4.0")),
                "deathRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("3.0 3.0")),
                "birthRateAmongDemes", new SkylineMatrixParameter(
                        null,
                        new RealParameter("0.0 0.0")),
                "migrationRate", new SkylineMatrixParameter(
                        new RealParameter("1.0"),
                        new RealParameter("0.1 0.1 0.2 0.2")),
                "samplingRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("1.5 1.5")),
                "removalProb", new SkylineVectorParameter(
                        null,
                        new RealParameter("1.0 1.0")),
                "rhoSampling", new TimedParameter(
                        originParam,
                        new RealParameter("0.0 0.0")));

		System.out.println("Number of intervals: " + parameterization.getTotalIntervalCount());

        Assert.assertEquals(3, parameterization.birthRates.length);
        Assert.assertEquals(3, parameterization.deathRates.length);
        Assert.assertEquals(3, parameterization.samplingRates.length);
        Assert.assertEquals(3, parameterization.removalProbs.length);
        Assert.assertEquals(3, parameterization.rhoValues.length);

        for (int interval=0; interval<3; interval++) {
            double migRate = interval < 1 ? 0.1 : 0.2;

            for (int state1 = 0; state1 < 2; state1++) {
                Assert.assertEquals(4.0, parameterization.birthRates[interval][state1], TOLERANCE);
                Assert.assertEquals(3.0, parameterization.deathRates[interval][state1], TOLERANCE);
                Assert.assertEquals(1.5, parameterization.samplingRates[interval][state1], TOLERANCE);
                Assert.assertEquals(1.0, parameterization.removalProbs[interval][state1], TOLERANCE);
                Assert.assertEquals(0.0, parameterization.rhoValues[interval][state1], TOLERANCE);

                for (int state2 = 0; state2 < 2; state2++) {
                    if (state2 == state1)
                        continue;

                    Assert.assertEquals(migRate, parameterization.migRates[interval][state1][state2], TOLERANCE);
                    Assert.assertEquals(0.0, parameterization.crossBirthRates[interval][state1][state2], TOLERANCE);
                }
            }
        }
    }
}
