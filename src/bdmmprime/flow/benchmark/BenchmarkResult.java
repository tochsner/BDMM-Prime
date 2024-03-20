package bdmmprime.flow.benchmark;

import bdmmprime.parameterization.Parameterization;
import beast.base.evolution.tree.Tree;

import java.util.StringJoiner;


public class BenchmarkResult {

    Parameterization parameterization;
    Tree tree;
    BenchmarkRun flowRun;
    BenchmarkRun bdmmRun;

    public BenchmarkResult(Parameterization parameterization, Tree tree, BenchmarkRun flowRun, BenchmarkRun bdmmRun) {
        this.parameterization = parameterization;
        this.tree = tree;
        this.flowRun = flowRun;
        this.bdmmRun = bdmmRun;
    }

    @Override
    public String toString() {
        StringJoiner joiner = new StringJoiner(",");

        joiner.add(Integer.toString(this.tree.getNodeCount()));
        joiner.add(Integer.toString(this.tree.getLeafNodeCount()));

        joiner.add(Integer.toString(this.parameterization.getNTypes()));

        joiner.add(Double.toString(this.flowRun.likelihood));
        joiner.add(Long.toString(this.flowRun.duration));

        joiner.add(Double.toString(this.bdmmRun.likelihood));
        joiner.add(Long.toString(this.bdmmRun.duration));

        return joiner.toString();
    }

    public static String getHeaders() {
        StringJoiner joiner = new StringJoiner(",");

        joiner.add("node_count");
        joiner.add("leaf_count");

        joiner.add("types_count");

        joiner.add("flow_likelihood");
        joiner.add("flow_duration");

        joiner.add("bdmm_likelihood");
        joiner.add("bdmm_duration");

        return joiner.toString();
    }
}
