package bdmmprime.util;

import beast.core.Input;
import beast.evolution.tree.TraitSet;

import java.util.stream.Collectors;

/**
 * TraitSet in which value input is optional and the values
 * are initialized to a stand-in value.  Used by the BDMM-Prime BEAUti template,
 * where the trait set must be specified before any value (or indeed
 * the taxa themselves) can be known.
 */
public class InitializedTraitSet extends TraitSet {

    public InitializedTraitSet() {
        traitsInput.setRule(Input.Validate.OPTIONAL);
    }

    @Override
    public void initAndValidate() {

        if (traitsInput.get() == null) {
            String value = taxaInput.get().getTaxaNames().stream()
                    .map(n -> n + "=NOT_SET")
                    .collect(Collectors.joining(","));

            traitsInput.setValue(value, this);
        }

        super.initAndValidate();
    }
}
