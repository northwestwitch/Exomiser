/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package de.charite.compbio.exomiser.core.model.pathogenicity;

import java.util.Collection;

import de.charite.compbio.jannovar.annotation.PutativeImpact;
import de.charite.compbio.jannovar.annotation.VariantEffect;

/**
 *
 * @author Jules Jacobsen <jules.jacobsen@sanger.ac.uk>
 * @author Manuel Holtgrewe <manuel.holtgrewe@charite.de>
 */
public abstract class VariantTypePathogenicityScores {

    /**
     * This is the pathogenicity value we will give to missense (nonsynonymous)
     * variants for which we cannot find values for mutationTaster, polyphen2,
     * or SIFT.
     */
    public static final float DEFAULT_MISSENSE_SCORE = 0.6f;
    /**
     * Pathogenicity score for a mutation class such as INTERGENIC where we are
     * entirely sure it is nonpathogenic (for the purposes of this software).
     */
    public static final float NON_PATHOGENIC_SCORE = 0.00f;
    /**
     * Assumed pathogenicity score for a frameshift mutation.
     */
    public static final float FRAMESHIFT_SCORE = 0.95f;
    /**
     * Assumed pathogenicity score for a nonframeshift indel mutation.
     */
    public static final float NONFRAMESHIFT_INDEL_SCORE = 0.85f;
    /**
     * Assumed pathogenicity score for a nonsense mutation.
     */
    public static final float NONSENSE_SCORE = 0.95f;
    /**
     * Assumed pathogenicity score for a splice site mutation.
     */
    public static final float SPLICING_SCORE = 0.90f;
    /**
     * Assumed pathogenicity score for a synonymous mutation.
     */
    public static final float SYNONYMOUS_SCORE = 0.10f;
    /**
     * Assumed pathogenicity score for a stoploss mutation.
     */
    public static final float STOPLOSS_SCORE = 0.70f;
    /**
     * Assumed pathogenicity score for a variant that causes the start codon to
     * be lost.
     */
    public static final float STARTLOSS_SCORE = 0.95f;

    public static final float getPathogenicityScoreOf(Collection<VariantEffect> effects) {
        if (effects.isEmpty()) {
            // no effect annotated
            return NON_PATHOGENIC_SCORE; 
        }
        // highest priority effect
        VariantEffect variantEffect = effects.iterator().next(); 

        // guard against the case that the highest-impact effect is neither high nor moderate
        if (variantEffect.getImpact().ordinal() > PutativeImpact.MODERATE.ordinal()) {
             // neither HIGH nor MODERATE
            return NON_PATHOGENIC_SCORE;
        }
        switch (variantEffect) {
            case MISSENSE_VARIANT:
                return DEFAULT_MISSENSE_SCORE;
            case SYNONYMOUS_VARIANT:
                return SYNONYMOUS_SCORE;
            case FRAMESHIFT_ELONGATION:
            case FRAMESHIFT_TRUNCATION:
            case FRAMESHIFT_VARIANT:
                return FRAMESHIFT_SCORE;
            case MNV:
            case FEATURE_TRUNCATION:
            case DISRUPTIVE_INFRAME_DELETION:
            case DISRUPTIVE_INFRAME_INSERTION:
            case INFRAME_DELETION:
            case INFRAME_INSERTION:
            case INTERNAL_FEATURE_ELONGATION:
            case COMPLEX_SUBSTITUTION:
                return NONFRAMESHIFT_INDEL_SCORE;
            case SPLICE_ACCEPTOR_VARIANT:
            case SPLICE_DONOR_VARIANT:
            case SPLICE_REGION_VARIANT:
                return SPLICING_SCORE;
            case START_LOST:
                return STARTLOSS_SCORE;
            case STOP_LOST:
                return STOPLOSS_SCORE;
            case STOP_GAINED:
                return NONSENSE_SCORE;
            default:
                return NON_PATHOGENIC_SCORE;
        }
    }

}