/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package de.charite.compbio.exomiser.core.util;

import de.charite.compbio.exomiser.core.model.Variant;
import de.charite.compbio.exomiser.core.model.Gene;
import de.charite.compbio.exomiser.core.model.VariantEvaluation;
import de.charite.compbio.jannovar.pedigree.CompatibilityCheckerException;
import de.charite.compbio.jannovar.pedigree.Genotype;
import de.charite.compbio.jannovar.pedigree.GenotypeListBuilder;
import de.charite.compbio.jannovar.pedigree.ModeOfInheritance;
import de.charite.compbio.jannovar.pedigree.Pedigree;
import de.charite.compbio.jannovar.pedigree.PedigreeDiseaseCompatibilityDecorator;
import de.charite.compbio.jannovar.reference.GenomeInterval;
import htsjdk.variant.variantcontext.Allele;

import java.util.ArrayList;
import java.util.EnumSet;
import java.util.List;
import java.util.Set;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import com.google.common.collect.ImmutableList;
import htsjdk.variant.variantcontext.VariantContext;

/**
 *
 * This class allows us to do segregation analysis for the variants supplied to
 * it i.e., to determine if they are compatible with autosomal recessive,
 * autosomal dominant, or X-linked recessive inheritance.
 *
 * @author Jules Jacobsen <jules.jacobsen@sanger.ac.uk>
 */
public class InheritanceModeAnalyser {

    private static final Logger logger = LoggerFactory.getLogger(InheritanceModeAnalyser.class);

    private final Pedigree pedigree;

    public InheritanceModeAnalyser(Pedigree pedigree) {
        this.pedigree = pedigree;
    }

    /**
     * Analyses the inheritance modes for a gene according to the 
     * @param gene
     * @return a Set of inheritance modes with which the gene is compatible with.
     */
    public Set<ModeOfInheritance> analyseInheritanceModesForGene(Gene gene) {
        return analyseInheritanceModes(gene.getPassedVariantEvaluations());
    }

    /**
     * Caution - this only really works if 
     * @param variantEvaluations
     * @return 
     */
    protected Set<ModeOfInheritance> analyseInheritanceModes(List<VariantEvaluation> variantEvaluations) {
        Set<ModeOfInheritance> inheritanceModes = EnumSet.noneOf(ModeOfInheritance.class);
        
        List<Variant> variantList = new ArrayList<>();
        
        for (VariantEvaluation variantEvaluation : variantEvaluations) {
            variantList.add(variantEvaluation.getVariant());
        }
        inheritanceModes.addAll(analyseInheritanceModesForVariants(variantList));
      
        return inheritanceModes;
    }

    private Set<ModeOfInheritance> analyseInheritanceModesForVariants(List<Variant> variants) {
        Set<ModeOfInheritance> inheritanceModes = EnumSet.noneOf(ModeOfInheritance.class);
           
        Variant firstVariant = variants.get(0);
        // Build list of genotypes from the given variants.
        String geneID = firstVariant.getGeneSymbol();
        // Use interval of transcript of first region, only used for the chromosome information anyway.
        GenomeInterval geneInterval = variants.get(0).getAnnotationList().get(0).transcript.txRegion;
        GenotypeListBuilder genotypeListBuilder = new GenotypeListBuilder(geneID, pedigree.getNames(), firstVariant.isXChromosomal());
        for (Variant variant : variants) {
            final int altAlleleID = variant.getAltAlleleID();
            VariantContext variantContext = variant.getVariantContext();
            final int numSamples = variantContext.getNSamples();
            ImmutableList.Builder<Genotype> gtBuilder = new ImmutableList.Builder<Genotype>();
            for (int i = 0; i < numSamples; ++i) {
                final String name = pedigree.members.get(i).name;
                final List<Allele> alleles = variantContext.getGenotype(name).getAlleles();
                if (alleles.size() != 2) {
                    gtBuilder.add(Genotype.NOT_OBSERVED);
                    continue;
                }
                
                final boolean isAlt0 = alleles.get(0).basesMatch(variantContext.getAlternateAllele(altAlleleID));
                final boolean isAlt1 = alleles.get(1).basesMatch(variantContext.getAlternateAllele(altAlleleID));
                if (!isAlt0 && !isAlt1) {
                    gtBuilder.add(Genotype.HOMOZYGOUS_REF);
                } else if ((isAlt0 && !isAlt1) || (!isAlt0 && isAlt1)) {
                    gtBuilder.add(Genotype.HETEROZYGOUS);
                } else {
                    gtBuilder.add(Genotype.HOMOZYGOUS_ALT);
                }
            }
            genotypeListBuilder.addGenotypes(gtBuilder.build());
        }

        // FIXME: there are more modes of inheritance implemented in Jannovar
        final ImmutableList<ModeOfInheritance> toCheck = ImmutableList.of(
                ModeOfInheritance.AUTOSOMAL_RECESSIVE,
                ModeOfInheritance.AUTOSOMAL_DOMINANT, 
                ModeOfInheritance.X_RECESSIVE);
        
        PedigreeDiseaseCompatibilityDecorator checker = new PedigreeDiseaseCompatibilityDecorator(pedigree);
        for (ModeOfInheritance mode : toCheck) {
            try {
                if (checker.isCompatibleWith(genotypeListBuilder.build(), mode)) {
                    inheritanceModes.add(mode);
                }
            } catch (CompatibilityCheckerException e) {
                throw new RuntimeException("Problem in the mode of inheritance checks!", e);
            }
        }

        return inheritanceModes;
    }

}
