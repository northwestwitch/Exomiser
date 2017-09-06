/*
 * The Exomiser - A tool to annotate and prioritize genomic variants
 *
 * Copyright (c) 2016-2017 Queen Mary University of London.
 * Copyright (c) 2012-2016 Charité Universitätsmedizin Berlin and Genome Research Ltd.
 * Copyright (c) 2017 Berlin Institute of Health.
 *
 * This program is free software: you can redistribute it and/or modify it under the terms of the
 * GNU Affero General Public License as published by the Free Software Foundation, either version 3
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
 * even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License along with this program.
 * If not, see <http://www.gnu.org/licenses/>.
 */

package org.monarchinitiative.exomiser.core.prioritisers;

import static java.util.stream.Collectors.toMap;

import com.github.phenomics.ontolib.base.OntoLibException;
import com.github.phenomics.ontolib.formats.hpo.HpoOntology;
import com.github.phenomics.ontolib.formats.hpo.HpoTerm;
import com.github.phenomics.ontolib.formats.hpo.HpoTermRelation;
import com.github.phenomics.ontolib.io.base.TermAnnotationParserException;
import com.github.phenomics.ontolib.io.obo.hpo.HpoGeneAnnotationParser;
import com.github.phenomics.ontolib.io.obo.hpo.HpoOboParser;
import com.github.phenomics.ontolib.io.scoredist.H2ScoreDistributionReader;
import com.github.phenomics.ontolib.ontology.algo.InformationContentComputation;
import com.github.phenomics.ontolib.ontology.data.ImmutableOntology;
import com.github.phenomics.ontolib.ontology.data.ImmutableTermId;
import com.github.phenomics.ontolib.ontology.data.Ontology;
import com.github.phenomics.ontolib.ontology.data.Term;
import com.github.phenomics.ontolib.ontology.data.TermAnnotation;
import com.github.phenomics.ontolib.ontology.data.TermAnnotations;
import com.github.phenomics.ontolib.ontology.data.TermId;
import com.github.phenomics.ontolib.ontology.data.TermRelation;
import com.github.phenomics.ontolib.ontology.scoredist.ObjectScoreDistribution;
import com.github.phenomics.ontolib.ontology.similarity.PrecomputingPairwiseResnikSimilarity;
import com.github.phenomics.ontolib.ontology.similarity.ResnikSimilarity;
import java.io.File;
import java.io.IOException;
import java.util.AbstractMap;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.Stream;
import org.monarchinitiative.exomiser.core.model.Gene;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Variant priotization with the Phenix algorithm.
 *
 * <p>
 * Variants are scored using the Phenomizer score which uses Resnik similarity in the HPO phenotypic
 * abnormality sub ontology and (optionally) p-value computations based on empirically observed
 * Resnik similarity score distributions.
 * </p>
 *
 * <p>
 * For the Phenomizer score, the <tt>hp.obo</tt>file is required along with the
 * <tt>ALL_SOURCES_ALL_FREQUENCIES_genes_to_phenotype.txt</tt> file. These two files are available
 * through the Human Phenotype Ontology Jenkins server. The file <tt>phenix.h2.mv.db</tt> is also
 * required and distributed with the Exomiser data ZIP file.
 * </p>
 *
 * <h4>Notes and Open Points</h4>
 *
 * <p>
 * Note that the genes will only be ranked by the normalized score. In the future, we should also
 * try to rank by p-value.
 * </p>
 *
 * @author Sebastian Koehler
 *         <a href= "dr.sebastian.koehler@gmail.com">dr.sebastian.koehler@gmail.com</a>
 * @author Manuel Holtgrewe <a href= "manuel.holtgrewe@bihealth.de">manuel.holtgrewe@bihealth.de</a>
 */
public final class PhenixPriority implements Prioritiser {

    /** Logger instance to use in this class. */
    private static final Logger LOGGER = LoggerFactory.getLogger(PhenixPriority.class);

    /**
     * {@link PriorityType} to label with.
     */
    private static final PriorityType PRIORITY_TYPE = PriorityType.PHENIX_PRIORITY;

    /**
     * The HPO as an {@link ImmutableOntology} object.
     */
    private HpoOntology hpo;

    /** The phenotypic abnormality sub ontology. */
    private Ontology<HpoTerm, HpoTermRelation> abnormalPhenoSubOntology;

    /** The semantic similarity measure used to calculate phenotypic similarity. */
    private ResnikSimilarity<HpoTerm, HpoTermRelation> resnikSimilarity;

    /** Labels of Entrez gene ID to positively associated terms. */
    private Map<Integer, Collection<TermId>> entrezGeneIdToTermIds;

    /** Reader for Ontolib score distributions. */
    private H2ScoreDistributionReader scoreDistReader;

    /** The score to use if computing a score fails (missing link). */
    private static final double DEFAULT_SCORE = 0.0;

    /**
     * The negative log of the p-value to use if computing a score fails (missing link).
     */
    private static final double DEFAULT_NEGATIVE_LOG_P = 0.0;

    /** Links between genes and phenotypes. */
    private List<TermAnnotation> termToGeneAnnotations;

    /** Whether or not to use symmetric score variant. */
    private boolean symmetric;

    /**
     * Path to the directory that has the files needed to calculate the score distribution.
     */
    private String pathPhenixData;

    /**
     * Create a new instance of the PhenixPriority.
     *
     * @param pathPhenixData Path to directory with Phenix data ({@code hp.obo},
     *        {@code ALL_SOURCES_ALL_FREQUENCIES_genes_to_phenotype.txt}, and
     *        {@code phenix.h2.mv.db} file.
     * @param symmetric Flag to indicate if the semantic similarity score should be calculated using
     *        the symmetric formula.
     */
    public PhenixPriority(String pathPhenixData, boolean symmetric) {
        LOGGER.info("Constructing Phenix Priotizer...");
        if (!pathPhenixData.endsWith(File.separator)) {
            pathPhenixData += File.separator;
        }
        this.pathPhenixData = pathPhenixData;
        this.symmetric = symmetric;

        final String pathPhenixDb = pathPhenixData + "phenix";
        try {
            this.scoreDistReader =
                    new H2ScoreDistributionReader(pathPhenixDb, "PHENIX_SCORE_DISTRIBUTION");
        } catch (OntoLibException e) {
            throw new PhenixException(
                    "Could not load phenix score distribution database " + pathPhenixDb, e);
        }

        this.loadHpoData();
        this.buildResnikSimilarity();
        LOGGER.info("Done constructing Phenix Priotizer.");
    }

    /**
     * Stub constructor, used for testing only.
     * 
     * @param symmetric
     */
    protected PhenixPriority(boolean symmetric) {
        this.symmetric = symmetric;
    }

    /**
     * Load the HPO data from <tt>hp.obo</tt> and
     * <tt>ALL_SOURCES_ALL_FREQUENCIES_genes_to_phenotype.txt</tt>.
     *
     * <p>
     * Will set {@link #hpo}, {@link #abnormalPhenoSubOntology}, and {@link #entrezGeneIdToTerms} as
     * a side effect. These are later used by {@link #buildResnikSimilarity()} for building the
     * Resnik similarity computation.
     * </p>
     */
    private void loadHpoData() {
        // Load HPO data from OBO file.
        final String pathHpoObo = this.pathPhenixData + "hp.obo";
        LOGGER.info("Loading OBO file {} ...", pathHpoObo);
        try {
            hpo = new HpoOboParser(new File(pathHpoObo)).parse();
        } catch (IOException e) {
            throw new PhenixException("Problem loading hp.obo file", e);
        }
        // Extract phenotypic abnormality sub ontology.
        abnormalPhenoSubOntology = hpo.getPhenotypicAbnormalitySubOntology();
        LOGGER.info("Done loading OBO file.");

        // Load file with gene-to-phenotype links.
        final String pathGeneToTermLink =
                this.pathPhenixData + "ALL_SOURCES_ALL_FREQUENCIES_genes_to_phenotype.txt";
        LOGGER.info("Loading gene-to-phenotype file {} ...", pathGeneToTermLink);
        final File fileGeneToTermLink = new File(pathGeneToTermLink);
        termToGeneAnnotations = new ArrayList<>();
        try (HpoGeneAnnotationParser parser = new HpoGeneAnnotationParser(fileGeneToTermLink)) {
            while (parser.hasNext()) {
                termToGeneAnnotations.add(parser.next());
            }
        } catch (IOException | TermAnnotationParserException e) {
            throw new PhenixException("Problem reading gene-to-phenotype file", e);
        }
        // Convert term-to-gene annotations to mapping with gene-to-terms information.
        entrezGeneIdToTermIds = TermAnnotations
                .constructTermLabelToAnnotationsMap(abnormalPhenoSubOntology, termToGeneAnnotations)
                .entrySet().stream().map(e -> {
                    final int intKey = Integer.parseInt(e.getKey().substring("ENTREZ:".length()));
                    return new AbstractMap.SimpleEntry<Integer, Collection<TermId>>(intKey,
                            e.getValue());
                }).collect(Collectors.toMap(e -> e.getKey(), e -> e.getValue()));
        LOGGER.info("Done loading gene-to-phenotype file");
    }

    /**
     * Build the Resnik similarity object and perform information content precomputation.
     */
    private void buildResnikSimilarity() {
        LOGGER.info("Performing IC precomputation...");
        final InformationContentComputation<? extends Term,
                ? extends TermRelation> icPrecomputation =
                        new InformationContentComputation<>(abnormalPhenoSubOntology);
        final Map<TermId, Collection<String>> termLabels =
                TermAnnotations.constructTermAnnotationToLabelsMap(abnormalPhenoSubOntology,
                        termToGeneAnnotations);
        final Map<TermId, Double> termToIc = icPrecomputation.computeInformationContent(termLabels);
        LOGGER.info("Done precomputing IC.");

        LOGGER.info("Performing Resnik precomputation...");
        final int numThreads = 1; // TODO: get from options?
        resnikSimilarity = new ResnikSimilarity<>(new PrecomputingPairwiseResnikSimilarity<>(
                abnormalPhenoSubOntology, termToIc, numThreads), symmetric);
        LOGGER.info("Done with Resnik precomputation.");
    }

    /**
     * @return flag to output results of filtering against Phenix.
     */
    @Override
    public PriorityType getPriorityType() {
        return PRIORITY_TYPE;
    }

    @Override
    public Stream<PhenixPriorityResult> prioritise(List<String> hpoIds, List<Gene> genes) {
        if (hpoIds.isEmpty()) {
            throw new PhenixException(
                    "Please supply some HPO terms. PhenIX is unable to prioritise genes without these.");
        }

        // Translate HPO term IDs to TermId objects, remove duplicates.
        LOGGER.info("Creating HPO query terms from {}", hpoIds);
        final List<TermId> hpoQueryTerms = makeHpoQueryTerms(hpoIds);
        LOGGER.info("Created HPO query terms {}", hpoQueryTerms);

        // Compute Resnik similarity scores.
        final Map<Gene, PhenixScore> geneScores =
                genes.stream().collect(toMap(Function.identity(), scoreGene(hpoQueryTerms)));

        // Compute normalization factor.
        final double maxSemSimScore = geneScores.values().stream()
                .mapToDouble(PhenixScore::getSemanticSimilarityScore).max().orElse(DEFAULT_SCORE);
        final double normalisationFactor = calculateNormalisationFactor(maxSemSimScore);
        LOGGER.info("Similarity computations in HPO done for {} genes", genes.size());

        // TODO(holtgrewe): Perform multiple testing correction correction to p values.
        return geneScores.entrySet().stream().map(entry -> {
            final Gene gene = entry.getKey();
            final PhenixScore phenixScore = entry.getValue();
            final double score = phenixScore.getSemanticSimilarityScore() * normalisationFactor;
            return new PhenixPriorityResult(gene.getEntrezGeneID(), gene.getGeneSymbol(), score,
                    phenixScore.getSemanticSimilarityScore(), phenixScore.getNegativeLogP());
        });
    }

    /**
     * Convert list of {@link String} HPO term IDs to a {@link List} of {@link TermId}s.
     *
     * @param hpoIds String HPO term IDS to convert.
     * @return {@link TermId}s after translation.
     */
    private List<TermId> makeHpoQueryTerms(List<String> hpoIds) {
        // Note that algorithm is O(n^2) but in practice small arrays are *fast*.
        final List<TermId> result = new ArrayList<>();
        for (String hpoId : hpoIds) {
            final TermId termId = hpo.getPrimaryTermId(ImmutableTermId.constructWithPrefix(hpoId));
            if (termId != null && !result.contains(termId)) {
                result.add(termId);
            }
        }
        return result;
    }

    /**
     * Build {@link Function} for computing gene similarities to the list of query {@link TermId}s.
     *
     * @param queryTermIds {@link TermId}s from the query.
     * @return {@link PhenixScore} with score information of this object.
     */
    private Function<Gene, PhenixScore> scoreGene(List<TermId> queryTermIds) {
        return gene -> {
            // Get the TermIds that are associated with the given Entrez Gene; short-circuit
            // if the gene is not annotated with any HPO terms.
            final Collection<TermId> geneTermIds =
                    entrezGeneIdToTermIds.get(gene.getEntrezGeneID());
            if (geneTermIds == null) {
                return new PhenixScore(DEFAULT_SCORE, DEFAULT_NEGATIVE_LOG_P);
            }

            // Compute the semantic similarity score using Resnik similarity.
            final double score = resnikSimilarity.computeScore(queryTermIds, geneTermIds);
            if (Double.isNaN(score)) {
                LOGGER.error("Score was NaN for geneID {} (query terms: {})",
                        gene.getEntrezGeneID(), queryTermIds);
            }

            // Interpolate the empirical p value.
            double negLogP = DEFAULT_NEGATIVE_LOG_P;
            try {
                final ObjectScoreDistribution scoreDist = scoreDistReader
                        .readForTermCountAndObject(queryTermIds.size(), gene.getEntrezGeneID());
                negLogP = -1 * Math.log(scoreDist.estimatePValue(score));
            } catch (OntoLibException e) {
                LOGGER.info("No score distribution for term count {} and gene ID {}.",
                        new Object[] {queryTermIds.size(), gene.getEntrezGeneID()});
                e.printStackTrace();
            }

            return new PhenixScore(score, negLogP);
        };
    }

    /**
     * The gene relevance scores are to be normalized to lie between zero and one.
     *
     * <p>
     * This function, which relies upon the variable {@link #maxSemSim} being set in
     * {@link #scoreGene}, divides each score by {@link #maxSemSim}, which has the effect of putting
     * the phenomizer scores in the range [0..1]. Note that this is not the same as rank
     * normalization!
     * </p>
     */
    private double calculateNormalisationFactor(double maxSemSimScore) {
        if (maxSemSimScore < 1) {
            return 1.0;
        } else {
            return 1.0 / maxSemSimScore;
        }
    }

    /**
     * Exception class thrown on problems in Phenix priotizer.
     */
    private static class PhenixException extends RuntimeException {

        /** UID for serialization. */
        private static final long serialVersionUID = 1L;

        /**
         * Construct with message.
         *
         * @param message Message text.
         */
        private PhenixException(String message) {
            super(message);
        }

        /**
         * Construct with message and causing {@link Throwable}.
         *
         * @param message Message text.
         * @param cause Causing {@link Throwable}.
         */
        private PhenixException(String message, Throwable cause) {
            super(message, cause);
        }

    }

    // TODO move this to the messages
    // /**
    // * @return an ul list with summary of phenomizer prioritization.
    // */
    // @Override
    // public String getHTMLCode() {
    // String s = String.format("Phenomizer: %d genes were evaluated; no phenotype
    // data available for %d of them",
    // this.analysedGenes, this.offTargetGenes);
    // String t = null;
    // if (symmetric) {
    // t = String.format("Symmetric Phenomizer query with %d terms was performed",
    // this.numberQueryTerms);
    // } else {
    // t = String.format("Asymmetric Phenomizer query with %d terms was performed",
    // this.numberQueryTerms);
    // }
    // String u = String.format("Maximum semantic similarity score: %.2f, maximum
    // negative log. of p-value: %.2f", this.maxSemSim, this.maxNegLogP);
    // return String.format("<ul><li>%s</li><li>%s</li><li>%s</li></ul>\n", s, t,
    // u);
    //
    // }

    /**
     * Helper record class for storing pair of score and neg. log. p-value.
     */
    private static final class PhenixScore {

        /** Semantic similarity score. */
        private final double semanticSimilarityScore;

        /** Negative logarithm of the p-value. */
        private final double negativeLogP;

        /**
         * Construct with the given values.
         *
         * @param semanticSimilarityScore Semantic similarity score.
         * @param negativeLogP Negative logarithm of p-value.
         */
        PhenixScore(double semanticSimilarityScore, double negativeLogP) {
            this.semanticSimilarityScore = semanticSimilarityScore;
            this.negativeLogP = negativeLogP;
        }

        /**
         * @return Semantic similarity score.
         */
        double getSemanticSimilarityScore() {
            return semanticSimilarityScore;
        }

        /**
         * @return Negative logarithm of the p-value.
         */
        double getNegativeLogP() {
            return negativeLogP;
        }

    }

    /*
     * Essentially, this is only present for testing...
     *
     * (non-Javadoc)
     * @see java.lang.Object#equals(java.lang.Object)
     */
    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        PhenixPriority that = (PhenixPriority) o;
        return symmetric == that.symmetric;
    }

}