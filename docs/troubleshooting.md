# Common issues

## A taxon I know to be detectable with my primer system is listed as not amplified?!

In-silico prediction of amplification can be affected by various factors, including the quality of the reference database. 
Generally, BarBeQuE provides indications, but should not be taken as evidence for amplification or lack thereof. 

## No results are returned when providing a specific taxon level

When using the `--taxon` argument, please make sure that the taxon corresponds to one of the supported taxonomic levels.

## My primer set is not predicted to amplify as expected

Many reasons could contribute to failed in-silico amplification. Please make sure that:
- You use a suitable reference database
- The min and max arguments of your primer set are reasonable

In addition, if you are using gene-specific databases (such as Midori), make sure that your primers bind within the target gene. If one or both primers 
bind up-/downstream of the gene, that binding motif will not be included in the gene-level database and therefore no amplification can be predicted. When 
in doubt, try using the `refseq_mito` database to produce some predicted amplicon sequences and blast these against a suitable database to see where they would map. 

## The Process HELPER_TAXONOMIC_COVERAGE times out

This process checks the taxa predicted to be amplified against all species belonging to the group you specified with `--taxon`. 
The aim is to identify which members are detected by the respective primer set, which are not and which are missing from the database. 
If this taxonomic group has too many members (= species), traversing the underlying phylogenetic tree will take a very long time - potentially 
longer than the run time that is given to this process by default. Nextflow will repeat the job should it exceed this initial limit by doubling 
the run time. Failing that, the runtime will be quadrupled. If the job still fails, the pipeline fails permanently. The only solution here 
would be to repeat the analysis for multiple, more restricted taxonomic groups. 