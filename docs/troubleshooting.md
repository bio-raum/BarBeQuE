# Common issues

## A taxon I know to be detectable with my primer system is listed as not amplified?!

In-silico prediction of amplification can be affected by various factors, including the quality of the reference database. Generally, BarBeQuE provides indications, but should not be taken as evidence for amplification or lack thereof. 

## No results are returned when providing a specific taxon level

When using the `--taxon` argument, please make sure that the taxon corresponds to one of the supported taxonomic levels.

## My primer set is not predicted to amplify as expected

Many reasons could contribute to failed in-silico amplification. Please make sure that:
- You use a suitable reference database
- The min and max arguments of your primer set are reasonable

In addition, if you are using gene-specific databases (such as Midori), make sure that your primers bind within the target gene. If one or both primers bind up-/downstream of the gene, that binding motif will not be included in the gene-level database and therefore no amplification can be predicted. When in doubt, try using the `refseq_mito` database to produce some predicted amplicon sequences and blast these against a suitable database to see where they would map. 

