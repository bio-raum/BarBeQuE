# Pipeline structure

The pipeline accepts one or several sets of metabarcoding primers and predicts the amplification status of species using one or several reference database(s). 

First, the primers are used to perform in-silico PCR on a given database. The results are then dereplicated and filtered (by size range) and used to calculate some useful metrics. In addition, if users specify a target taxon to profile, BarBeQue will check - for all species known to belong to that group via the NCBI taxonomy database - if that species is amplified, not amplified or missing from the database. For all hits, the users is furthermore provided with a list of any additional taxa with an identical barcode sequence (ie. redundant) - information that can be used to decide if a primer set is suitable to resolve a particular taxonomic group or not. 

The main use is in metabarcoding to test newly designed primers or check the suitability of a specific primer set for a given question (i.e. "can this be amplified with this particular primer?"). 
