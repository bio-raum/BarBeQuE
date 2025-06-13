# Outputs 

## Reports

<details markdown=1>
<summary>reports</summary>

`consensus` - consensus call of raw taxon assignments

`raw` 

- `crabs` -  raw output from CRABS (insilico pcr, rereplication, filtering)

- `vsearch` - raw output from VSEARCH (clustering of in-silico OTUs)

Additional , when `--taxon`is specified:

`subset` - CRABS result for a specific taxonomic subset (tsv, png)

`tax_coverage` -  Comparison of taxonomic assignments and the taxa represented in the database (tsv, nwk)

`pipeline_info` - traces and logs as well as the input sample sheet

</details>

<details markdown=1>
<summary>TreeViewer support</summary>
BarBeQuE can produce files that are compatible with [TreeViewer](https://treeviewer.org/) - if the pipeline is run with the `--taxon` argument. 
The relevant files are located in the `tax_coverage` subfolder. 

- Open TreeViewer and load the tree file (.nwk)
- Under the "Attachment" tab, select "Add attachment" and select the data file (.tsv)
- Under Modules, select "Add further transformation" and select "Parse node states". 
- Under Modules, select "Further transformations", 
  - Select the "Parse node states" transformation you just created
  - Select the previously attached file as "Data file"
  - As separator, enter "\t" (tab)
  - Under "New attribute", select "Use first row as header"
  - Apply
  
When you zoom in the tree, you should now see colored branches indicating the state of the attached taxon
- green: will likely be amplified
- brown: will likely not be amplified
- grey: taxon missing from the database
</details>

## Pipeline run metrics

<details markdown=1>
<summary>pipeline_info</summary>

This folder contains the pipeline run metrics

- pipeline_dag.svg - the workflow graph (only available if GraphViz is installed)
- pipeline_report.html - the (graphical) summary of all completed tasks and their resource usage
- pipeline_report.txt - a short summary of this analysis run in text format
- pipeline_timeline.html - chronological report of compute tasks and their duration
- pipeline_trace.txt - Detailed trace log of all processes and their various metrics

</details>
