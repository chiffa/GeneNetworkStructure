In main: 
- Extract the parsing of the weights into its own function, then module
- Replace Database weight with TRRUST
    - Problem - TRRUST does not have data for yeast
    - Problem #2 - all the other databases we found with Joel are for humans and mice
    - we could technically purify the direct relationships by performing correlation mapping
    and then structure clearing for yeast TF network
- We can as well replace the structural topology with the TFs from Yeast or pure P2P

While attempting to perform the structural data matching with the effect, run into the issue
of translating identifiers - from ENSEMBL stable IDs to Gene names - Biomart only
has a small fraction of mappings. Have to spin up a BioFlow instance.