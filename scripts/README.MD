#plinkMerge 
plinkMerge is a helper python script to manage merging several plink bed files on intersection of snps across all the bed files that need to be merged
plinkMerge requires only one dependency
* plink 1.90 - has to be in the path- the script will check and exit with an error if not found.
plinkMerge requires the following arguments 
1. -bedList ( list of plink files no file extension to be merged - one per line)
2. -out a string output of the common snp list to be written out
3. -bedOut (a string output of the merged output)
```python ```