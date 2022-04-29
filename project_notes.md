# Notes for ind2pop Niches Project

Code for producing empirical examples associated with draft ms using mixture distributions to characterize the relationship between individual- and population scale niches. (Lu et al *in prep*)

**Contact:**  
  Scott Yanco, PhD  
  scott.yanco@yale.edu

## TODO:  
  * Bi-variate Plots
  * Climate Vulnerability:   
    * Take mods to elephants and implement for gadwall - this will complete "past annotation"
    * Need to do future annotation - think about how to do this - essentially need a single offset value
  * Gather data permissions
  * Get nicer colors...
    


## Activity log

|Date|Activity|
|:-|:------------|
|2021-11-02|Built directory/repo, importedbreezy template, imported niche functions script|
|2021-12-01|Simplifying approach...|
|2021-12-09|Modifying workflow, working with cranes and elephants, annos on HPC|
|2022-01-11|New year!  Picking this back up, swapped storks for Gadwall and making figs|
|2022-01-12|Worked up trees, uploaded figs to draft ms|
|2022-01-13|Wrote out conda envs, cleaning project docs|

## Notes
*  Should be able to do this more simply than breezy-style or with mosey db.  Thus trying to build entire workflow as a single r script (calling a couple external function scripts)

*No STOAT pre-annotations for trees....  Request batch?  Also, though, is FIA in mol.org? Might skip the plants

*Might try LST - WJ suggested that EVI may be particularity amenable to this...

*April 13 2022 notes*  
Swicthed annotation to pull max temp from CMIP5 dataset on GEE.  This allows forwards-backward consistency in modeling.  See TODO for remaining stpes to implementation.