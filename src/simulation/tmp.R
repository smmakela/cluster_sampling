#!/usr/bin/env Rscript

'usage: my_prog.R --really_long_option_name_one --another_really_long_option_name --last_long_name

options: 
 
 -r --really_long_option_name_one      a very long descritipion that goes on to more than one line [default: this/file/path]
 -a --another_really_long_option_name       Add
 -l --last_long_name        Long
' -> doc

#library(docopt, lib.loc = "/vega/stats/users/smm2253/rpackages")
#opts <- docopt(doc)
#str(opts)



# load the docopt library
library(docopt, lib.loc = "/vega/stats/users/smm2253/rpackages")
opts <- docopt(doc)
str(opts)

