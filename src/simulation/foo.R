#!/usr/bin/env Rscript

'
Usage: my_prog2.R [-a -r -m <msg>]

Options:
 -a        Add
 -r        Remote
 -m <msg>  Message
' -> doc

library(docopt, lib.loc = "/vega/stats/users/smm2253/rpackages")
opts <- docopt(doc)
str(opts)

