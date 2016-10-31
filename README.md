[![Build Status](https://travis-ci.org/naturalis/supersmart.svg?branch=master)](https://travis-ci.org/naturalis/supersmart)

![SUPERSMART](http://www.supersmart-project.org/images/logo200x480.png "SUPERSMART")

INFO
----
This is the source code repository of the Self-Updating Platform for the Estimation 
of Rates of Speciation, Migration, And Relationships of Taxa (SUPERSMART), a 
pipeline for estimating time-calibrated species trees from public data. Constantly 
updated information about the pipeline can be retrieved from the project website at 
www.supersmart-project.org. An article about SUPERSMART has been published in the journal [Systematic Biology](http://sysbio.oxfordjournals.org/content/early/2016/10/13/sysbio.syw066.abstract).

INSTALLATION
------------
Installation depends on vagrant and Virtualbox, see INSTALL.md, the installation guide on the [wiki page](https://github.com/naturalis/supersmart/wiki/Installation-instructions) and the [project website](http://www.supersmart-project.org).

COMMANDS
--------
The functionality of the pipeline is exposed via the [smrt](script/smrt) command.
The 'smrt' command offers [twelve subcommands](https://github.com/naturalis/supersmart/wiki/Analysis-run-through), 
representing each step of the pipeline. 

CONFIGURATION
-------------
Basic configuration is done via options passed to the 'smrt' command. 
Advanced configuration of the pipeline API is possible in the 
[conf/supersmart.ini](conf/supersmart.ini) file. The [smrt-config](script/smrt-config) 
command provides an interface for changing the configuration file.

EXAMPLES
--------
The examples/* directories contain instructional examples to demonstrate the input
file formats.

BUILD STATUS
------------
Currently, the build status at Travis is:

[![Build Status](https://travis-ci.org/naturalis/supersmart.svg?branch=master)](https://travis-ci.org/naturalis/supersmart)

