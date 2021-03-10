* Version 4.5.1:
    * fix small issues (broken links, etc.) in 4.5.0
* Version 4.5.0:
    * TIPP has moved to a new repo: https://github.com/TeraTrees/TIPP
* Version 4.4.0:
    * TIPP has moved to TIPP2
* Version 4.3.21:
    * one more small bug fix
* Version 4.3.20:
    * Added `--ignore-overlap` option to ignore overlap. 
* Version 4.3.19:
    * Enable `--symfrac` from config file
    * Enable config files that overwrite the default values only when provided
* Version 4.3.18:
    * UPP speed improvement for cases with super gappy backbone alignment (more can be done)
    * Slightly better logging
* Version 4.3.17:
    * Attempt to make tests python 3.8 friendly
    * Add UPP test
* Version 4.3.16:
	* Refactor code slightly. Reduces warnings
* Version 4.3.15:
	* some py3.8 changes and some java changes  from pgrt
* Version 4.3.14:
	* Added option `-R "Nmin Nmax"` to specify the range of full length sequences to UPP
* Version 4.3.13:
    * Added backtranslation functionality to UPP for amino acid sequences.
* Version 4.3.12:
	* Added `-rt` to remove temp files	
* Version 4.3.11:
	* Make sure PASTA user options can be passed
* Version 4.3.10:
     * fix issue #70 a bug when hmmsearch fake jobs where not piped
* Version 4.3.9:
     * fix info file path in `run_sepp.sh` script
* Version 4.3.7:
     * Clean up `run-sepp.sh` to expose more options
     * Beginning steps for building more robust builds
* Version 4.3.6:
     * raise a ValueError if names of fragments to be inserted contain whitespace characters ` ` or `\t`.
* Version 4.3.5:
     * raise a ValueError if names of fragments to be inserted collide with names of reference sequences.
* Version 4.3.4:
     * Cleaning up the default (info) log  a bit
* Version 4.3.3:
     * avoid pipes for hmmsearch to deal with memory issues
     * Increase the default `-F` to 20000
* Version 4.3.2: more fixes to UPP
     * added a  seed number
     * Fixed PASTA alignment and tree file names
     * degap input files when already aligned
     * avoid checking file names twice
* Version 4.3.1: fixes to UPP
     * In the absence of fragments, outputs are generated in a consistent fashion
     * PASTA alignment and trees are given a name prefixed by the output prefix naem
* Version 4.3.0:
     * Added (Uyen Mai) an option `-M` and `-S` midpoint to break by midpoint and to stop by diameter
* Version 4.2.3:
     * Added a hard-coded --groups 10 to the pplacer runs to help with memory
* Version 4.2.2:
     * Add max chunk size
     * Fix logging messages
     * Fix a bug for the `-D` option after python 3 transition
     * Changed the pplacer step so that each run is on one fragment chunk (not all fragments), as specified using `-F`
* Version 4.2.1:
     * Fixed a bug with `run_tipp_tools.py` path
* Version 4.2.0:
     * Major change (trying again): TIPP should not include fragments that are fully unclassified in its profile. Fixed this
* Version 4.1.0:
     * Reduce verbosity by moving some logs from info to debug
* Verstion 4.0.1:
     * Small bug fixes in UPP after transition to python3
* Version 4.0.0
     * Several small bug fixes, including, the location of the config file for TIPP and UPP
     * Major change: TIPP should not include fragments that are fully unclassified in its profile. Fixed this
* Version 3.2.2:
     * Bug fix: check for info file for SEPP
* Version 3.2.1:
     * Can now handle internal node names
* Version 3.2.0:
     * Can now have a contained installation
* Version 3.1:
     * Bug fix for TIPP for python 3
* Version 3:0:
    * Moved to python3
* Version 1.1:
    * Added -d option for directories
    * Bundled all requied software
    * Changed installation so that dependencies are automatically installed
    * Moved to latest pplacer
    * Misc. bug fixes
