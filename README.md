# SpatialAutocorrelation
Software for assessing whether microbial taxa, genes, and traits have structured spatial distributions.

This software checks whether the taxa, genes, and traits in a set of geographically distributed metagenomic samples are distributed randomly, or if they are aggregated or segregated (i.e., are autocorrelated). Starting with a taxon, gene, or trait table in HDF5 BIOM format (http://biom-format.org) and metadata on the locations where samples were collected, the software calculates the Moran's I statistic for each taxon, gene, or trait, and the statistical significance of that Moran's I statistic using a randomization test. The randomization test uses a Markov chain Monte Carlo algorithm and is highly optimized to allow large data tables to be processed. The shell script 'SpatialAutocorrelationTest.sh' runs an example data set and check for errors with the installation. Addition documentation on the software is available by running

<pre><code>
java -cp Autocorrelation.jar edu.ucsf.SpatialAutocorrelation.SpatialAutocorrelationLauncher --help | less
</code></pre>

from the bin directory.

## Dependencies
This software requires Java Runtime Environment (JRE) version 1.6 or greater. On many Apple systems, even if JRE 1.6 or greater is installed, the default version for running applications may be 1.5. The Java version can be checked by typing 'java -version' into a terminal. If an updated version is installed but is not being used, a few updates will need to be made; namely you might try typing the following commands in a terminal:

<pre><code>
sudo rm /usr/bin/java
sudo ln -s /Library/Internet\ Plug-Ins/JavaAppletPlugin.plugin/Contents/Home/bin/java /usr/bin
</code></pre>

Additional information on correcting the Java version can be found here: http://stackoverflow.com/questions/12757558/installed-java-7-on-mac-os-x-but-terminal-is-still-using-version-6.
