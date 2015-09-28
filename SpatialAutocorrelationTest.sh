#!/bin/bash

#loading local variables
sDir=$HOME/Documents/Research/Java/Distribution/SpatialAutocorrelation

#running test
echo ''
echo 'Running tests...'
echo ''
cd $sDir
java -cp bin/Autocorrelation.jar edu.ucsf.SpatialAutocorrelation.SpatialAutocorrelationLauncher --sDataPath=$sDir/test/TestData.biom --sOutputPath=$sDir/test/TestOutput.csv --iMCMCChains=10 --iMCMCIterations=10000 --rgsDistanceNeighborhoods='0-150' --bOutputData=true --iPrevalenceMinimum=0 --bNormalize=false

#checking that output is identical
iDataDifferences=`diff test/TestOutput.data.correct.csv test/TestOutput.data.csv | wc -l`

#TODO need to update code with random seed for randomizations so randomizations are replicable
cut -d\, -f1-2,4-7,10-12 test/TestOutput.data.correct.csv > test/temp.1.csv
cut -d\, -f1-2,4-7,10-12 test/TestOutput.data.csv > test/temp.2.csv
iOutputDifferences=`diff test/temp.1.csv test/temp.2.csv | wc -l`

echo ''
echo '*************************************'
if [ $iDataDifferences == '0' ]
then
	echo 'Data output test passed.'
else
	echo 'ERROR: Data output test failed.'	
fi
if [ $iOutputDifferences == '0' ]
then
	echo 'Analysis test passed.'
else
	echo 'ERROR: Analysis test failed.'	
fi
if [ $iDataDifferences == '0' ] && [ $iOutputDifferences == '0' ]
then
	echo 'All tests passed.'
else
	echo 'ERROR: One or more tests failed.'
fi
echo '*************************************'
echo ''

#cleaning up
rm test/temp.*.csv
rm test/TestOutput.csv
rm test/TestOutput.data.csv


