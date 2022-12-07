# Copy the code
# Exchange the lines in ONCOCNV.sh below "Set path arguments" until "Run commands"
# in version 6.9 these are the lines 41-58
# Path are then set automatically during the TMLflow

TOOLDIR=$(realpath $5)            #set a path to the directory with CNV detection scripts. In your case, it can be, for example, ONCOCNV.v5.5/
echo $TOOLDIR
DATADIR=$(realpath $6)             #set a path to the directory with .BAM files
echo $DATADIR
OUTDIR=$(realpath $4)        #set a path to the *existing* directory to write the output
echo $OUTDIR
GENOME=$(realpath $3)            #set a path to the .fasta file with the reference genome
echo $GENOME

targetBed="$7"      #.bed file with start and end position of amplicons (amplicons overlaping for more than 75% will be automatically merged)
echo $targetBed

#set control and test samples (located in $DATADIR):

controls="$2"
echo $controls
tests="$1"
echo $tests
