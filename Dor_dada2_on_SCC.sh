####################################
#Preparations
####################################
#connect to the SCC with port forwarding
ssh efadeev@scc1.bu.edu

#activate environment
conda activate dada2-v1.14

export PATH=/usr2/collab/efadeev/.local/bin:$PATH
module load R/3.6.2

# set variables, program and script locations
FILELOCATION="/projectnb/npsegre/efadeev_playground/For_DS"
NSAMPLE="96"

cd $FILELOCATION

###step 0: shorten file names
# save original file names
ls -1v ./Original/*R1_001.fastq.gz > ./originalR1
ls -1v ./Original/*R2_001.fastq.gz > ./originalR2

mkdir Renamed

# copy original files to new file names
new=1
for i in $(ls -1v ./Original/*R1_001.fastq.gz)
do
cp ${i} ./Renamed/${new}"_R1.fastq.gz"
((new++))
done

new=1
for i in $(ls -1v ./Original/*R2_001.fastq.gz)
do
cp ${i} ./Renamed/${new}"_R2.fastq.gz"
((new++))
done

# check that the renaming was done correctly
ls -1v ./Renamed/[0-9]*_R1.fastq.gz > ./renamedR1
ls -1v ./Renamed/[0-9]*_R2.fastq.gz  > ./renamedR2

paste ./originalR1 ./renamedR1 > ./fileID_R1
paste ./originalR2 ./renamedR2 > ./fileID_R2

#the following commands schould not give any output
while read line ; do
diff $(echo "$line")
done < ./fileID_R1

while read line ; do
diff $(echo "$line")
done < ./fileID_R2

###step 1: primer clipping 
# bacterial primer V4-V5
FBC=^GTGYCAGCMGCCGCGGTAA # forward primer not anchored due to the mixed versions
RBC=^CCGYCAATTYMTTTRAGTTT # reverse primer
OFWD=18 # length of forward primer (17) - 1
OREV=19 # length of reverse primer (21) - 1
ERROR=0.16

mkdir Clipped
mkdir Logfiles

#run primers clipping
qsub -wd ${FILELOCATION} -t 1-${NSAMPLE} -o ${FILELOCATION}/Logfiles \
-e ${FILELOCATION}/Logfiles -v FBC=${FBC},RBC=${RBC},OFWD=${OFWD},OREV=${OREV},ERROR=${ERROR} \
../primers_clipping.sh 