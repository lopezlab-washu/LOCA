#!/bin/bash

set -e

usage=$'\n\t\033[1m*** LOng read Copy-back viral genome Analysis ***\033[0m\n
Identification of cbVG reads based on blast alignment against virus reference genome.\n
Usage:\n\t\033[1m'$(basename "$0")$'\033[0m [-h] \033[1m-v VIRUS_DB -t VIRUS_trailer_DB -f FASTQ_FILE -g GENOME_LENGTH\033[0m [-d nsVG_TYPE] [-n NUMBER]
where:
\t \033[1m-v VIRUS_DB  full length virus (blast db)\033[0m
\t \033[1m-t VIRUS_trailer_DB  trailer sequence of the virus (blast db)\033[0m
\t \033[1m-f FASTQ_FILE  raw fastq file\033[0m
\t \033[1m-g GENOME_LENGTH  standard virus genome length\033[0m
\t -h  show this help text and exit\n'

DVGTYPE="CB"
N=5

options=':hH:v:t:f:g:n:d:'
while getopts $options option; do
  case "$option" in
    h) echo "$usage"; exit;;
    v) VIRUSDB=$OPTARG;;
    t) VIRUSDBtrailer=$OPTARG;;
    f) FASTQ=$OPTARG;;
    g) GL=$OPTARG;;
    :) echo "$usage"; echo -e "\033[5m\033[4mERROR\033[0m\033[5m:\033[0m missing argument for $OPTARG\n"; exit 1;;
   \?) echo "$usage"; echo -e "\033[5m\033[4mERROR\033[0m\033[5m:\033[0m illegal option -$OPTARG\n"; exit 1;;
  esac
done

# mandatory arguments
if [ ! "$VIRUSDB" ] | [ ! "$VIRUStrailer" ] | [ ! "$FASTQ" ] | [ ! "$GL" ]; then
  echo "$usage"
  echo -e "\033[5m\033[4mERROR\033[0m\033[5m:\033[0m options -H, -v, -t, f and -g must be provided.\n"
  exit 1
fi

echo
echo -e "\t* BLAST ${DVGTYPE}-VG validation analysis *"
echo

sample=$(basename ${FASTQ} .fastq.gz)

# Fastq to Fasta
seqtk seq -a ${FASTQ} > ${sample}.fasta

# RUN BLAST VS VIRUS TRAILER
blastn -db ${VIRUSDBtrailer} -query ${sample}.fasta -out ${sample}_trailer.blast -outfmt "6 qseqid qstart qend sseqid sstart send sstrand"

# Extract aligned reads
cut -f 1 ${sample}_trailer.blast | sort | uniq > ${sample}_trailer_read_ids.txt
seqtk subseq ${FASTQ} ${sample}_trailer_read_ids.txt > ${sample}_trailer.fastq && gzip ${sample}_trailer.fastq

# RUN BLAST VS VIRUS (full length)
blastn -db ${VIRUSDB} -query ${sample}.fasta -out ${sample}_Virus.blast -outfmt "6 qseqid qstart qend sseqid sstart send sstrand"

# Extract unaligned reads
cut -f 1 ${sample}_Virus.blast | sort | uniq > ${sample}_Virus_read_ids.txt
grep '^>' ${sample}.fasta | sed 's/>//' > ${sample}_read_ids.txt
grep -v -f ${sample}_Virus_read_ids.txt ${sample}_read_ids.txt > ${sample}_nonVirus_read_ids.txt
seqtk subseq ${FASTQ} ${sample}_nonVirus_read_ids.txt > ${sample}_nonVirus.fastq && gzip ${sample}_nonVirus.fastq

# CHECK BLAST OUTPUT
echo -e "Checking BLAST results..."
# 1. sort blast output
sort -k1,1 -k5,5n ${sample}_Virus.blast | tr -d '' > ${sample}.blast.sorted
# 2. count number of ranges (e.g. nsVG with 3 ranges " 3 blast_data")
cut -f 1 ${sample}.blast.sorted | uniq -c > ${sample}.blast.counts
# 3. check list of ranges values (1,2,3,4,17,etc.). If no value > 1, then report 0 DVG found.
if [[ ! $(less ${sample}.blast.counts | grep -Eoh "^ +[0-9]+ " | tr -d " " | sort -n | uniq | grep -v ^1$) ]]; then
  # print report
  echo " * Found 0 cbVG reads *"
  # delete intermediate files
  rm ${sample}.blast.*
  rm ${samplefa} ${SAMPLETXT}
  echo
  echo "Finished."
  echo
  exit 0
fi
# 4. Store list of range values (1,2,3,4,17,etc.). Reads with only 1 range are removed from list.
less ${sample}.blast.counts | grep -Eoh "^ +[0-9]+ " | tr -d " " | sort -n | uniq | grep -v ^1$ > ${sample}.blast.uniq.ranges
# 5. Extract reads with at least 2 ranges
less ${sample}.blast.counts | sed  "s/^ *//" | sort -n -k 1 | grep -v "^1 " > ${sample}.blast.junctions
# 6. Split reads list according to number of ranges (e.g. reads with 3 ranges are saved into a file named sample.blast.junctions.3)
while read i ; do grep -E "^${i} " ${sample}.blast.junctions | cut -d " " -f 2 > ${sample}.blast.junctions.${i} ; done < ${sample}.blast.uniq.ranges 
# 7. Keep only reads with 2 ranges
sort ${sample}.blast.junctions.2 > ${sample}.blast.junctions.2.sorted
# 8. Fetch strand orientation info for reads with 2 ranges
join -j1 -o "0,1.7" ${sample}.blast.sorted ${sample}.blast.junctions.2.sorted > ${sample}.blast.junctions.2.strands
# 9. Check ref strand orientation (sstrand values for range1,range2 are plus,minus or minus,plus for CB (n=1) and plus,plus or minus,minus for DEL (n=2))
uniq -c ${sample}.blast.junctions.2.strands | sed  "s/^ *//" | grep "^${n} " | cut -d " " -f 2 | uniq > ${sample}.blast.junctions.2.checked
# 10. Extract reads with 2 ranges and right strand orientation from blast output
join -j1 ${sample}.blast.sorted ${sample}.blast.junctions.2.checked | tr " " "	" > ${sample}.blast.data

### NEW TABLE OUTPUT ###
# cbVG
if  [ ${DVGTYPE} = "CB" ]; then
  i=0
  # going over each line from the blast validated data
  while read li ; do
    # range 1 (r1)
    if [[ $(( $i % 2 )) == 0 ]]; then
      # store read start and end positions for control
      r1s=$(echo ${li} | cut -d " " -f 5);
      r1e=$(echo ${li} | cut -d " " -f 6);
      # r1 start and end positions
      range1start=$(( ${r1s}<${r1e} ? ${r1s} : ${r1e} ));
      range1end=$(( ${r1s}>${r1e} ? ${r1s} : ${r1e} ));
      # check cbVG orientation
      strand1=$(echo ${li} | cut -d " " -f 7)
    # range 2 (r2)
    else
      # check read start and end positions
      r2s=$(echo ${li} | cut -d " " -f 5);
      r2e=$(echo ${li} | cut -d " " -f 6);
      # r2 start and end positions
      range2start=$(( ${r2s}<${r2e} ? ${r2s} : ${r2e} ));
      range2end=$(( ${r2s}>${r2e} ? ${r2s} : ${r2e} ));
      # strand
      strand2=$(echo ${li} | cut -d " " -f 7)
      if [ ${strand2} != ${strand1} ]; then
          BREAK=$(( ${range1start}<${range2start} ? ${range1start} : ${range2start} ))
          REJOIN=$(( ${range1start}>${range2start} ? ${range1start} : ${range2start} ))
          # get cbVG size and extra info
          readid=$(echo ${li} | cut -d " " -f 1)
          cbVGsize=$(awk "BEGIN{print ${GL}-${BREAK}+${GL}-${REJOIN}+2}")
          loopsize=$(awk "BEGIN{print ${REJOIN}-${BREAK}}")
          stemsize=$(awk "BEGIN{print (${GL}-${REJOIN}+1)*2}")
          perstem=$(awk "BEGIN{print ${stemsize}/${cbVGsize}*100}" )
          # write info to file
          echo -e "${readid}\t${cbVGsize}\t${BREAK}\t${REJOIN}\t${loopsize}\t${stemsize}\t${perstem}"
      fi
    fi
    i=$((i+1))
  done < ${sample}.blast.data > ${sample}.${DVGTYPE}.info.N${N}.tmp
fi

# fetch list of nsVG read sequences
# 1. get list of read IDs
cut -f 1 ${sample}.${DVGTYPE}.info.N${N}.tmp > ${sample}.${DVGTYPE}.readids.tmp
# 2. fetch sequences from FASTQ file
seqtk subseq -t ${FASTQ} ${sample}.${DVGTYPE}.readids.tmp | cut -f 1,3 > ${sample}.${DVGTYPE}.tab.tmp
# 3. merge with report
join <(sort ${sample}.${DVGTYPE}.info.N${N}.tmp) <(sort ${sample}.${DVGTYPE}.tab.tmp) > ${sample}.${DVGTYPE}.info.N${N}.seq.tmp

# get species occurence/break/rejoin/size/ (species is defined by SIZE and B/R) and sort by SIZE and B
sort -k2n -k3n ${sample}.${DVGTYPE}.info.N${N}.seq.tmp > ${sample}.${DVGTYPE}.info.N${N}.sorted.tmp
# initialize variables: occurence, break and rejoin
# flag variable is used for flagging equivalent B and R in order to aggregate them
# o1/O: junction occurence, b1/B/r1/R: BREAK and REJOIN positions, s1/S: DVG size
b1=0; r1=0; s1=0; species=(); flag=0
while read i; do
  S=$(echo ${i} | cut -d ' ' -f 2)
  B=$(echo ${i} | cut -d ' ' -f 3)
  R=$(echo ${i} | cut -d ' ' -f 4)
  if [[ ${s1} == ${S} ]]; then # same size
    # start new species group
    if [[ ${flag} == 0 ]]; then # prev species is new
      if (( ${B} <= $(( ${b1} + ${N} )) )); then # current species is similar to prev
        # store new B/R/S and species list
        b1=${B}; r1=${R}; s1=${S}; species+=("${i}")
        # flag updated species
        flag=1
      else # current species is different
        # print prev species info to file
        OLDIFS=$IFS; IFS=""; for j in "${species[@]}"; do
          echo -e "${j} ${b1}_${r1}"
        done; IFS=$OLDIFS
        # store current species info
        b1=${B}; r1=${R}; s1=${S}; species=("${i}")
      fi
    # update previous species group
    elif (( ${B} <= $(( ${b1} + ${N} )) )); then # current species is similar to prev
      # store new B/R/S and sum of occurences
      b1=${B}; r1=${R}; s1=${S}; species+=("${i}")
    else # current species is different
      # print prev species info to file
      OLDIFS=$IFS; IFS=""; for j in "${species[@]}"; do
        echo -e "${j} ${b1}_${r1}"
      done; IFS=$OLDIFS
      # store current species info
      b1=${B}; r1=${R}; s1=${S}; species=("${i}")
      # flag new species
      flag=0
    fi
  else # current species is different
    #print prev species info to file
    OLDIFS=$IFS; IFS=""; for j in "${species[@]}"; do
      echo -e "${j} ${b1}_${r1}"
    done; IFS=$OLDIFS
    # store current species info
    b1=${B}; r1=${R}; s1=${S}; species=("${i}")
    flag=0
  fi
done < ${sample}.${DVGTYPE}.info.N${N}.sorted.tmp >> ${sample}.${DVGTYPE}.info.N${N}.txt

# Don't forget the last species (not printed at the end of the loop)
lastrow=$(tail -n 1 ${sample}.${DVGTYPE}.info.N${N}.sorted.tmp)
B=$(echo ${lastrow} | cut -d ' ' -f 3)
R=$(echo ${lastrow} | cut -d ' ' -f 4)
# print to file
OLDIFS=$IFS; IFS=""; for j in "${species[@]}"; do
  echo -e "${j} ${B}_${R}" >> ${sample}.${DVGTYPE}.info.N${N}.txt
done; IFS=$OLDIFS
# add header
sed -i "1i READ_ID SIZE BREAK_POSITION REJOIN_POSITION LOOP_SIZE STEM_SIZE PERC_STEM SEQUENCE SPECIES" ${sample}.${DVGTYPE}.info.N${N}.txt

# Generate species mode and B_R plot
Rscript nsVG_species_plot_fix.R ${sample}.${DVGTYPE}.info.N${N}.txt ${GL} ${sample} ${DVGTYPE}

# delete intermediate files
rm ${sample}.blast.* *tmp #${sample}.${DVGTYPE}.info.N${N}.txt

# Extract unaligned reads
cut -f 1 ${sample}_Virus.blast | sort | uniq > ${sample}_Virus_read_ids.txt
grep '^>' ${sample}.fasta | sed 's/>//' > ${sample}_read_ids.txt
grep -v -f ${sample}_Virus_read_ids.txt ${sample}_read_ids.txt > ${sample}_nonVirus_read_ids.txt

# Extract cbVG reads
cut -d " " -f 1 ${sample}.${DVGTYPE}.info.N${N}.txt | sort | uniq > ${sample}_cbVG_read_ids.txt
seqtk subseq ${FASTQ} ${sample}_cbVG_read_ids.txt > ${sample}_cbVG.fastq && gzip ${sample}_cbVG.fastq

# Extract FL reads
grep -v -f ${sample}_cbVG_read_ids.txt ${sample}_Virus_read_ids.txt > ${sample}_VirusFL_read_ids.txt
seqtk subseq ${FASTQ} ${sample}_VirusFL_read_ids.txt > ${sample}_VirusFL.fastq && gzip ${sample}_VirusFL.fastq


echo "nb_reads nb_nonvirus nb_trailer nb_Virus nb_cbVG_read nb_cbVG_species" > ${sample}_report_read_counts.txt
nbreads=$(wc -l ${sample}_read_ids.txt | cut -d " " -f 1)
nbnonvirus=$(wc -l ${sample}_nonVirus_read_ids.txt | cut -d " " -f 1)
nbtrailer=$(wc -l ${sample}_trailer_read_ids.txt | cut -d " " -f 1)
nbvirus=$(wc -l ${sample}_VirusFL_read_ids.txt | cut -d " " -f 1)
nbcbvgread=$(wc -l ${sample}_cbVG_read_ids.txt | cut -d " " -f 1)
nbcbvgspecies=$(cut -f 10 ${sample}.${DVGTYPE}.info.N${N}_mode.txt | grep -v mode | sort | uniq | wc -l)

echo "${nbreads} ${nbnonvirus} ${nbtrailer} ${nbvirus} ${nbcbvgread} ${nbcbvgspecies}" >> ${sample}_report_read_counts.txt

echo
echo "Finished."
echo
