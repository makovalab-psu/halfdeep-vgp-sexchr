#!/bin/bash
#SBATCH --nodes=1
#SBATCH --partition=cpu                 
#SBATCH --time=48:00:00          
#SBATCH --mem=100G              
#SBATCH --ntasks=90              
START_TIME=$(date +%s)

# 0. Get parameteres
# 1. download fasta
# 2. download reads
# 3. read mapping
# 4. sam to bam
# 5. halfdeep structure
# 6. halfdeep calling
# 7. Backup to somewhere
# 8. Cleaning ()
# manual parameter: cores, minimap path, halfdeep path, genodsp path, memory

print_duration() {
  local END_TIME=$(date +%s)
  local DIFF=$((END_TIME - $1))
  local SEC=$((DIFF % 60))
  local MIN=$(((DIFF % 3600) / 60))
  local HOUR=$((DIFF / 3600))
  echo "Consuming: ${HOUR} hour ${MIN} min ${SEC} sec"
}



mkdir -p "1.data"

#### static parameters ####
cores=90
memory=90
minimap_exe="/u/bko/Program/minimap2/minimap2"



eval "$(mamba shell hook --shell bash)"
mamba activate halfdeep_env
export PATH=~/Program/halfdeep:${PATH}
export PATH=~/Program/genodsp:${PATH}

## Get parameters
if [ "$#" -ne 8 ]; then
    echo "Usage: $0 <data.No> <species name (2 words)> <accession> <ncbi asm file name> <VGP assembly id> <read DB (AWS/NCBI)> <read tech> <read file name list>"
    exit 1
fi

data_no="$1"
species="$2" # Bos_Taurus
accession="$3" # GCA_000111222.1
asm_name="$4" # 
asm_id="$5" # mBosTau1
db="$6" # AWS | NCBI
tech="$7" # HiFi, Hifiasm Hi-C phasing ~ CLR I, TrioCanu
read_list="$8" # AWSread1,AWSread2.. | NCBI SRA Accession

echo "Data No: $1"
echo "Species: $2"
echo "Accession: $3"
echo "Assembly Name: $4"
echo "VGP ID: $5"
echo "DB: $6"
echo "Tech: $7"
echo "Read List: $8"
echo -e "\n"

sleep 10
mkdir -p "1.data/$data_no.$species"
mkdir -p "fail_alert"

## get assembly fasta from NCBI
haplotype_url="${asm_name%_genomic.fna}"
asm_url="https://ftp.ncbi.nlm.nih.gov/genomes/all/${accession:0:3}/${accession:4:3}/${accession:7:3}/${accession:10:3}/${haplotype_url}/${asm_name}.gz"

echo "wget $asm_url"

mkdir -p "1.data/$data_no.$species/ncbi_assembly"

max_attempts=20
attempt=0
while [[ ! -f "1.data/$data_no.$species/ncbi_assembly/${species}.fasta.gz" && $attempt -lt $max_attempts ]]; do
    echo "[$species] NCBI report download attempt ($((attempt+1))/$max_attempts)..."
    echo "assembly download attempt of ${accession} for ${asm_name}.gz"
    ncbi_cmd="wget -O 1.data/$data_no.$species/ncbi_assembly/${species}.fasta.gz $asm_url"
    echo "$ncbi_cmd"
    eval "$ncbi_cmd" 2>/dev/null 
    sleep 5
    ((attempt++))
done

touch "1.data/$data_no.$species/ncbi_assembly/${species}.fasta.gz"
echo "${species}.fasta.gz timestamp updated to current time."
echo -e "\n"

read_path="1.data/$data_no.$species/reads/"
mkdir -p "$read_path"

print_duration $START_TIME


## Get reads
## AWS mode (Master code have get read file name by deciding HiFi fastq -> HiFi bam -> CLR bam) -> add later

if [ "$db" == "AWS" ]; then
    if [[ "$tech" =~ [Hh]i[Ff]i ]]; then
        # tech contains 'hifi', use pacbio_hifi
        folder="pacbio_hifi"
    else
        # tech does not contain 'hifi', use pacbio
        folder="pacbio"
    fi
    
    base_s3_url="s3://genomeark/species/$species/$asm_id/genomic_data/$folder/"
    IFS=',' read -ra FILE_ARRAY <<< "$read_list"
    
    # individual read or bam download by individual aws request
    for file in "${FILE_ARRAY[@]}"; do
	if [[ "$file" =~ \.pbi$ ]]; then
            continue
        fi
        echo "Downloading $file..."
	attempts=0

	while [ $attempts -lt 10 ]; do
            aws_cmd="aws s3 cp ${base_s3_url}${file} $read_path --no-sign-request"
            echo "$aws_cmd"
	    eval "$aws_cmd" > /dev/null 2>&1


	    if [ $? -eq 0 ]; then
                echo "$file downloaded successfully."
		touch "${read_path}$file"
                echo "$file timestamp updated to current time."
                
		# .pbi download
		if [[ "$file" == *.bam ]]; then
                    pbi_file="${file}.pbi"
                    echo "Downloading $pbi_file..."
                    aws_cmd="aws s3 cp ${base_s3_url}$pbi_file $read_path --no-sign-request"
                    echo "$aws_cmd"
                    eval "$aws_cmd" > /dev/null 2>&1
                    
		    # only twice trial for pbi
                    if [ $? -eq 0 ]; then
                        echo "$pbi_file downloaded successfully."
			touch "${read_path}$pbi_file"
			echo "$pbi_file timestamp updated to current time."

                    else
			echo "$aws_cmd"
			eval "$aws_cmd" > /dev/null 2>&1
			if [ $? -eq 0 ]; then
		            echo "$pbi_file downloaded successfully."
			    touch "${read_path}$pbi_file"
                            echo "$pbi_file timestamp updated to current time."
			else
                            echo "[[WARNING]] Failed to download $pbi_file by trying twice."
			    
			    if [ "$folder" == "pacbio" ]; then
				echo "[INFO] Attempting to generate .pbi index locally using pbindex..."
				pbindex_cmd="pbindex -j $cores ${read_path}${file}"
				echo "$pbindex_cmd"
                                eval "$pbindex_cmd"
				if [ $? -eq 0 ]; then
				    echo "[INFO] .pbi index successfully generated with pbindex."
				else
				    echo "[[ERROR]] pbindex failed for $file."
				fi
		            fi
			fi
                    fi
                fi
		break
	    else
		attempts=$((attempts + 1))
		echo "Attempt $attempts failed. Retrying..."
		sleep 5
	    fi
	done
	# read down failure alert
	if [ $attempts -ge 10 ]; then
	    echo "[[FAILED]] to download $file after 10 attempts."
	    echo "[[FAILED]] to download $file after 10 attempts at $(date +%s)." >> "fail_alert/$data_no.$species.fail"
	    echo "EXIT by read not downloaded"
	    exit 1
	fi
    done
elif [ "$db" == "SRA" ]; then # bug fix 20260119 AWS -> SRA (i may changed else to elif for running WGET using this pipeline, then just copy above (AWS) for this (must be SRA), but forgot to put SRA
    if [[ "$tech" =~ [Hh]i[Ff]i ]]; then
        # tech contains 'hifi', use pacbio_hifi
        folder="pacbio_hifi"
    else
        # tech does not contain 'hifi', use pacbio
        folder="pacbio"
    fi
    # read_path="1.data/$data_no.$species/reads/"
    echo "[NCBI mode] Starting prefetch + fasterq-dump"
    IFS=',' read -ra SRR_ARRAY <<< "$read_list"

    for SRR in "${SRR_ARRAY[@]}"; do
        echo "[SRR] Processing $SRR..."

        # 1. prefetch
        sra_file="${read_path}${SRR}.sra"
        cmd="prefetch --max-size 500G $SRR --output-file $sra_file"
        echo "$cmd"
        eval "$cmd"
        if [ ! -f "$sra_file" ]; then
            echo "[[ERROR]] $sra_file not found after prefetch."
            echo "$SRR failed at prefetch." >> "fail_alert/$data_no.$species.fail"
            exit 1
        fi

	cmd="vdb-validate $sra_file"
        echo "$cmd"
        eval "$cmd"
	if [ $? -ne 0 ]; then
            echo "[[ERROR]] $sra_file is invalid or corrupted."
            echo "$SRR failed at vdb-validate." >> "fail_alert/$data_no.$species.fail"
            exit 1
        fi

        # 2. fasterq-dump
        tmp_dir="${read_path}/tmp_${SRR}"
        cmd="fasterq-dump $sra_file -e $cores --mem ${memory}000000000 -t $tmp_dir -O $read_path"
        echo "$cmd"
        eval "$cmd"
        fastq_file="${read_path}/${SRR}.fastq"

        if [ $? -eq 0 ] && [ -s "$fastq_file" ]; then
            echo "[SRR] fasterq-dump succeeded for $SRR."
            cmd="rm -rf $tmp_dir"
            echo "$cmd"
            eval "$cmd"

            if [[ "$tech" =~ [Hh]i[Ff]i ]]; then
                filtered_fastq="${read_path}/${SRR}.filtered.fastq"
                fasta_out="${read_path}/${SRR}.fasta"
                cmd="seqkit seq -Q 20 $fastq_file > $filtered_fastq" # extracthifi also >= Q20. Don't know mean or min
                echo "$cmd"
                eval "$cmd" || { echo "[[Error]] seqkit filtering failed"; exit 1; }

                cmd="seqkit fq2fa $filtered_fastq -o $fasta_out"
                echo "$cmd"
                eval "$cmd" || { echo "[[Error]] FASTQ to FASTA conversion failed"; exit 1; }
		cmd="rm $fastq_file; rm $filtered_fastq"
		echo "$cmd"
                eval "$cmd"
	    else
		fasta_out="${read_path}/${SRR}.fasta"
                cmd="seqkit fq2fa $fastq_file > $fasta_out"
                echo "$cmd"
                eval "$cmd" || { echo "[[Error]] FASTQ to FASTA conversion failed for CLR"; exit 1; }
		cmd="rm $fastq_file"
		echo "$cmd"
                eval "$cmd"
	    fi

            continue
	fi
	echo "[$SRR] fasterq-dump failed or no FASTQ." >> "fail_alert/$data_no.$species.fail"
	exit 1
    done    
else
    if [[ "$tech" =~ [Hh]i[Ff]i ]]; then
        # tech contains 'hifi', use pacbio_hifi
        folder="pacbio_hifi"
    else
        # tech does not contain 'hifi', use pacbio
        folder="pacbio"
    fi
    echo "Skipping read download"
fi

print_duration $START_TIME


## Bam2Fastq (hifi.bam or CLR.bam)
if [ -z "$(ls -A $read_path)" ]; then
    echo "[[Error]] No files found in $read_path"  >> "fail_alert/$data_no.$species.fail"
    exit 1
fi

if ls "$read_path" | grep -E '\.(fastq|fastq\.gz|fq|fq\.gz|fasta|fasta\.gz)$' > /dev/null; then
    echo "Found fastx. Skipping bam2fasta."
else
    echo "No fastq/fasta found."

    if [[ "$tech" =~ [Hh]i[Ff]i ]]; then
        echo ".bam for HiFi detected. Running extracthifi step first..."

	IFS=',' read -ra FILE_ARRAY <<< "$read_list"

	for file in "${FILE_ARRAY[@]}"; do

	    if [[ "$file" =~ \.pbi$ ]]; then
                continue
            fi

	    ccs_bam="${read_path}${file}"
            hifi_bam="${read_path}${file%.bam}.hifi.bam"
           
            extracthifi_cmd="extracthifi -j $cores $ccs_bam $hifi_bam"
            echo "$extracthifi_cmd"
            eval "$extracthifi_cmd"

            if [ $? -ne 0 ]; then
                echo "[[Error]] extracthifi failed for $ccs_bam" >> "fail_alert/$data_no.$species.fail"
                exit 1
            fi

	    echo "Creating .pbi index for $hifi_bam"
	    pbindex_cmd="pbindex -j $cores $hifi_bam"
	    echo "$pbindex_cmd"
            eval "$pbindex_cmd"

            if [ $? -ne 0 ]; then
                echo "[[Error]] pbindex failed for $hifi_bam" >> "fail_alert/$data_no.$species.fail"
                exit 1
            fi
        done

	echo "Running bam2fasta on extracted HiFi BAMs..."
	bam2fax_cmd="bam2fasta -u -j $cores -o ${read_path}${species} ${read_path}*.hifi.bam"
	echo "$bam2fax_cmd"
        eval "$bam2fax_cmd"

	if [ $? -ne 0 ]; then
            echo "[[Error]] bam2fasta failed on extracted hifi.bam" >> "fail_alert/$data_no.$species.fail"
            exit 1
        fi
    else
	echo "Not HiFi. Running bam2fasta directly on BAM files..."
        bam2fax_cmd="bam2fasta -u -j $cores -o ${read_path}${species} ${read_path}*.bam" # 
        echo "$bam2fax_cmd"
        eval "$bam2fax_cmd"

	if [ $? -ne 0 ]; then
            echo "[[Error]] bam2fasta failed on bam files" >> "fail_alert/$data_no.$species.fail"
            exit 1
        fi
    fi
fi

print_duration $START_TIME


## mapping reads
asm_path="1.data/$data_no.$species/ncbi_assembly/${species}"
fastx_path="1.data/$data_no.$species/reads/"
mapping_path="1.data/$data_no.$species/mapping/"

mkdir -p "$mapping_path"

indexing_cmd="$minimap_exe -d $asm_path.mmi ${asm_path}.fasta.gz"
echo "$indexing_cmd"
eval "$indexing_cmd"

shopt -s nullglob
read_files=(1.data/"$data_no"."$species"/reads/*.{fasta,fastq,fq,fasta.gz,fastq.gz,fq.gz})

if [ "$folder" = "pacbio_hifi" ]; then
    mapping_cmd="$minimap_exe -x map-hifi -a -t $cores $asm_path.mmi ${read_files[@]} > ${mapping_path}${species}.sam"
    #mapping_cmd="$minimap_exe -x map-hifi -a -t $cores $asm_path.mmi ${fastx_path}*.{fasta,fastq,fasta.gz,fastq.gz} > ${mapping_path}${species}.sam"
elif [ "$folder" = "pacbio" ]; then
    mapping_cmd="$minimap_exe -x map-pb -a -t $cores $asm_path.mmi ${read_files[@]} > ${mapping_path}${species}.sam"
    #mapping_cmd="$minimap_exe -x map-pb -a -t $cores $asm_path.mmi ${fastx_path}*.{fasta,fastq,fasta.gz,fastq.gz} > ${mapping_path}${species}.sam"
else
    echo "Error: Unexpected read type '$folder' it would be a critical error on data or parsing"
    exit 1
fi

echo "$mapping_cmd"
eval "$mapping_cmd"

if [ $? -ne 0 ]; then
    echo "[[Error]] minimap2 failed: Non-zero exit status" >> "fail_alert/$data_no.$species.fail"
    exit 1
fi

print_duration $START_TIME


## sorting
sam_file="${mapping_path}${species}.sam"
sorted_bam_file="${mapping_path}${species}.bam"

sort_cmd="sambamba view -t $cores -l 1 -S -f bam $sam_file | sambamba sort -t $cores -l 1 -m ${memory}G -o $sorted_bam_file /dev/stdin"
echo "$sort_cmd"
eval "$sort_cmd"
if [ $? -ne 0 ]; then
    echo "[[Error]] sambamba view or sort failed: Non-zero exit status" >> "fail_alert/$data_no.$species.fail"
    exit 1
fi

if [ ! -s "$sorted_bam_file" ]; then
    echo "[[FAILED]] BAM sorting failed: Sorted BAM file not created." >> "fail_alert/$data_no.$species.fail"
    exit 1
fi

echo "[INFO] Removing intermediate SAM file: $sam_file"
rm -f "$sam_file"

print_duration $START_TIME


## depth calling
depth_path="1.data/$data_no.$species/Merged.depth.dat.gz"
depth_cmd="samtools depth -Q 1 -@ $cores $sorted_bam_file | pigz -p $cores > $depth_path"
echo "$depth_cmd"
eval "$depth_cmd"

if [ $? -ne 0 ]; then
    echo "[[Error]] samtools depth or pigz failed: Non-zero exit status" >> "fail_alert/$data_no.$species.fail"
    exit 1
fi

print_duration $START_TIME


## halfdeep structure
cd "1.data/$data_no.$species"

asm_inner_path="ncbi_assembly/${species}.fasta.gz"
ref=`echo $asm_inner_path | sed 's/.fasta$//g' | sed 's/.fa$//g' | sed 's/.fsa_nt$//g' | sed 's/.fasta.gz$//g' | sed 's/.fa.gz$//g' | sed 's/.fsa_nt.gz$//g'`
refbase=`basename $ref`


mkdir halfdeep
mkdir halfdeep/$refbase

hd_base="halfdeep/${refbase}/"
mv "Merged.depth.dat.gz" $hd_base


## halfdeep running
hd_cmd="halfdeep_ko.sh $asm_inner_path"
echo "$hd_cmd"
eval "$hd_cmd"

if [ $? -ne 0 ]; then
    echo "[[Error]] halfdeep.sh pipeline error: Non-zero exit status" >> "fail_alert/$data_no.$species.fail"
    exit 1
fi

if [ ! -s "${hd_base}halfdeep.dat" ]; then
    echo "[[FAILED]] halfdeep running failed by some reason." >> "fail_alert/$data_no.$species.fail"
    exit 1
fi

cd ../../


echo "wrapper script all DONE successfully"
print_duration $START_TIME



END_TIME=$(date +%s)
DIFFERENCE=$((END_TIME - START_TIME))
SECONDS=$((DIFFERENCE % 60))
HOURS=$((DIFFERENCE / 3600))
MINUTES=$(((DIFFERENCE % 3600) / 60))

echo "$HOURS:$MINUTES:$SECONDS"
