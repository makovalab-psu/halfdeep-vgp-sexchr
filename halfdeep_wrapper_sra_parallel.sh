#!/bin/bash
#SBATCH --nodes=1
#SBATCH --partition=cpu                 
#SBATCH --time=48:00:00         
#SBATCH --mem=150G              
#SBATCH --ntasks=100              
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
cores=100
memory=140
splitread=10
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
#while [[ ! -f "1.data/$data_no.$species/ncbi_assembly/${species}.fasta.gz" && $attempt -lt $max_attempts ]]; do
#    echo "[$species] NCBI report download attempt ($((attempt+1))/$max_attempts)..."
#    echo "assembly download attempt of ${accession} for ${asm_name}.gz"
#    ncbi_cmd="wget -O 1.data/$data_no.$species/ncbi_assembly/${species}.fasta.gz $asm_url"
#    echo "$ncbi_cmd"
#    eval "$ncbi_cmd" 2>/dev/null 
#    sleep 5
#    ((attempt++))
#done

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
else
    if [[ "$tech" =~ [Hh]i[Ff]i ]]; then
        # tech contains 'hifi', use pacbio_hifi
        folder="pacbio_hifi"
    else
        # tech does not contain 'hifi', use pacbio
        folder="pacbio"
    fi
    echo "[NCBI mode] Starting prefetch + fasterq-dump"
    IFS=',' read -ra SRR_ARRAY <<< "$read_list"

    #for SRR in "${SRR_ARRAY[@]}"; do
    #    echo "[SRR] Processing $SRR..."

    #    # 1. prefetch
    #    sra_file="${read_path}${SRR}.sra"
    #    cmd="prefetch --max-size 500G $SRR --output-file $sra_file"
    #    echo "$cmd"
    #    eval "$cmd"
    #    if [ ! -f "$sra_file" ]; then
    #        echo "[[ERROR]] $sra_file not found after prefetch."
    #        echo "$SRR failed at prefetch." >> "fail_alert/$data_no.$species.fail"
    #        exit 1
    #    fi

    #    cmd="vdb-validate $sra_file"
    #    echo "$cmd"
    #    eval "$cmd"
    #    if [ $? -ne 0 ]; then
    #        echo "[[ERROR]] $sra_file is invalid or corrupted."
    #        echo "$SRR failed at vdb-validate." >> "fail_alert/$data_no.$species.fail"
    #        exit 1
    #    fi

    #    # 2. fasterq-dump
    #    tmp_dir="${read_path}/tmp_${SRR}"
    #    cmd="fasterq-dump $sra_file -e $cores --mem ${memory}000000000 -t $tmp_dir -O $read_path"
    #    echo "$cmd"
    #    eval "$cmd"
    #    fastq_file="${read_path}/${SRR}.fastq"

    #    if [ $? -eq 0 ] && [ -s "$fastq_file" ]; then
    #        echo "[SRR] fasterq-dump succeeded for $SRR."
    #        cmd="rm -rf $tmp_dir"
    #        echo "$cmd"
    #        eval "$cmd"

    #        if [[ "$tech" =~ [Hh]i[Ff]i ]]; then
    #            filtered_fastq="${read_path}/${SRR}.filtered.fastq"
    #            fasta_out="${read_path}/${SRR}.fasta"
    #            cmd="seqkit seq -Q 20 $fastq_file > $filtered_fastq" # extracthifi also >= Q20. Don't know mean or min
    #            echo "$cmd"
    #            eval "$cmd" || { echo "[[Error]] seqkit filtering failed"; exit 1; }

    #            cmd="seqkit fq2fa $filtered_fastq -o $fasta_out"
    #            echo "$cmd"
    #            eval "$cmd" || { echo "[[Error]] FASTQ to FASTA conversion failed"; exit 1; }
    #            cmd="rm $fastq_file; rm $filtered_fastq"
    #            echo "$cmd"
    #            eval "$cmd"
    #        else
    #            fasta_out="${read_path}/${SRR}.fasta"
    #            cmd="seqkit fq2fa $fastq_file > $fasta_out"
    #            echo "$cmd"
    #            eval "$cmd" || { echo "[[Error]] FASTQ to FASTA conversion failed for CLR"; exit 1; }
    #            cmd="rm $fastq_file"
    #            echo "$cmd"
    #            eval "$cmd"
    #        fi
    
    #        continue
    #    fi
    #    echo "[$SRR] fasterq-dump failed or no FASTQ." >> "fail_alert/$data_no.$species.fail"
    #    exit 1
    #done
fi

print_duration $START_TIME


## Bam2Fasta (hifi.bam or CLR.bam)
b2f_merged=0
if [ -z "$(ls -A $read_path)" ]; then
    echo "[[Error]] No files found in $read_path"  >> "fail_alert/$data_no.$species.fail"
    exit 1
fi

if ls "$read_path" | grep -E '\.(fastq|fastq\.gz|fq|fq\.gz)$' > /dev/null; then # BUG finx 07022025
    echo "Found fastq. Skipping bam2fasta. Only merge"

    merge_fastq_cmd="zcat \"$read_path\"/*.fastq.gz \"$read_path\"/*.fq.gz | seqtk seq -A - > \"${read_path}${species}.fasta\""
    echo "$merge_fastq_cmd"
    eval "$merge_fastq_cmd"

    merge_fastq_cmd2="cat \"$read_path\"/*.fastq \"$read_path\"/*.fq | seqtk seq -A - >> \"${read_path}${species}.fasta\""
    echo "$merge_fastq_cmd2"
    eval "$merge_fastq_cmd2"
    fmerged=1
elif ls "$read_path" | grep -E '\.(fasta|fasta\.gz)$' > /dev/null; then
    echo "Found fasta. Merging fasta files..."

    merge_fasta_cmd="zcat \"$read_path\"/*.fasta.gz > \"${read_path}${species}.fasta\""
    echo "$merge_fasta_cmd"
    eval "$merge_fasta_cmd"

    merge_fasta_cmd2="cat \"$read_path\"/*.fasta >> \"${read_path}${species}.fasta\""
    echo "$merge_fasta_cmd2"
    eval "$merge_fasta_cmd2"
    fmerged=1
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
        fmerged=1

	if [ $? -ne 0 ]; then
            echo "[[Error]] bam2fasta failed on extracted hifi.bam" >> "fail_alert/$data_no.$species.fail"
            exit 1
        fi
    else
	echo "Not HiFi. Running bam2fasta directly on BAM files..."
        bam2fax_cmd="bam2fasta -u -j $cores -o ${read_path}${species} ${read_path}*.bam" # bam2fasta
        echo "$bam2fax_cmd"
        eval "$bam2fax_cmd"
        fmerged=1

	if [ $? -ne 0 ]; then
            echo "[[Error]] bam2fasta failed on bam files" >> "fail_alert/$data_no.$species.fail"
            exit 1
        fi
    fi
fi

print_duration $START_TIME


# split read ## 4 big file
if [ "$fmerged" -eq 1 ]; then 
    total=$(seqkit stats ${read_path}${species}.fasta | awk 'NR==2{print $4}' | tr -d ',')
    echo $total
    chunk_size=$(( (total + $splitread - 1) / $splitread ))
    echo $chunk_size
    mkdir -p ${read_path}${species}_split/
    for i in $(seq 1 $splitread); do
        start=$(( (i - 1) * chunk_size + 1 ))
        end=$(( i * chunk_size ))
        if [ $end -gt $total ]; then
            end=$total
        fi
        seqkit_cmd="seqkit range -r $start:$end ${read_path}${species}.fasta -o ${read_path}${species}_split/chunk_${i}.fasta"
        echo $seqkit_cmd
        eval $seqkit_cmd
    done


    if [ $? -eq 0 ]; then
        echo "Splitting done. Removing original FASTA file: ${read_path}${species}.fasta"
        rm -f "${read_path}${species}.fasta"
    else
        echo "[[Error]] split reads failed with merged fasta file" >> "fail_alert/$data_no.$species.fail"
        echo "Splitting failed. Original FASTA file retained."
        exit 1
    fi
fi

print_duration $START_TIME


## mapping reads
asm_path="1.data/$data_no.$species/ncbi_assembly/${species}"
fastx_path="1.data/$data_no.$species/reads/"
mapping_path="1.data/$data_no.$species/mapping/"

mkdir -p "$mapping_path"

#indexing_cmd="$minimap_exe -I 90G -d $asm_path.mmi ${asm_path}.fasta.gz"
#echo "$indexing_cmd"
#eval "$indexing_cmd"

shopt -s nullglob
read_files=(1.data/${data_no}.${species}/reads/${species}_split/*.fasta)

jobid_list=()

for read_file in "${read_files[@]}"; do
    base_name=$(basename "$read_file" .fasta)
    sam_out="${mapping_path}${base_name}.sam"
    log_out="${mapping_path}${base_name}.log"

    if [ "$folder" = "pacbio_hifi" ]; then
        map_type="map-hifi"
    elif [ "$folder" = "pacbio" ]; then
        map_type="map-pb"
    else
        echo "Error: Unexpected read type '$folder' it would be a critical error on data or parsing"
        exit 1
    fi
    
    sbatch_cmd="sbatch --job-name=map_${base_name} \
        --account=begy-delta-cpu \
        --partition=cpu \
        --cpus-per-task=$cores \
        --mem=180G \
        --time=48:00:00 \
        --output=$log_out \
        --wrap=\"$minimap_exe -x $map_type -a -t $cores $asm_path.mmi $read_file > $sam_out\""

    echo "$sbatch_cmd"
    submission_output=$(eval "$sbatch_cmd" 2>&1)
    echo "$submission_output"
    jobid=$(echo "$submission_output" | awk '/Submitted batch job/ {print $4}')
    if [[ -z "$jobid" || ! "$jobid" =~ ^[0-9]+$ ]]; then
        echo "[[Error]] sbatch failed or returned unexpected output for $read_file" >> "fail_alert/$data_no.$species.fail"
        echo "$submission_output" >> "fail_alert/$data_no.$species.fail"
    else
        jobid_list+=("$jobid")
    fi
done

if [ "${#jobid_list[@]}" -gt 0 ]; then
    dep_str=$(IFS=:; echo "${jobid_list[*]}")
    echo "Waiting for jobs to finish: $dep_str"

    wait_cmd="sbatch --dependency=afterok:$dep_str \
        --account=begy-delta-cpu \
        --partition=cpu \
        --time=00:10:00 \
        --mem=1G \
        --cpus-per-task=1 \
        --job-name=merge_${species} \
        --output=${mapping_path}${species}_merge.log \
        --wrap=\"echo 'All mappings done for $species'; touch ${mapping_path}${species}.done\""


    echo "$wait_cmd"
    eval "$wait_cmd"
else
    echo "[[Error]] No valid mapping jobs submitted for $species." >> "fail_alert/$data_no.$species.fail"
    exit 1
fi

timeout=129600  # 24 hours = 86400 seconds
elapsed=0

while [ ! -f "${mapping_path}${species}.done" ]; do
    echo "Waiting for mapping jobs to finish for $species..."
    sleep 60
    elapsed=$((elapsed + 60))
    if [ $elapsed -ge $timeout ]; then
        echo "[[Error]] Timeout waiting for .done file for $species" >> "fail_alert/$data_no.$species.fail"
        exit 1
    fi
done

echo "All mapping jobs done for $species. Proceeding."


print_duration $START_TIME


## sorting

sorted_bam_dir="${mapping_path}sorted_bams"
mkdir -p "$sorted_bam_dir"
sorted_bams=()
sorted_bam_file="${mapping_path}${species}.bam"  

for sam in "${mapping_path}"*.sam; do
    base=$(basename "$sam" .sam)
    bam="${sorted_bam_dir}/${base}.bam"

    echo "[INFO] Processing $sam → $bam"

    sort_cmd="sambamba view -t $cores -S -f bam \"$sam\" | sambamba sort -t $cores -l 1 -m ${memory}G -o \"$bam\" /dev/stdin"
    echo "$sort_cmd"
    eval "$sort_cmd"
    if [ $? -ne 0 ]; then
        echo "[[Error]] $sam failed" >> "fail_alert/$data_no.$species.fail"
        exit 1
    fi

    rm -f "$sam"
    sorted_bams+=("$bam")
done

merge_cmd="sambamba merge -t $cores \"$sorted_bam_file\" ${sorted_bams[*]}"
echo "$merge_cmd"
eval "$merge_cmd"

if [ $? -ne 0 ]; then
    echo "[[Error]] BAM merge failed for $species" >> "fail_alert/$data_no.$species.fail"
    exit 1
fi

echo "[INFO] Removing intermediate split BAM files for $species"
rm -f "${sorted_bams[@]}"


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
