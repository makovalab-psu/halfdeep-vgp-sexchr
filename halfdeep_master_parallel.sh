#!/bin/bash

# 1. Get meta data (manually divided) information
# 2. Parsing parameters from the manual information
# 3. Run halfdeep_wrapper.sh with paramdeters

MAX_JOBS=2

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <metadata_file>"
    exit 1
fi

METADATA_FILE="$1"
sed -E ':a; s/\t\t/\t \t/g; ta' $METADATA_FILE > "${METADATA_FILE}.nd2space.txt"

while IFS=$'\t' read -r data_no order lineage scientific_name species english_name assembly_id main_hap accession sex sex_chromosomes size tech contigs check_num underbar_name check_gca fasta_name hifi_num hifi_vol hifi_files bam_num bam_vol bam_files hifiFQ_num hifiFQ_vol hifiFQ_files asm_num asm_vol asm_files read_on_ncbi fasta_manual steps additional_comment sex_chr sex_found_inmain read_db SRA_reads
do
    if [[ "$data_no" == "Data No." ]]; then
        continue
    fi

    [[ -z "$data_no" ]] && data_no="NULL"
    [[ -z "$order" ]] && order="NULL"
    [[ -z "$lineage" ]] && lineage="NULL"
    [[ -z "$scientific_name" ]] && scientific_name="NULL"
    [[ -z "$species" ]] && species="NULL"
    [[ -z "$english_name" ]] && english_name="NULL"
    [[ -z "$assembly_id" ]] && assembly_id="NULL"
    [[ -z "$main_hap" ]] && main_hap="NULL"
    [[ -z "$accession" ]] && accession="NULL"
    [[ -z "$sex" ]] && sex="NULL"
    [[ -z "$sex_chromosomes" ]] && sex_chromosomes="NULL"
    [[ -z "$size" ]] && size="NULL"
    [[ -z "$tech" ]] && tech="NULL"
    [[ -z "$contigs" ]] && contigs="NULL"
    [[ -z "$check_num" ]] && check_num="NULL"
    [[ -z "$underbar_name" ]] && underbar_name="NULL"
    [[ -z "$check_gca" ]] && check_gca="NULL"
    [[ -z "$fasta_name" ]] && fasta_name="NULL"
    [[ -z "$hifi_num" ]] && hifi_num="NULL"
    [[ -z "$hifi_vol" ]] && hifi_vol="NULL"
    [[ -z "$hifi_files" ]] && hifi_files="NULL"
    [[ -z "$bam_num" ]] && bam_num="NULL"
    [[ -z "$bam_vol" ]] && bam_vol="NULL"
    [[ -z "$bam_files" ]] && bam_files="NULL"
    [[ -z "$hifiFQ_num" ]] && hifiFQ_num="NULL"
    [[ -z "$hifiFQ_vol" ]] && hifiFQ_vol="NULL"
    [[ -z "$hifiFQ_files" ]] && hifiFQ_files="NULL"
    [[ -z "$asm_num" ]] && asm_num="NULL"
    [[ -z "$asm_vol" ]] && asm_vol="NULL"
    [[ -z "$asm_files" ]] && asm_files="NULL"
    [[ -z "$read_on_ncbi" ]] && read_on_ncbi="NULL"
    [[ -z "$fasta_manual" ]] && fasta_manual="NULL"
    [[ -z "$steps" ]] && steps="NULL"
    [[ -z "$additional_comment" ]] && additional_comment="NULL"
    [[ -z "$sex_chr" ]] && sex_chr="NULL"
    [[ -z "$sex_found_inmain" ]] && sex_found_inmain="NULL"
    [[ -z "$read_db" ]] && read_db="NULL"
    [[ -z "$SRA_reads" ]] && SRA_reads="NULL"

    split_array=($species)
    species="${split_array[0]} ${split_array[1]}"
    species=$(echo "$species" | sed 's/ /_/g')
    
    mkdir -p "1.data/"
    mkdir -p "2.log"
    
    LOGFILE="2.log/${data_no}.${species}.log"
    
    file_list=""
    if [[ "$read_db" == "SRA" ]]; then
        file_list="$SRA_reads"
    else
        if [[ "$hifiFQ_files" != "None" && -n "$hifiFQ_files" ]]; then
            file_list="$hifiFQ_files"
        elif [[ "$hifi_files" != "None" && -n "$hifi_files" ]]; then
            file_list="$hifi_files"
        else
            file_list="$bam_files"
        fi
    fi
    # Usage: $0 <data.No> <species name (2 words)> <accession> <ncbi asm file name> <VGP assembly id> <read DB (AWS/NCBI)> <read tech> <read file name list>"   
    #wrapper_cmd="sbatch -o \"$LOGFILE\" -e \"$LOGFILE\" halfdeep_sambamba.sh \"$data_no\" \"$species\" \"$accession\" \"$fasta_name\" \"$assembly_id\" \"$read_db\" \"$tech\" $file_list"
    
    
    wrapper_cmd="sbatch -o \"$LOGFILE\" -e \"$LOGFILE\" halfdeep_wrapper_sra_parallel.sh \"$data_no\" \"$species\" \"$accession\" \"$fasta_name\" \"$assembly_id\" \"$read_db\" \"$tech\" $file_list"
    
    #while [[ $(squeue -u bko --states=R,PD | grep halfdeep_w | wc -l) -ge 2 ]]; do
    #    echo "Waiting: halfdeep jobs >2"
    #    sleep 60
    #done

    while true; do
        count=$(squeue -u bko --states=R,PD 2>/dev/null | grep halfdeep_w | wc -l)

        if [[ -z "$count" || ! "$count" =~ ^[0-9]+$ ]]; then
            echo "squeue failed or returned unexpected output. Retrying in 60 sec..."
            sleep 60
            continue
        fi

        if [[ $count -ge 1 ]]; then
            echo "Waiting: halfdeep jobs >= 1 ($count running/pending)"
            sleep 3600
        else
            break
        fi
    done

    echo $wrapper_cmd
    jid=$(eval "$wrapper_cmd" | awk '{print $4}')
    
    echo $jid

    BKUP_LOGFILE="2.log/${data_no}.${species}_bkup.log"
    bkup_cmd="sbatch --dependency=afterok:$jid -o \"$BKUP_LOGFILE\" -e \"$BKUP_LOGFILE\" halfdeep_bkup_osn.sh \"$data_no\" \"$species\" \"$accession\" \"$fasta_name\" \"$assembly_id\" \"$read_db\" \"$tech\" $file_list"
    echo $bkup_cmd
    eval "$bkup_cmd"


done < "${METADATA_FILE}.nd2space.txt"

