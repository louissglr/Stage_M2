#!/bin/bash

#SBATCH --job-name=star_solo          # Nom du job
#SBATCH --comment="STARsolo" # Description
#SBATCH --output=%x_%j.out            # Fichier de sortie
#SBATCH --error=%x_%j.err             # Fichier d'erreur
#SBATCH --time=0-24:00:00              # Temps limite
#SBATCH --nodes=1                     # Nombre de nœuds
#SBATCH --ntasks=1                     # Nombre de tâches
#SBATCH --cpus-per-task=25              # Nombre de CPU
#SBATCH --mem=50G                      # Mémoire totale
#SBATCH --qos=short                    # Priorité

# Ajout des alertes par mail
#SBATCH --mail-type=END,FAIL   # Alerte lorsque le job se termine ou échoue
#SBATCH --mail-user=louis.sanglier@univ-nantes.fr  # Adresse e-mail pour les notifications

# Paths
STAR_EXEC="/home/sanglier-l@univ-nantes.fr/software/STAR/bin/Linux_x86_64_static/STAR"
GENOME_DIR="/home/sanglier-l@univ-nantes.fr/Align/ref2/ref_mm10_td_sv/index_mm10_td_sv"
DATA_DIR="/home/sanglier-l@univ-nantes.fr/data/scrnaseq"
OUTPUT_DIR="/scratch/nautilus/users/sanglier-l@univ-nantes.fr/align/STARsolo_all"

# Création du répertoire output
mkdir -p $OUTPUT_DIR
chmod u+x $OUTPUT_DIR

# Boucle sur les fichiers R1 pour trouver les paires
for R1_FILE in ${DATA_DIR}/*_R1_001.fastq.gz; do
    SAMPLE_NAME=$(basename $R1_FILE _R1_001.fastq.gz)
    R2_FILE="${DATA_DIR}/${SAMPLE_NAME}_R2_001.fastq.gz"

    if [[ -f "$R2_FILE" ]]; then
        SAMPLE_OUT_DIR="${OUTPUT_DIR}/${SAMPLE_NAME}"
        mkdir -p $SAMPLE_OUT_DIR

        # STARsolo
        $STAR_EXEC \
            --genomeDir $GENOME_DIR \
            --readFilesIn $R2_FILE $R1_FILE \
            --readFilesCommand zcat \
            --runThreadN 25 \
            --soloType CB_UMI_Simple \
            --soloCBwhitelist None \
            --soloBarcodeReadLength 0 \
            --soloCBstart 1 --soloCBlen 12 \
            --soloUMIstart 13 --soloUMIlen 9 \
            --soloFeatures Gene \
            --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
            --soloCellFilter EmptyDrops_CR \
            --outSAMattributes CR UR CY UY CB UB \
            --outSAMtype BAM SortedByCoordinate \
            --soloMultiMappers EM \
            --soloUMIdedup 1MM_CR \
            --soloUMIfiltering MultiGeneUMI_CR \
            --outFileNamePrefix ${SAMPLE_OUT_DIR}/
    else
        echo "Fichier R2 manquant pour $SAMPLE_NAME, skipping."
    fi
done

echo "Alignement STARsolo terminé pour tous les échantillons."
