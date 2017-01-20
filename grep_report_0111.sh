#!/bin/bash
#SBATCH -o grep_report_0111
#SBATCH -e grep_report_0111
#SBATCH -J grep_report_0111
#SBATCH --mail-type=all
#SBATCH --mail-user=kohsuanchen803@gmail.com
grep -Pn -i "elongation|woronin|ubiquitin|actin|tubulin|polyubiquitin|ribosylation factor|cyclophilin|cytochrome|histone|Hsp70|cell morphogenesis|EF hand|EF-hand|Ribosomal|Ribosom|Histidine|Elongation factor|rhodopsin|protein folding|RasGAP|ATP synthase|G-beta repeat|KH domain|Mitochondrial carrier|ribosylation|FKBP|Thioredoxin|ABC transporter|bZIP transcription factor|Cyclophilin|DEAH box |ATP synthase |Tyrosinase|NTF2|SAC3|HMG_box|Autophagy|Tubulin|RNA recognition motif|Histone|Rad60|leucine zipper|NAD_binding|ATPase activity|tyrosine|Rho-GDI|zinc-finger|Na_H_Exchanger|PfkB|unfolded protein binding|methyltransferase activity|proteasome|Av_adeno_fibre|Zinc knuckle|mitochondrial |RMI1_MOUSE|HUMAN|Transcription initiation factor|bioluminescence|profilin|tyrosine |transcription factor TFIIA |Ras family|AMP-binding|Rotamase|Zinc knuckle|PALI_NEUCR" report| cut -f 2 > trinotate_based_filter_Trinityname_0111

