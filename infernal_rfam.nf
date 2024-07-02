// Wrapping the process described here:
// https://docs.rfam.org/en/latest/genome-annotation.html
//
// See also:
// tRNAscan-SE for tRNAs
// RNAMMER for rRNA
// snoscan for snoRNAs
// SRPscan for SRP RNA

include { gt_tidysort } from './git-submodules/nf_utility/nf_utility.nf'

nextflow.enable.dsl=2

process download_rfam
{
conda 'bioconda::hmmer'
output:
 path "Rfam.cm"
 path "Rfam.clanin"

shell:
'''
wget ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.clanin
gunzip Rfam.cm.gz
'''
}

process cmpress_rfam
{
conda 'bioconda::infernal'
input:
 path "Rfam.cm"
output:
 path "Rfam.cm*"

shell:
'''
cmpress Rfam.cm
'''
}

process calc_genome_Z
{
conda 'bioconda::hmmer'
input:
 path genome
output:
 env INFERNAL_Z_SCORE
shell:
'''
resnum=$(esl-seqstat !{genome} | grep "Total # residues" | tr -s " " | cut -f 4 -d " ")
INFERNAL_Z_SCORE=$(awk "BEGIN {printf(\\"%.6f\\", $resnum*2/1000000)}") ## Since bash can't do floating point, using awk
'''
}

process cmscan_rfam
{
conda 'bioconda::hmmer'
publishDir "results/${task.process}", pattern: "", mode: 'link',overwrite:'true'
cache 'deep'
cpus 6
input:
 tuple path(genome),path(rfam_cm),path(rfam_cm_pressed1),path(rfam_cm_pressed2),path(rfam_cm_pressed3),path(rfam_cm_pressed4),path(rfam_clanin),val(infernal_z_score)
output:
 path("${genome}.tblout")
 path("${genome}.cmscan")
shell:
'''
##esl-seqstat !{genome} | grep "Total # residues" | tr -s " " | cut -f 4 -d " "
##(Total # of residues * 2) / 1000000 = 187.138164 ##For example

cmscan -Z !{infernal_z_score} --cpu !{task.cpus} --cut_ga --rfam --nohmmonly --tblout !{genome}.tblout --fmt 2 --clanin !{rfam_clanin} !{rfam_cm} !{genome} > !{genome}.cmscan
'''
}

process publish_results
{
publishDir "results/${task.process}", pattern: "", mode: 'link',overwrite:'true'
input:
 path merged_cmscan
 path merged_tbout
 path merged_gff
output:
 path "$merged_cmscan"
 path "$merged_tbout"
 path "$merged_gff"
shell:
'''
##Just a dummy process
'''
}

process download_infernal_tblout2gff
{
output:
 path "infernal-tblout2gff.pl"
shell:
'''
## See https://github.com/EddyRivasLab/infernal/issues/16
wget https://raw.githubusercontent.com/nawrockie/jiffy-infernal-hmmer-scripts/master/infernal-tblout2gff.pl
'''
}

process infernal_tblout2gff
{
publishDir "results/${task.process}", pattern: "", mode: 'link',overwrite:'true'
input:
 tuple path(perlscript),path(tblfile)
output:
 path "${tblfile}.gff"
shell:
'''
##perl infernal-tblout2gff.pl --cmscan !{tblfile} > !{tblfile}.gff
cat !{tblfile} | grep -v "#" | awk {'print $4 "\t"  $2 "\t" $3 "\t" $10 "\t" $11 "\t" $18 "\t" $12 "\t" "." "\t" "Clan="$6}' > !{tblfile}.gff
'''
}

workflow {
  genome_ch = channel.fromPath(params.genome)  
  calc_genome_Z(genome_ch)
  download_rfam()
  download_infernal_tblout2gff()
  cmpress_rfam(download_rfam.out[0])
  
  split_genome = genome_ch.splitFasta(size:'1.MB',file:true)
  split_cmds = split_genome.combine(download_rfam.out[0]).combine(cmpress_rfam.out.collect()).combine(download_rfam.out[1]).combine(calc_genome_Z.out)
  cmscan_rfam(split_cmds)

  tblout2gff_cmds = download_infernal_tblout2gff.out.combine(cmscan_rfam.out[0])
  infernal_tblout2gff(tblout2gff_cmds)

  cmscan_ch = cmscan_rfam.out[1].collectFile(keepHeader:true,skip:16,name:params.genome+"_merged.cmscan")
  tblout_ch = cmscan_rfam.out[0].collectFile(name:params.genome+"_merged.tblout")
  gff_ch = infernal_tblout2gff.out.collectFile(name:params.genome+"_merged.gff")
  gt_tidysort(gff_ch)
  publish_results(cmscan_ch,tblout_ch,gff_ch)
}
