This workflow is based on:  

bwa/0.7.12 on 1000G phase II GRCh37 genome
GATK/3.5-0 with Haplotype Caller and hard filters (using GATK and bcbio best practices for exome)

NISC_laneBam_bwaRealign_callGVCF.sh is the wrapper script. It is fed a matrix file*, the sbatch job name, and the location of the exome bait bed file (specific to certain library prep kits). 

The script calls ~/bin/exome_workflow_v02/create_scp_and_sbatch_jobs_for_NISC_laneBams.py which creates scripts that are executed to scp the NISC lane bam files from Trek (1), process thems with BWA (2), and creates the script to call raw genotypes (1, 3). 

1. create_scp_and_sbatch_jobs_for_NISC_laneBams.py
2. realign_NISC_laneBams_with_bwa.py
3. process_and_callGVCF.sh 

Afterwards we have GVCF (raw genonotype) files for each exome sample. The next step is to run these in together, giving GVCF_to_hardFilteredVCF.sh the GVCFs in a \n separated list of the GVCFs you want to call together. 

Then you have a single VCF containing the genotypes for your cohort. The next steps are to run the Gemini workflow to prioritize variants for examination. 





* The matrix file contains the 'common' and the lane bam file name, which can be used to find the exact file on Trek

```bash
sqlite3 ~/NISC_laneBam.sqlite3.db "SELECT NISC_LaneBams.Sample, LaneBam_File FROM NISC_LaneBams INNER JOIN Sample_Info ON NISC_LaneBams.Sample=Sample_Info.Sample WHERE Sample_Info.Project='DDL' AND Sample_Info.DateAdded='2016-05-04'" | sort | cut -d"|" -f1,2,3 --output-delimiter=' '
1232 160229_OPTIMUS_C855NANXX.8.11961896
1232 160311_OPTIMUS_C8E3LANXX.3.11961896
1232 160315_YOSHI_C8E8CANXX.1.11961896
1232 160315_YOSHI_C8E8CANXX.2.11961896
1233 160229_OPTIMUS_C855NANXX.8.11961899
1233 160311_OPTIMUS_C8E3LANXX.3.11961899
1233 160315_YOSHI_C8E8CANXX.1.11961899
1233 160315_YOSHI_C8E8CANXX.2.11961899
1238 160229_OPTIMUS_C855NANXX.8.11961880
1238 160311_OPTIMUS_C8E3LANXX.3.11961880
1238 160315_YOSHI_C8E8CANXX.1.11961880
1238 160315_YOSHI_C8E8CANXX.2.11961880
1313 160324_YOSHI_C8RWMANXX.1.12099201
1313 160324_YOSHI_C8RWMANXX.2.12099201
1313 160325_OPTIMUS_C8RT6ANXX.1.12099201
1313 160325_OPTIMUS_C8RT6ANXX.2.12099201
1313 160325_OPTIMUS_C8RT6ANXX.3.12099201
1484 160229_OPTIMUS_C855NANXX.8.11961900
1484 160311_OPTIMUS_C8E3LANXX.3.11961900
1484 160315_YOSHI_C8E8CANXX.1.11961900
1484 160315_YOSHI_C8E8CANXX.2.11961900
1501 160413_ROSIE_C8T1PANXX.4.12173415
1501 160416_OPTIMUS_C8RT2ANXX.1.12173415
1501 160416_OPTIMUS_C8RT2ANXX.2.12173415
1501 160416_OPTIMUS_C8RT2ANXX.3.12173415
1522 160413_ROSIE_C8T1PANXX.4.12173419
1522 160416_OPTIMUS_C8RT2ANXX.1.12173419
1522 160416_OPTIMUS_C8RT2ANXX.2.12173419
1522 160416_OPTIMUS_C8RT2ANXX.3.12173419
LV111315 160229_OPTIMUS_C855NANXX.8.11961883
LV111315 160311_OPTIMUS_C8E3LANXX.3.11961883
LV111315 160315_YOSHI_C8E8CANXX.1.11961883
LV111315 160315_YOSHI_C8E8CANXX.2.11961883
```
