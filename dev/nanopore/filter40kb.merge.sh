#filter reads > 40kb, merge files

#samtools view -h --threads 4 D15.xad.bam chrX | awk 'length($10) > 40000 || $1 ~ /^@/' | samtools view -bS - > X.40k.xad.bam
#sambamba index -t 4 X.40k.xad.bam X.40k.xad.bam.bai

samtools view -h --threads 4 D15.xac.bam chrX | awk 'length($10) > 40000 || $1 ~ /^@/' | samtools view -bS - > X.40k.xac.bam
sambamba index -t 4 X.40k.xac.bam X.40k.xac.bam.bai

samtools view -h --threads 4 D15.xag.bam chrX | awk 'length($10) > 40000 || $1 ~ /^@/' | samtools view -bS - > X.40k.xag.bam
sambamba index -t 4 X.40k.xag.bam X.40k.xag.bam.bai

samtools view -h --threads 4 D15.xae.bam chrX | awk 'length($10) > 40000 || $1 ~ /^@/' | samtools view -bS - > X.40k.xae.bam
sambamba index -t 4 X.40k.xae.bam X.40k.xae.bam.bai

samtools view -h --threads 4 D15.xah.bam chrX | awk 'length($10) > 40000 || $1 ~ /^@/' | samtools view -bS - > X.40k.xah.bam
sambamba index -t 4 X.40k.xah.bam X.40k.xah.bam.bai

samtools view -h --threads 4 D15.xaf.bam chrX | awk 'length($10) > 40000 || $1 ~ /^@/' | samtools view -bS - > X.40k.xaf.bam
sambamba index -t 4 X.40k.xaf.bam X.40k.xaf.bam.bai

samtools view -h --threads 4 D15.xab.bam chrX | awk 'length($10) > 40000 || $1 ~ /^@/' | samtools view -bS - > X.40k.xab.bam
sambamba index -t 4 X.40k.xab.bam X.40k.xab.bam.bai

samtools view -h --threads 4 D15.xai.bam chrX | awk 'length($10) > 40000 || $1 ~ /^@/' | samtools view -bS - > X.40k.xai.bam
sambamba index -t 4 X.40k.xai.bam X.40k.xai.bam.bai

sambamba merge -t 8 D15.X.40k.bam X.40k.xad.bam X.40k.xac.bam X.40k.xag.bam X.40k.xae.bam X.40k.xah.bam X.40k.xaf.bam X.40k.xab.bam X.40k.xai.bam

samtools view -h --threads 4 D15.X.40k.bam | awk 'length($10) > 60000 || $1 ~ /^@/' | samtools view -bS - > D15.X.60k.bam
sambamba index -t 4 D15.X.60k.bam D15.X.60k.bam.bai

sambamba merge -t 8 D15.bam D15.xab.bam D15.xac.bam D15.xad.bam D15.xag.bam D15.xae.bam D15.xah.bam D15.xaf.bam D15.xai.bam

samtools view -b D1172_001.bam chr3:27132710-27704513 chr10:69010502-69661371 -o D1172_001.mini.bam
sambamba index -t 1 D1172_001.mini.bam D1172_001.mini.bam.bai
samtools view -b --threads 4 D1186_01.bam chr3:27132710-27704513 chr10:69010502-69661371 -o D1186_01.mini.bam
sambamba index -t 1 D1186_01.mini.bam D1186_01.mini.bam.bai

