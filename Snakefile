import os

DS_NAME = "autism_noncoding"
REGIONS = "probe_regions.bed"
PARAM_FILE = "params.txt"

REF = "/net/eichler/vol2/eee_shared/assemblies/human_1kg_v37/human_1kg_v37.fasta"
OUTDIR = "GATK_depth"

SAMPLE_DIR = "/net/eichler/vol21/projects/autism_noncoding_capture/nobackups/analysis/bam"

XHMM_PATH = "/net/eichler/vol4/home/tychele/software/xhmm/statgen-xhmm-3c57d886bc96/xhmm"

SAMPLES = {}

for file in os.listdir(SAMPLE_DIR):
    if file.endswith(".bam"):
        sn = file.replace(".final.sort.bam", "")
        SAMPLES[sn] = "%s/%s" % (SAMPLE_DIR, file)

def get_bam_from_sample(wildcards):
    return SAMPLES[wildcards.sample]

rule all:
    input: expand("%s/{sample}.DATA.{ext}" % OUTDIR, sample = SAMPLES.keys(), ext = ["sample_interval_summary", "sample_interval_statistics"])

rule discover_cnvs:
    input: "%s.PCA_norm.txt" % DS_NAME, "autism_noncoding.same_filtered.RD.txt"
    output: "%s.xcnv" % DS_NAME, "%s.aux_xcnv" % DS_NAME
    params: sge_opts = "-l mfree=20G"
    shell:
        "{XHMM_PATH} --discover -p {PARAM_FILE} -r {input[0]} -R {input[1]} -c {output[0]} -a {output[1]} -s {DS_NAME}"

rule filter_rd_data:
    input: "merged_GATK_depths.txt", 
           "%s.filtered_centered.RD.filtered_targets.txt" % DS_NAME, "%s.PCA_norm.filt_center.excluded_targets.txt" % DS_NAME,
           "%s.filtered_centered.RD.filtered_samples.txt" % DS_NAME, "%s.PCA_norm.filt_center.excluded_samples.txt" % DS_NAME
    output: "autism_noncoding.same_filtered.RD.txt"
    params: sge_opts = "-l mfree=20G"
    shell:
        "{XHMM_PATH} --matrix -r {input[0]} --excludeTargets {input[1]} --excludeTargets {input[2]} --excludeSamples {input[3]} --excludeSamples {input[4]} -o {output}"

rule pca_filt_center:
    input: "%s.PCA_norm.txt" % DS_NAME
    output: "%s.PCA_norm.filt_center.txt" % DS_NAME, "%s.PCA_norm.filt_center.excluded_targets.txt" % DS_NAME, "%s.PCA_norm.filt_center.excluded_samples.txt" % DS_NAME
    params: sge_opts = "-l mfree=20G"
    shell:
        "{XHMM_PATH} --matrix -r {input[0]}  --centerData --centerType sample --zScoreData -o {output[0]} "
        "--outputExcludedTargets {output[1]} --outputExcludedSamples {output[2]} --maxSdTargetRD 30"

rule pca_normalize:
    input: "%s.filtered_centered.RD.txt" % DS_NAME, "%s.RD_PCA" % DS_NAME
    output: "%s.PCA_norm.txt" % DS_NAME
    params: sge_opts = "-l mfree=20G"
    shell:
        "{XHMM_PATH} --normalize -r {input[0]} --PCAfiles {input[1]} --normalizeOutput {output[0]} --PCnormalizeMethod PVE_mean --PVE_mean_factor 0.7"

rule run_pca:
    input: "%s.filtered_centered.RD.txt" % DS_NAME
    output: "%s.RD_PCA" % DS_NAME
    params: sge_opts = "-l mfree=20G"
    shell:
        "{XHMM_PATH} --PCA -r {input} --PCAfiles {output}"

rule filt_center:
    input: "merged_GATK_depths.txt"
    output: "%s.filtered_centered.RD.txt" % DS_NAME, "%s.filtered_centered.RD.filtered_targets.txt" % DS_NAME, "%s.filtered_centered.RD.filtered_samples.txt" % DS_NAME
    params: sge_opts = "-l mfree=8G"
    shell:
        "{XHMM_PATH} --matrix -r {input[0]} --centerData --centerType target -o {output[0]} "
        "--outputExcludedTargets {output[1]} --outputExcludedSamples {output[2]} "
        "--minTargetSize 10 --maxTargetSize 10000 --minMeanTargetRD 10 --maxMeanTargetRD 500 --minMeanSampleRD 25 --maxMeanSampleRD 200 --maxSdSampleRD 150"

rule merge_GATK_depths:
    input: "GATK_depths.txt"
    output: "merged_GATK_depths.txt"
    params: sge_opts = "-l mfree=12G"
    shell:
        "{XHMM_PATH} --mergeGATKdepths -o {output} --GATKdepthsList {input}"

rule make_GATK_depth_list:
    input: expand("%s/{sample}.DATA.sample_interval_summary" % OUTDIR, sample = SAMPLES.keys())
    output: "GATK_depths.txt"
    params: sge_opts = ""
    run:
        with open(output[0], "w") as writer:
            for file in input:
                writer.write(file + "\n")

rule get_GATK_depth:
    input: get_bam_from_sample
    output: "%s/{sample}.DATA.sample_interval_summary" % OUTDIR, "%s/{sample}.DATA.sample_interval_statistics" % OUTDIR
    params: sge_opts = "-l mfree=21G -N GATK_{sample}", sample = "{sample}", prefix = "%s/{sample}.DATA" % OUTDIR, tmpdir = "$TMPDIR"
    shell:
        "module load java/7u17 GATK/3.3-0; "
        "java -Xmx10G -jar $GATK_DIR/GenomeAnalysisTK.jar -T DepthOfCoverage -I {input} -L {REGIONS} -R {REF} "
        "-dt BY_SAMPLE -dcov 5000 -l INFO --omitDepthOutputAtEachBase --omitLocusTable --minBaseQuality 0 --minMappingQuality 20 "
        "--start 1 --stop 5000 --nBins 200 --includeRefNSites --countType COUNT_FRAGMENTS -o {params.prefix}"
