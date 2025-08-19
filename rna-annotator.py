import subprocess
import os
import sys
import argparse
import pandas as pd
import time
import requests
import csv 
import gzip
import tempfile 

# --- IGV Configuration Constants ---
BED_DIR = "results/bed_tracks"
IGV_PORT = 60151

# --- Core Analysis Functions ---

def create_directories():
    results_dir = "results"
    bed_tracks_dir = os.path.join(results_dir, "bed_tracks")
    detailed_results_dir = os.path.join(results_dir, "detailed_results")

    if not os.path.exists(results_dir):
        os.makedirs(results_dir)
    if not os.path.exists(bed_tracks_dir):
        os.makedirs(bed_tracks_dir)
    if not os.path.exists(detailed_results_dir):
        os.makedirs(detailed_results_dir)

    return bed_tracks_dir, detailed_results_dir

# --- HELPER FUNCTION FOR CROSS-PLATFORM BEDTOOLS ---
def run_bedtools_intersect(input_file, chrom, start, end):
    """A robust, cross-platform way to run bedtools intersect."""
    temp_region_filename = None
    try:
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.bed') as temp_region_file:
            temp_region_file.write(f"{chrom}\t{start}\t{end}\n")
            temp_region_filename = temp_region_file.name
        
        cmd = ["bedtools", "intersect", "-a", input_file, "-b", temp_region_filename, "-wa"]
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        return result
    finally:
        if temp_region_filename and os.path.exists(temp_region_filename):
            os.remove(temp_region_filename)

def extract_functional_elements(chrom, start, end, bed_files, bed_tracks_dir, detailed_results_dir):
    merged_bed9_file = os.path.join(bed_tracks_dir, f"RefSeq_Functional_Elements__{chrom}_{start}_{end}.bed")
   
    for bed_file in bed_files:
        if not os.path.exists(bed_file):
            print(f"❌ Error: File not found: {bed_file}")
            continue
        try:
           
            result = run_bedtools_intersect(bed_file, chrom, start, end)
            
            if not result.stdout.strip():
                print(f"⚠️ No data found in {bed_file} for the specified region.")
                continue 
            
            
            with open(merged_bed9_file, "a") as merged_bed_out:
                tsv_output = os.path.join(detailed_results_dir, f"{os.path.basename(bed_file).replace('.bed', f'_{chrom}_{start}_{end}.tsv')}")
                
                with open(tsv_output, "w") as tsv_out:
                    for line in result.stdout.strip().split("\n"):
                        cols = line.split("\t")
                        if len(cols) < 3: continue
                        tsv_out.write("\t".join(cols) + "\n")
                        
                        if os.path.basename(bed_file) == "FEbiolregions_AR110_GRCh38.p14.bed" and len(cols) > 10:
                            name = cols[10]
                        elif len(cols) > 3:
                            name = cols[3]
                        else:
                            name = "."
                        
                        chrom_out, start_out, end_out = cols[0], cols[1], cols[2]
                        bed9_line = f"{chrom_out}\t{start_out}\t{end_out}\t{name}\t0\t.\t{start_out}\t{end_out}\t0,0,255\n"
                        merged_bed_out.write(bed9_line)
                print(f"📄 TSV written: {tsv_output}")
        except subprocess.CalledProcessError as e:
            print(f"❌ Error processing {bed_file}:\n{e.stderr.strip()}")
    print(f"\n✅ ✅ Combined BED9 track ready: {merged_bed9_file}")

def process_eclips_data(chrom, start, end, bed_tracks_dir, detailed_results_dir):
    eclips_data = "resources_data_sets/eclips_data.tsv"
    print(f"Processing eclips data from {eclips_data}...")
    
    eclips_bed9_file = os.path.join(bed_tracks_dir, f"eclips__{chrom}_{start}_{end}.bed")
    try:
        
        result = run_bedtools_intersect(eclips_data, chrom, start, end)

        if not result.stdout.strip():
            print("⚠️ No eclips_data found in this region.")
            return
        
        tsv_output = os.path.join(detailed_results_dir, f"eclips_{chrom}_{start}_{end}.tsv")
        with open(tsv_output, "w") as tsv_out, open(eclips_bed9_file, "w") as merged_bed_out:
            for line in result.stdout.strip().split("\n"):
                cols = line.split("\t")
                tsv_out.write("\t".join(cols) + "\n")
                chrom_out, start_out, end_out = cols[0], cols[1], cols[2]
                name = cols[10]
                bed9_line = f"{chrom_out}\t{start_out}\t{end_out}\t{name}\t0\t.\t{start_out}\t{end_out}\t225,0,0\n"
                merged_bed_out.write(bed9_line)
        print(f"✅ BED track created: {eclips_bed9_file}")
    except subprocess.CalledProcessError as e:
        print(f"❌ Error processing {eclips_data}:\n{e.stderr.strip()}")
        
        
ENSEMBL_API = "https://rest.ensembl.org/overlap/region/human"

def fetch_snps(chrom, start, end,bed_tracks_dir, detailed_results_dir):
    
    region = f"{chrom}:{start}-{end}"
    url = f"{ENSEMBL_API}/{region}?feature=variation"
    headers = {"Content-Type": "application/json"}
    response = requests.get(url, headers=headers)
    if not response.ok:
        raise Exception(f"API error: {response.status_code} - {response.text}")
    snp_data = response.json()
    if not snp_data:
        print("⚠️ No SNPs found in this region.")
        return
    out_base = f"snp_annotations_{chrom}_{start}_{end}"
    csv_file = os.path.join(detailed_results_dir, f"{out_base}.csv")
    bed_file = os.path.join(bed_tracks_dir,f"{out_base}.bed")
    df = pd.json_normalize(snp_data)
    df.to_csv(csv_file, index=False)
    print(f"✅ CSV written: {csv_file}")
    with open(bed_file, "w") as bed:
        for snp in snp_data:
            chrom_out = f"chr{snp['seq_region_name']}"
            bed_start = snp['start'] - 1
            bed_end = snp['end']
            alleles = snp.get('alleles', [])
            allele_str = "/".join(alleles) if isinstance(alleles, list) else str(alleles)
            name = f"{allele_str}"
            bed.write(f"{chrom_out}\t{bed_start}\t{bed_end}\t{name}\t0\t.\t{bed_start}\t{bed_end}\t255,165,0\n")
    print(f"✅ BED written: {bed_file}")

def fetch_miRNA(chrom, start, end, bed_tracks_dir, detailed_results_dir):
    try:
        
        result = run_bedtools_intersect("resources_data_sets/hsa.gff3", chrom, start, end)

        mirna_tsv_output = os.path.join(detailed_results_dir, f"miRNA_{chrom}_{start}_{end}.tsv")
        mirna_bed9_file = os.path.join(bed_tracks_dir, f"miRNA_{chrom}_{start}_{end}.bed")
        if not result.stdout.strip():
            print("⚠️ No miRNA found in this region.")
            return
        
        with open(mirna_tsv_output, "w") as mirna_tsv , open(mirna_bed9_file, "w") as mirna_bed_out:
            for line in result.stdout.strip().split("\n"):
                cols = line.split("\t")
                mirna_tsv.write("\t".join(cols) + "\n")
                if len(cols) < 9:
                    continue
                chrom_out = cols[0]
                start_out = int(cols[3]) - 1
                end_out = int(cols[4]) 
                name= cols[8].split(";")[2].split("=")[1]
                bed9_line = f"{chrom_out}\t{start_out}\t{end_out}\t{name}\t0\t.\t{start_out}\t{end_out}\t0,128,128\n" # Removed typo .t
                mirna_bed_out.write(bed9_line)
        print(f"✅ BED track created: {mirna_bed9_file}")
    except subprocess.CalledProcessError as e:
        print(f"❌ Error processing hsa.gff3 file:\n{e.stderr.strip()}")

def fetch_chemmodi(chrom, start, end, bed_tracks_dir, detailed_results_dir):
    try:
       
        result = run_bedtools_intersect("resources_data_sets/All_Modifications_hg38.csv", chrom, start, end)

        chemmodi_tsv_output = os.path.join(detailed_results_dir, f"chemical_modifications_{chrom}_{start}_{end}.tsv")
        chemmodi_bed9_file = os.path.join(bed_tracks_dir, f"chemical_modifications_{chrom}_{start}_{end}.bed")
        header_string = "chromosome\tmodStart\tmodEnd\tmodId\tscore\tstrand\tmodName\tmodType\tsupportNum\tsupportList\tpubmedIds\tgeneName\tgeneType\tregion\tsequence\tmotif_score\n"
        if not result.stdout.strip():
            print("⚠️ No chemical modifications found in this region.")
            return
    
        with open(chemmodi_tsv_output, "w") as chemmodi_tsv , open(chemmodi_bed9_file, "w") as chemmodi_bed_out:
            chemmodi_tsv.write(header_string)
            for line in result.stdout.strip().split("\n"):
                cols = line.split("\t")
                chemmodi_tsv.write("\t".join(cols) + "\n")
                if len(cols) < 9:
                    continue
                chrom_out = cols[0]
                start_out = int(cols[1]) - 1
                end_out = int(cols[2]) 
                name = cols[7] if cols[3].startswith("otherMod_site") else cols[3]
                bed9_line = f"{chrom_out}\t{start_out}\t{end_out}\t{name}\t0\t.\t{start_out}\t{end_out}\t0,0,225\n"
                chemmodi_bed_out.write(bed9_line)
        print(f"✅ BED track created: {chemmodi_bed9_file}")
    except subprocess.CalledProcessError as e:
        print(f"❌ Error processing chemical modifications file:\n{e.stderr.strip()}")

def polyA(chrom, start, end, bed_tracks_dir, detailed_results_dir):
    try:
        
        result = run_bedtools_intersect("resources_data_sets/polyA._sites.bed", chrom, start, end)
        
        polyA_tsv_output = os.path.join(detailed_results_dir, f"polyA_sites_{chrom}_{start}_{end}.tsv")
        polyA_bed9_file = os.path.join(bed_tracks_dir, f"polyA_sites_{chrom}_{start}_{end}.bed")
        header_string = "chromosome\tStart\tEnd\tunique cluster ID\taverage expression\tstrand\tpercentage of samples that support the particular cluster\tnumber of different 3' end sequencing protocols that support the particular cluster\taverage expression (tags per million, tpm) across all samples\tcluster annotation\tmore information\n"
        if not result.stdout.strip():
            print("⚠️ No polyA sites found in this region.")
            return
        with open(polyA_tsv_output, "w") as polyA_tsv , open(polyA_bed9_file, "w") as polyA_bed_out:
            polyA_tsv.write(header_string)    
            for line in result.stdout.strip().split("\n"):
                cols = line.split("\t")
                polyA_tsv.write("\t".join(cols) + "\n")
                if len(cols) < 3:
                    continue
                chrom_out = cols[0]
                start_out = int(cols[1]) - 1
                end_out = int(cols[2])
                name_map = {"TE": "terminal exon", "EX": "exonic", "IN": "intronic", "DS": "1,000 nt downstream of an annotated terminal exon", "AE": "anti-sense to exon", "AI": "anti-sense to intron", "AU": "1,000 nt upstream in anti-sense direction of a transcription start site", "IG": "intergenic"}
                code = cols[9]
                name = name_map.get(code, code)
                bed9_line = f"{chrom_out}\t{start_out}\t{end_out}\t{name}\t0\t.\t{start_out}\t{end_out}\t0,128,0\n"
                polyA_bed_out.write(bed9_line)
        print(f"✅ BED track created: {polyA_bed9_file}")
    except subprocess.CalledProcessError as e:
        print(f"❌ Error processing polyA_sites.bed file:\n{e.stderr.strip()}")     

def repeated_element(chrom, start, end, bed_tracks_dir, detailed_results_dir):
    
    Dfam_API = "https://dfam.org/api/annotations"
    params = {"assembly": "hg38", "chrom": chrom, "start": start, "end": end, "nrph": "true"}
    response = requests.get(Dfam_API, params=params)
    if response.status_code == 200:
        data = response.json()
        hits = data.get('hits', [])
        if not hits:
            print("⚠️ No repeated elements found in this region.")
            return
        out_base = f"pfam_repeted_elements_{chrom}_{start}_{end}"
        csv_file = os.path.join(detailed_results_dir, f"{out_base}.csv")
        bed_file = os.path.join(bed_tracks_dir,f"{out_base}.bed")
        df = pd.json_normalize(hits)
        df.to_csv(csv_file, index=False)
        print(f"✅ CSV written: {csv_file}")
        with open(bed_file, "w") as bed:
            for snp in hits:
                chrom_out = snp["sequence"]
                bed_start = int(snp['ali_start']) - 1  if snp['strand'] == "+" else int(snp['ali_end']) - 1
                bed_end = int(snp['ali_end']) if snp['strand'] == "+" else int(snp['ali_start'])
                name = snp["type"]
                bed.write(f"{chrom_out}\t{bed_start}\t{bed_end}\t{name}\t0\t.\t{bed_start}\t{bed_end}\t165,42,42\n")
        print(f"✅ BED written: {bed_file}")

def chemical_prop(chrom, start, end, bed_tracks_dir, detailed_results_dir, chemical_prop_file):
    for chem_file in chemical_prop_file:
        merged_bed9_file = os.path.join(bed_tracks_dir, f"{os.path.basename(chem_file).replace('.wig.gz', f'_{chrom}_{start}_{end}.wig')}")
        if not os.path.exists(chem_file):
            print(f"❌ Error: File not found: {chem_file}")
            continue
        try:
           
            result = run_bedtools_intersect(chem_file, chrom, start, end)

            if not result.stdout.strip():
                print(f"⚠️ No data found in {chem_file} for the specified region.")
                continue
            with open(merged_bed9_file, "w") as merged_bed_out:
                with gzip.open(chem_file, "rt") as f:
                    header_line = f.readline()
                merged_bed_out.write(header_line)
                merged_bed_out.write(result.stdout)
                tsv_output = os.path.join(detailed_results_dir, f"{os.path.basename(chem_file).replace('.wig.gz', f'_{chrom}_{start}_{end}.tsv')}")
                with open(tsv_output, "w") as tsv_out:
                    for line in result.stdout.strip().split("\n"):
                        cols = line.split("\t")
                        if len(cols) < 3: continue
                        tsv_out.write("\t".join(cols) + "\n")
                    print(f"📄 TSV written: {tsv_output}")
                print(f"📄 TSV written: {tsv_output}")
        except subprocess.CalledProcessError as e:
            print(f"❌ Error processing {chem_file}:\n{e.stderr.strip()}")
        print(f"\n✅ ✅ Combined WIG track ready: {merged_bed9_file}")

def clinvar(chrom, start, end, bed_tracks_dir, detailed_results_dir, clinvar_file):
    try:
        chrom_for_bedtools = chrom.removeprefix('chr')
        
       
        result = run_bedtools_intersect(clinvar_file, chrom_for_bedtools, start, end)
        
        clinvar_tsv_output = os.path.join(detailed_results_dir, f"clinvar_{chrom}_{start}_{end}.tsv")
        clinvar_bed9_file = os.path.join(bed_tracks_dir, f"clinvar_{chrom}_{start}_{end}.bed")
        if not result.stdout.strip():
            print("⚠️ No clinvar data found in this region.")
            return
        with open(clinvar_tsv_output, "w") as clinvar_tsv, open(clinvar_bed9_file, "w") as clinvar_bed_out:
            for line in result.stdout.strip().split("\n"):
                if line.startswith('#'): continue
                cols = line.split("\t")
                clinvar_tsv.write("\t".join(cols) + "\n")
                chrom_out = cols[0] if cols[0].startswith('chr') else f"chr{cols[0]}"
                start_out = int(cols[1]) - 1
                ref_allele = cols[3]
                end_out = start_out + len(ref_allele)
                significance = 'N/A'
                info_parts = cols[7].split(';')
                for part in info_parts:
                    if part.startswith("CLNSIG="):
                        significance = part.split('=')[1]
                        break
                name = f"{cols[3]}/{cols[4]}:{significance}"
                bed9_line = f"{chrom_out}\t{start_out}\t{end_out}\t{name}\t0\t.\t{start_out}\t{end_out}\t0,0,0\n" 
                clinvar_bed_out.write(bed9_line)
        print(f"✅ BED track created: {clinvar_bed9_file}")
    except subprocess.CalledProcessError as e:
        print(f"❌ Error processing clinvar file:\n{e.stderr.strip()}") 

def target_scan(chrom, start, end, bed_tracks_dir, detailed_results_dir, target_scan_file):
    merged_bed9_file = os.path.join(bed_tracks_dir,  f'Target_scan_miRNA_binding_prediction_{chrom}_{start}_{end}.bed')
    with open(merged_bed9_file, "a") as merged_bed_out: 
        for target_file in target_scan_file:
            if not os.path.exists(target_file):
                print(f"❌ Error: File not found: {target_file}")
                continue
            try:
                
                result = run_bedtools_intersect(target_file, chrom, start, end)

                if not result.stdout.strip():
                    print(f"⚠️ No data found in {target_file} for the specified region.")
                    continue
                merged_bed_out.write(result.stdout)
                tsv_output = os.path.join(detailed_results_dir, f"{os.path.basename(target_file).replace('.bed.gz', f'_{chrom}_{start}_{end}.tsv')}")
                with open(tsv_output, "w") as tsv_out:
                    tsv_out.write(result.stdout)
                print(f"📄 TSV written: {tsv_output}")
            except subprocess.CalledProcessError as e:
                print(f"❌ Error processing {target_file}:\n{e.stderr.strip()}")
    print(f"\n✅ ✅ Combined BED track ready: {merged_bed9_file}")

def phastCons(chrom, start, end, bed_tracks_dir, detailed_results_dir, phastCons_file):
    output_bedgraph_file = os.path.join(bed_tracks_dir, f"phastCons_{chrom}_{start}_{end}.bedGraph")
    index_file = f"{phastCons_file}.tbi"
    if not os.path.exists(phastCons_file):
        print(f"❌ Error: Data file not found: {phastCons_file}")
        return
    if not os.path.exists(index_file):
        print(f"❌ FATAL ERROR: Index file not found: {index_file}")
        print(f"   -> Please index the file with: tabix -p bed {phastCons_file}")
        return
    try:
        region = f"{chrom}:{start}-{end}"
        
        cmd = f"tabix {phastCons_file} {region}"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True, check=True)
        if not result.stdout.strip():
            print(f"⚠️ No phastCons data found in this region.")
            return
        with open(output_bedgraph_file, "w") as bedgraph_out:
            track_header = f'track type=bedGraph name="Coevolution Conservation Score" description="Vertebrate Conservation" visibility=full graphType=heatmap autoScale=on color=50,150,50\n'
            bedgraph_out.write(track_header)
            bedgraph_out.write(result.stdout)
        print(f"✅ bedGraph track created: {output_bedgraph_file}")
        tsv_output = os.path.join(detailed_results_dir, f"phastCons_{chrom}_{start}_{end}.tsv")
        with open(tsv_output, "w") as tsv_out:
            tsv_out.write(result.stdout)
        print(f"📄 TSV written: {tsv_output}")
    except subprocess.CalledProcessError as e:
        print(f"❌ Error processing {phastCons_file} with tabix:\n{e.stderr.strip()}")

# --- IGV Launcher Functions ---

def launch_igv(igv_app_path):
    """
    Launches the IGV application in a cross-platform way.
    """
    print("🚀 Launching IGV...")
    if sys.platform == "win32":
        print(f"  (Detected Windows OS)")
        os.startfile(igv_app_path)
    elif sys.platform == "darwin":
        print(f"  (Detected macOS)")
        subprocess.Popen(["open", "-a", igv_app_path])
    else: # Linux
        print(f"  (Detected Linux OS)")
        subprocess.Popen([igv_app_path])
    print("⏳ Waiting 15 seconds for IGV to start and listen on its port...")
    time.sleep(15)

def send_to_igv(file_path, region=None):

    if sys.platform == "win32":
        print(f"⚠️ Warning: Automatic loading in IGV with 'nc' is not supported on Windows.")
        print(f"   Please load this file manually in IGV: {file_path}")
        return

    abs_path = os.path.abspath(file_path)
    load_cmd = f'echo \'load "{abs_path}"\' | nc localhost {IGV_PORT}'
    try:
        subprocess.run(load_cmd, shell=True, check=True, capture_output=True, text=True)
        print(f"✅ Sent to IGV: {abs_path}")
    except subprocess.CalledProcessError as e:
        print(f"❌ Failed to send 'load' command to IGV. Is IGV running?")
        print(f"   Error: {e.stderr.strip()}")
        return # Stop if we can't load the file

    if region:
        goto_cmd = f'echo "goto {region}" | nc localhost {IGV_PORT}'
        try:
            subprocess.run(goto_cmd, shell=True, check=True)
            print(f"📍 Moved IGV view to: {region}")
        except subprocess.CalledProcessError:
            print(f"❌ Failed to send 'goto' command to IGV.")

def open_bed_files_in_igv(region, igv_app_path):
    if not os.path.exists(BED_DIR):
        print(f"❌ Directory not found: {BED_DIR}")
        return
    bed_files = [f for f in os.listdir(BED_DIR) if f.endswith((".bed", ".wig", ".bedGraph"))]
    if not bed_files:
        print("❌ No track files found in the results directory to load!")
        return
    
    launch_igv(igv_app_path)

    for i, bed_file in enumerate(bed_files):
        bed_path = os.path.join(BED_DIR, bed_file)
        if i == len(bed_files) - 1:
            send_to_igv(bed_path, region)
        else:
            send_to_igv(bed_path)

# --- Main Execution ---

def main():
    bed_tracks_dir, detailed_results_dir = create_directories()
    parser = argparse.ArgumentParser(description="Run genomic analysis and optionally view results in IGV.")
    
    parser.add_argument('genomic_coordinates', 
                        type=str, 
                        help="Required. Genomic coordinates in 'chr:start-end' format (e.g., 'chr1:632108-632403').")
    parser.add_argument('-refseq_functional', 
                        action='store_true', 
                        help="Extracts RefSeq functional element annotations (e.g., biological regions).")
    parser.add_argument('-eclips', 
                        action='store_true', 
                        help="Extracts eCLIP-seq peaks to identify RNA-Binding Protein (RBP) sites.")
    parser.add_argument('-SNP', 
                        action='store_true', 
                        help="Fetches all known SNPs for the region from the Ensembl REST API.")
    parser.add_argument('-miRNA', 
                        action='store_true', 
                        help="Extracts known miRNA annotations from a local GFF3 file.")
    parser.add_argument('-chem_mod', 
                        action='store_true', 
                        help="Fetches known RNA chemical modifications from the RMBase database.")
    parser.add_argument('-polyA', 
                        action='store_true', 
                        help="Extracts polyadenylation sites from the PolyASite 2.0 database.")
    parser.add_argument('-repeated-element', 
                        action='store_true', 
                        help="Fetches repetitive element annotations (e.g., SINEs, LINEs) from the Dfam API.")
    parser.add_argument('-chemical_prop', 
                        action='store_true', 
                        help="Extracts local chemical probing data (e.g., icSHAPE) from local WIG files.")
    parser.add_argument('-clinvar', 
                        action='store_true', 
                        help="Extracts clinical variants (SNPs and indels) from a local ClinVar VCF file.")
    parser.add_argument('-target_scan', 
                        action='store_true', 
                        help="Extracts predicted miRNA binding sites from local TargetScan BED files.")
    parser.add_argument('-phastCons', 
                        action='store_true', 
                        help="Fetches evolutionary conservation scores (phastCons 100-way vertebrates).")
    parser.add_argument('-igv', 
                        action='store_true', 
                        help="After analysis, automatically launch IGV and load all generated track files.")
    parser.add_argument('--igv-path', 
                        type=str, 
                        help="Required if using -igv. Provide the full path to your IGV application executable.")
   

    args = parser.parse_args()

    try:
        if args.genomic_coordinates.startswith("chr"):
            chrom, start_end = args.genomic_coordinates.split(":")
            start, end = map(int, start_end.split("-"))
        elif args.genomic_coordinates.startswith("ENST"):
           trans_cor_map=pd.read_csv("resources_data_sets/mart_export.csv",sep=',',header=0)
           trans_cor_map=trans_cor_map[trans_cor_map['Transcript stable ID']== args.genomic_coordinates] if "." not in args.genomic_coordinates else  trans_cor_map[trans_cor_map['Transcript stable ID version']== args.genomic_coordinates]
           chrom,start,end=f"chr{trans_cor_map['Chromosome/scaffold name'].iloc[0]}",f"{trans_cor_map['Transcript start (bp)'].iloc[0]}",f"{trans_cor_map['Transcript end (bp)'].iloc[0]}"
        else :
            raise ValueError("Genomic coordinates must be either coordinate(e.g chr:xxxx-yyyy) or valid ensembl transcript  id.")
            
    except ValueError:
        print("❌ Error: Invalid genomic coordinates format. Please use 'chr:start-end or valid ensembl transcript  id")
        sys.exit(1)
    
    region_string = f"{chrom}:{start}-{end}"
    print(f"🧬 Processing region: {region_string}")


    if args.chemical_prop: 
        chemical_prop_file = [os.path.join("resources_data_sets", f) for f in os.listdir("resources_data_sets") if f.endswith(".wig.gz")]
        if not chemical_prop_file:
            print("❌ No chemical property files found in 'resources_data_sets' directory.")
        else:
            chemical_prop(chrom, start, end, bed_tracks_dir, detailed_results_dir, chemical_prop_file)


    if args.refseq_functional:
        bed_files = ["resources_data_sets/FEbiolregions_AR110_GRCh38.p14.bed", "resources_data_sets/FEfeats_AR110_GRCh38.p14.bed", "resources_data_sets/FErecombpartners_AR110_GRCh38.p14.inter.bed", "resources_data_sets/FEregintxns_AR110_GRCh38.p14.inter.bed"]
        extract_functional_elements(chrom, start, end, bed_files, bed_tracks_dir, detailed_results_dir)
    if args.eclips:
        process_eclips_data(chrom, start, end, bed_tracks_dir, detailed_results_dir)
    if args.SNP:
        fetch_snps(chrom, start, end, bed_tracks_dir, detailed_results_dir)
    if args.miRNA:
        fetch_miRNA(chrom, start, end, bed_tracks_dir, detailed_results_dir)
    if args.chem_mod:
        fetch_chemmodi(chrom, start, end, bed_tracks_dir, detailed_results_dir)
    if args.polyA:
        polyA(chrom, start, end, bed_tracks_dir, detailed_results_dir)
    if args.repeated_element:
        repeated_element(chrom, start, end, bed_tracks_dir, detailed_results_dir)
    if args.clinvar:
        clinvar_file = "resources_data_sets/clinvar.vcf.gz"
        if not os.path.exists(clinvar_file):
            print(f"❌ Error: ClinVar file not found: {clinvar_file}")
        else:
            clinvar(chrom, start, end, bed_tracks_dir, detailed_results_dir, clinvar_file)
    if args.target_scan:
        target_scan_file = [os.path.join("resources_data_sets", f) for f in os.listdir("resources_data_sets") if f.startswith("Target")]
        if not target_scan_file:
            print("❌ No target scan files found in 'resources_data_sets' directory.")
        else:
            target_scan(chrom, start, end, bed_tracks_dir, detailed_results_dir, target_scan_file)
    if args.phastCons:
        phastCons_file ="resources_data_sets/hg38.phastCons100way.bedGraph.gz"
        phastCons(chrom, start, end, bed_tracks_dir, detailed_results_dir,phastCons_file)

    if args.igv:
        
        if not args.igv_path:
            print("❌ Error: To use the -igv feature, you must provide the path with --igv-path")
        elif not os.path.exists(args.igv_path):
            print(f"❌ Error: IGV application not found at the path you provided: {args.igv_path}")
        else:
            print("\n--- Starting IGV process ---")
            try:
                open_bed_files_in_igv(region_string, args.igv_path)
            except Exception as e:
                print(f"❌ Error during IGV process: {e}")

if __name__ == '__main__':
    main()