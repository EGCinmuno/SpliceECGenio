"""
Script de Descarga de Estructura Genómica (Exones, Intrones y Secuencia)
Compatible con Ensembl REST API y NCBI Entrez.

Uso:
  python download.py
"""

import json
import re
import sys
import urllib.request
import urllib.error

# Configurar stdout para evitar errores de codificación en consolas Windows
if sys.platform == 'win32':
    try:
        sys.stdout.reconfigure(encoding='utf-8')
    except Exception:
        pass

# Opcional: si requests y biopython están instalados, los usa; de lo contrario usa urllib (stdlib)
try:
    import requests
    HAS_REQUESTS = True
except ImportError:
    HAS_REQUESTS = False

try:
    from Bio import Entrez, SeqIO
    Entrez.email = 'tu_correo@ejemplo.com'
    HAS_BIOPYTHON = True
except ImportError:
    HAS_BIOPYTHON = False

NCBI_CHROMOSOMES = {
    "1": "NC_000001.11", "2": "NC_000002.12", "3": "NC_000003.12", "4": "NC_000004.12",
    "5": "NC_000005.10", "6": "NC_000006.12", "7": "NC_000007.14", "8": "NC_000008.11",
    "9": "NC_000009.12", "10": "NC_000010.11", "11": "NC_000011.10", "12": "NC_000012.12",
    "13": "NC_000013.11", "14": "NC_000014.9", "15": "NC_000015.10", "16": "NC_000016.10",
    "17": "NC_000017.11", "18": "NC_000018.10", "19": "NC_000019.10", "20": "NC_000020.11",
    "21": "NC_000021.9", "22": "NC_000022.11", "X": "NC_000023.11", "Y": "NC_000024.10"
}

def http_get_json(url):
    """Realiza una petición GET y retorna el JSON parseado."""
    if HAS_REQUESTS:
        try:
            r = requests.get(url, headers={"Content-Type": "application/json"})
            if r.ok:
                return r.json()
        except Exception:
            pass
        return None
    else:
        req = urllib.request.Request(url, headers={"Content-Type": "application/json", "User-Agent": "NonTool/1.0"})
        try:
            with urllib.request.urlopen(req) as response:
                return json.loads(response.read().decode('utf-8'))
        except Exception:
            return None

def http_get_text(url, accept="text/plain"):
    """Realiza una petición GET y retorna el texto (FASTA)."""
    if HAS_REQUESTS:
        try:
            r = requests.get(url, headers={"Accept": accept})
            if r.ok:
                return r.text
        except Exception:
            pass
        return None
    else:
        req = urllib.request.Request(url, headers={"Accept": accept, "User-Agent": "NonTool/1.0"})
        try:
            with urllib.request.urlopen(req) as response:
                return response.read().decode('utf-8')
        except Exception:
            return None

def obtener_enst(nm_id):
    """Traduce un NM_ a ENST usando Ensembl."""
    nm_base = nm_id.split('.')[0]
    data = http_get_json(f"https://rest.ensembl.org/xrefs/symbol/homo_sapiens/{nm_base}?object_type=transcript")
    if data and isinstance(data, list):
        for xref in data:
            if xref.get('type') == 'transcript' and xref.get('id', '').startswith('ENST'):
                return xref['id']
    return None

def analizar_variante(variante_input):
    server = "https://rest.ensembl.org"
    partes = variante_input.split()
    
    match = re.match(r"(?:chr)?([0-9a-zA-Z]+):(\d+)", partes[0], re.IGNORECASE)
    if not match:
        print("[ERROR] Formato no valido. Usa ej: '1:234607399 A>G' o '5:169670615 A>G'")
        return
    chrom, pos = match.group(1), match.group(2)
    
    transcript_solicitado = partes[2] if len(partes) >= 3 else None
    print(f">> Iniciando analisis para chr{chrom}:{pos}...")

    transcripto_canonico = None
    gene_id, gene_name, strand = None, "Desconocido", 1

    # --- 1. MOTOR DE ESTRUCTURA ---
    if transcript_solicitado:
        print(f">> Analizando ID ingresado: {transcript_solicitado}")
        enst_id = transcript_solicitado if transcript_solicitado.startswith("ENST") else obtener_enst(transcript_solicitado)
        
        if enst_id:
            transcripto_canonico = http_get_json(server + f"/lookup/id/{enst_id}?expand=1")
            if transcripto_canonico:
                gene_id = transcripto_canonico.get('Parent')
                r_gen = http_get_json(server + f"/lookup/id/{gene_id}")
                if r_gen:
                    gene_name = r_gen.get('display_name', gene_id)
                    strand = r_gen.get('strand', 1)
        else:
            print("[WARN] No se pudo traducir a un ID genomico.")
            return

    else:
        print(">> Buscando gen solapante por defecto...")
        ext_overlap = f"/overlap/region/human/{chrom}:{pos}-{pos}?feature=gene"
        genes = http_get_json(server + ext_overlap)
        
        if genes and len(genes) > 0:
            gen = genes[0]
            gene_id = gen['id']
            gen_data = http_get_json(server + f"/lookup/id/{gene_id}?expand=1")
            if gen_data:
                gene_name = gen_data.get('display_name', gene_id)
                strand = gen_data.get('strand', 1)
                transcriptos = gen_data.get('Transcript', [])
                transcripto_canonico = next((t for t in transcriptos if t.get('is_canonical') == 1), transcriptos[0] if transcriptos else None)
        else:
            print("[ERROR] No se encontro ningun gen en esa coordenada exacta.")
            return

    # --- 2. IMPRIMIR ESTRUCTURA Y COORDENADAS ---
    if transcripto_canonico:
        print(f"\n[OK] Gen: {gene_name} | Transcripto Genomico: {transcripto_canonico['id']}")
        exones = sorted(transcripto_canonico.get('Exon', []), key=lambda x: x['start'])
        print(f"[*] Exones: {len(exones)} | Intrones: {max(0, len(exones) - 1)}")
        
        for i, exon in enumerate(exones):
            num_exon = i + 1 if strand == 1 else len(exones) - i
            print(f"  -> Exon {num_exon}: {exon['start']} - {exon['end']}")
            if i < len(exones) - 1:
                print(f"     -- Intron: {exon['end']+1} - {exones[i+1]['start']-1}")
        
        start = transcripto_canonico['start']
        end = transcripto_canonico['end']
        
        # --- 3. MOTOR DUAL DE DESCARGA ---
        print(f"\n>> Descargando secuencia GENOMICA (con intrones) de {gene_name}...")
        nombre_archivo = f"{gene_name}_genomic_chr{chrom}_{start}_{end}.fasta"
        descarga_exitosa = False

        # INTENTO 1: ENSEMBL
        print("   [1/2] Conectando con Ensembl API...")
        ext_region = f"/sequence/region/human/{chrom}:{start}..{end}:{strand}"
        seq_text = http_get_text(server + ext_region, accept="text/x-fasta")
        
        if seq_text:
            with open(nombre_archivo, "w", encoding="utf-8") as f:
                f.write(seq_text)
            print(f"   [OK] Guardado en: {nombre_archivo}")
            descarga_exitosa = True
        else:
            print("   [WARN] Ensembl no retorno secuencia.")

        # INTENTO 2: NCBI
        if not descarga_exitosa and HAS_BIOPYTHON:
            print("   [2/2] Activando NCBI (Entrez)...")
            ncbi_chr_id = NCBI_CHROMOSOMES.get(str(chrom))
            if ncbi_chr_id:
                try:
                    ncbi_strand = 1 if strand == 1 else 2
                    handle = Entrez.efetch(
                        db="nuccore", id=ncbi_chr_id, rettype="fasta", retmode="text",
                        seq_start=start, seq_stop=end, strand=ncbi_strand
                    )
                    record = SeqIO.read(handle, "fasta")
                    handle.close()
                    record.id = f"{gene_name}_Genomic_Region"
                    record.description = f"chr{chrom}:{start}-{end}"
                    SeqIO.write(record, nombre_archivo, "fasta")
                    print(f"   [OK] Guardado desde NCBI en: {nombre_archivo}")
                except Exception as e:
                    print(f"   [ERROR] NCBI fallo: {e}")

# Ejecución de prueba
if __name__ == "__main__":
    print("--- Probando Variante 1: SMN1 (Chr5) ---")
    analizar_variante("5:169670615 A>G NM_000344.4")
    
    print("\n" + "="*50 + "\n")
    print("--- Probando Variante 2: TARBP1 (Chr1) ---")
    analizar_variante("1:234607399 A>G")
