/**
 * Genetic Code — Standard codon table, amino acid properties, and translation utilities.
 */

export const CODON_TABLE = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L',
    'CTT':'L','CTC':'L','CTA':'L','CTG':'L',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M',
    'GTT':'V','GTC':'V','GTA':'V','GTG':'V',
    'TCT':'S','TCC':'S','TCA':'S','TCG':'S',
    'CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'ACT':'T','ACC':'T','ACA':'T','ACG':'T',
    'GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*',
    'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K',
    'GAT':'D','GAC':'D','GAA':'E','GAG':'E',
    'TGT':'C','TGC':'C','TGA':'*','TGG':'W',
    'CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'AGT':'S','AGC':'S','AGA':'R','AGG':'R',
    'GGT':'G','GGC':'G','GGA':'G','GGG':'G'
};

export const STOP_CODONS = ['TAA','TAG','TGA'];

// All codons that are NOT stop codons
export const SENSE_CODONS = Object.keys(CODON_TABLE).filter(c => CODON_TABLE[c] !== '*');

export const AA_PROPS = {
    'A':{ name:'Alanine',   code3:'Ala', color:'#10b981' },
    'R':{ name:'Arginine',  code3:'Arg', color:'#3b82f6' },
    'N':{ name:'Asparagine',code3:'Asn', color:'#8b5cf6' },
    'D':{ name:'Aspartate', code3:'Asp', color:'#ef4444' },
    'C':{ name:'Cysteine',  code3:'Cys', color:'#f59e0b' },
    'Q':{ name:'Glutamine', code3:'Gln', color:'#8b5cf6' },
    'E':{ name:'Glutamate', code3:'Glu', color:'#ef4444' },
    'G':{ name:'Glycine',   code3:'Gly', color:'#6b7280' },
    'H':{ name:'Histidine', code3:'His', color:'#3b82f6' },
    'I':{ name:'Isoleucine',code3:'Ile', color:'#10b981' },
    'L':{ name:'Leucine',   code3:'Leu', color:'#10b981' },
    'K':{ name:'Lysine',    code3:'Lys', color:'#3b82f6' },
    'M':{ name:'Methionine',code3:'Met', color:'#06b6d4' },
    'F':{ name:'Phe',       code3:'Phe', color:'#ec4899' },
    'P':{ name:'Proline',   code3:'Pro', color:'#f97316' },
    'S':{ name:'Serine',    code3:'Ser', color:'#a855f7' },
    'T':{ name:'Threonine', code3:'Thr', color:'#a855f7' },
    'W':{ name:'Tryptophan',code3:'Trp', color:'#ec4899' },
    'Y':{ name:'Tyrosine',  code3:'Tyr', color:'#ec4899' },
    'V':{ name:'Valine',    code3:'Val', color:'#10b981' },
    '*':{ name:'Stop',      code3:'Ter', color:'#dc2626' },
    'X':{ name:'Unknown',   code3:'Xaa', color:'#9ca3af' }
};

export function aa1to3(aa) {
    if (!aa) return 'Xaa';
    if (aa === '*') return 'Ter';
    return AA_PROPS[aa]?.code3 || aa;
}

export function aaFullName(aa) {
    return AA_PROPS[aa]?.name || 'Unknown';
}

/** Translate a raw DNA sequence (CDS) starting at offset=0 */
export function translateCDS(dna) {
    const s = (dna || '').toUpperCase().replace(/[^ACGT]/g, 'N');
    let protein = '';
    const codons = [];
    for (let i = 0; i + 2 < s.length; i += 3) {
        const codon = s.substring(i, i + 3);
        const aa = CODON_TABLE[codon] || (codon.includes('N') ? 'X' : '?');
        codons.push({ codon, aa, pos: codons.length + 1 });
        protein += aa;
        // Don't stop — include stop codons in output for analysis
    }
    return { protein, codons };
}
