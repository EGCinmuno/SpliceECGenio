/**
 * NCBI / RefSeq Data Service
 *
 * Provides pre-curated transcript models for 7 clinical genes.
 * Each transcript has:
 *   - id, gene, chromosome, strand
 *   - cdsStart / cdsEnd  (1-based positions in the mRNA sequence)
 *   - exons: array of { exonNum, txStart, txEnd } — positions in the mRNA
 *   - sampleVariants: array of { hgvs, type, desc }
 *
 * All coordinates are 1-based, in mRNA (transcript) space.
 * We generate a synthetic but deterministic mRNA sequence so that
 * translation always returns real amino acids.
 */

import { SENSE_CODONS, translateCDS } from './geneticCode.js';

// ─── Transcript Database ──────────────────────────────────────────────────────

export const TRANSCRIPTS = {
    BRCA1: {
        id: 'NM_007294.4', gene: 'BRCA1',
        chromosome: '17', strand: '-',
        genomicStart: 43044290,
        description: 'Breast/ovarian cancer susceptibility — RING, BRCT domains',
        // mRNA total length ~7088 nt; CDS: nt 233-5824
        cdsStart: 233, cdsEnd: 5824,
        // Exons in mRNA coordinates (1-based)
        exons: [
            { exonNum:1,  txStart:1,    txEnd:100  },
            { exonNum:2,  txStart:101,  txEnd:299  },
            { exonNum:3,  txStart:300,  txEnd:353  },
            { exonNum:4,  txStart:354,  txEnd:431  },
            { exonNum:5,  txStart:432,  txEnd:520  },
            { exonNum:6,  txStart:521,  txEnd:660  },
            { exonNum:7,  txStart:661,  txEnd:766  },
            { exonNum:8,  txStart:767,  txEnd:872  },
            { exonNum:9,  txStart:873,  txEnd:918  },
            { exonNum:10, txStart:919,  txEnd:4344 },
            { exonNum:11, txStart:4345, txEnd:4433 },
            { exonNum:12, txStart:4434, txEnd:4560 },
            { exonNum:13, txStart:4561, txEnd:4687 },
            { exonNum:14, txStart:4688, txEnd:4878 },
            { exonNum:15, txStart:4879, txEnd:5189 },
            { exonNum:16, txStart:5190, txEnd:5277 },
            { exonNum:17, txStart:5278, txEnd:5355 },
            { exonNum:18, txStart:5356, txEnd:5396 },
            { exonNum:19, txStart:5397, txEnd:5480 },
            { exonNum:20, txStart:5481, txEnd:5564 },
            { exonNum:21, txStart:5565, txEnd:5638 },
            { exonNum:22, txStart:5639, txEnd:5712 },
            { exonNum:23, txStart:5713, txEnd:5773 },
            { exonNum:24, txStart:5774, txEnd:7088 }
        ],
        sampleVariants: [
            { hgvs:'c.68_69delAG', type:'frameshift', desc:'Patogénica — p.Glu23Valfs*17' },
            { hgvs:'c.5266dupC',   type:'frameshift', desc:'Exón 20 — p.Gln1756Profs*74' },
            { hgvs:'c.181T>G',     type:'missense',   desc:'RING domain — p.Cys61Gly' },
            { hgvs:'c.5123C>A',    type:'nonsense',   desc:'BRCT domain — p.Ser1708*' }
        ]
    },

    CFTR: {
        id: 'NM_000492.4', gene: 'CFTR',
        chromosome: '7', strand: '+',
        genomicStart: 117465780,
        description: 'Cystic fibrosis transmembrane conductance regulator',
        cdsStart: 133, cdsEnd: 4575,
        exons: [
            { exonNum:1,  txStart:1,    txEnd:182  },
            { exonNum:2,  txStart:183,  txEnd:289  },
            { exonNum:3,  txStart:290,  txEnd:396  },
            { exonNum:4,  txStart:397,  txEnd:609  },
            { exonNum:5,  txStart:610,  txEnd:700  },
            { exonNum:6,  txStart:701,  txEnd:828  },
            { exonNum:7,  txStart:829,  txEnd:947  },
            { exonNum:8,  txStart:948,  txEnd:1047 },
            { exonNum:9,  txStart:1048, txEnd:1230 },
            { exonNum:10, txStart:1231, txEnd:1422 },
            { exonNum:11, txStart:1423, txEnd:1614 },
            { exonNum:12, txStart:1615, txEnd:1720 },
            { exonNum:13, txStart:1721, txEnd:1827 },
            { exonNum:14, txStart:1828, txEnd:2424 },
            { exonNum:15, txStart:2425, txEnd:2541 },
            { exonNum:16, txStart:2542, txEnd:2638 },
            { exonNum:17, txStart:2639, txEnd:2795 },
            { exonNum:18, txStart:2796, txEnd:2977 },
            { exonNum:19, txStart:2978, txEnd:3121 },
            { exonNum:20, txStart:3122, txEnd:3267 },
            { exonNum:21, txStart:3268, txEnd:3377 },
            { exonNum:22, txStart:3378, txEnd:3512 },
            { exonNum:23, txStart:3513, txEnd:3672 },
            { exonNum:24, txStart:3673, txEnd:3820 },
            { exonNum:25, txStart:3821, txEnd:3936 },
            { exonNum:26, txStart:3937, txEnd:4056 },
            { exonNum:27, txStart:4057, txEnd:6128 }
        ],
        sampleVariants: [
            { hgvs:'c.1521_1523delCTT', type:'in-frame',  desc:'p.Phe508del — DeltaF508, más frecuente' },
            { hgvs:'c.1624G>T',          type:'nonsense',   desc:'p.Gly542* — G542X' },
            { hgvs:'c.1585-1G>A',        type:'splicing',   desc:'Splicing aceptor canónico intrón 11' }
        ]
    },

    TP53: {
        id: 'NM_000546.6', gene: 'TP53',
        chromosome: '17', strand: '-',
        genomicStart: 7668400,
        description: 'Tumor suppressor p53 — guardián del genoma',
        cdsStart: 187, cdsEnd: 1368,
        exons: [
            { exonNum:1,  txStart:1,    txEnd:186  },
            { exonNum:2,  txStart:187,  txEnd:293  },
            { exonNum:3,  txStart:294,  txEnd:315  },
            { exonNum:4,  txStart:316,  txEnd:594  },
            { exonNum:5,  txStart:595,  txEnd:778  },
            { exonNum:6,  txStart:779,  txEnd:891  },
            { exonNum:7,  txStart:892,  txEnd:1001 },
            { exonNum:8,  txStart:1002, txEnd:1138 },
            { exonNum:9,  txStart:1139, txEnd:1212 },
            { exonNum:10, txStart:1213, txEnd:1319 },
            { exonNum:11, txStart:1320, txEnd:2591 }
        ],
        sampleVariants: [
            { hgvs:'c.637C>T',    type:'nonsense',  desc:'p.Arg213* — hotspot nonsense' },
            { hgvs:'c.743G>A',    type:'missense',   desc:'p.Arg248Gln — dominant negative' },
            { hgvs:'c.375+1G>A',  type:'splicing',   desc:'Splicing donante canónico exón 4' }
        ]
    },

    HBB: {
        id: 'NM_000518.5', gene: 'HBB',
        chromosome: '11', strand: '-',
        genomicStart: 5225460,
        description: 'Beta-globin — Beta-Thalassemia, Sickle Cell Disease',
        cdsStart: 51, cdsEnd: 494,
        exons: [
            { exonNum:1, txStart:1,   txEnd:142 },
            { exonNum:2, txStart:143, txEnd:365 },
            { exonNum:3, txStart:366, txEnd:626 }
        ],
        sampleVariants: [
            { hgvs:'c.315+1G>A', type:'splicing',  desc:'Beta-0 Thalassemia — donor canónico' },
            { hgvs:'c.118C>T',   type:'nonsense',   desc:'p.Gln40* — Beta-0 Codon 39' },
            { hgvs:'c.20A>T',    type:'missense',   desc:'p.Glu6Val — HbS Sickle Cell' }
        ]
    },

    SMN1: {
        id: 'NM_000344.4', gene: 'SMN1',
        chromosome: '5', strand: '+',
        genomicStart: 169670560,
        description: 'Survival motor neuron 1 — Spinal Muscular Atrophy (SMA)',
        cdsStart: 85, cdsEnd: 969,
        exons: [
            { exonNum:1, txStart:1,   txEnd:140  },
            { exonNum:2, txStart:141, txEnd:270  },
            { exonNum:3, txStart:271, txEnd:380  },
            { exonNum:4, txStart:381, txEnd:490  },
            { exonNum:5, txStart:491, txEnd:600  },
            { exonNum:6, txStart:601, txEnd:720  },
            { exonNum:7, txStart:721, txEnd:774  },
            { exonNum:8, txStart:775, txEnd:1450 }
        ],
        sampleVariants: [
            { hgvs:'c.774+1G>A', type:'splicing',  desc:'Exón 7 splicing donante — Chr5:169670615' },
            { hgvs:'c.840C>T',   type:'synonymous', desc:'SMN1 vs SMN2 differentiator exón 7' }
        ]
    },

    DMD: {
        id: 'NM_004006.3', gene: 'DMD',
        chromosome: 'X', strand: '+',
        genomicStart: 31119220,
        description: 'Dystrophin (79 exons) — Duchenne / Becker Muscular Dystrophy',
        cdsStart: 210, cdsEnd: 11264,
        exons: [
            { exonNum:1,  txStart:1,    txEnd:261  },
            { exonNum:2,  txStart:262,  txEnd:437  },
            { exonNum:3,  txStart:438,  txEnd:585  },
            { exonNum:4,  txStart:586,  txEnd:735  },
            { exonNum:5,  txStart:736,  txEnd:921  },
            { exonNum:6,  txStart:922,  txEnd:1023 },
            { exonNum:7,  txStart:1024, txEnd:1132 },
            { exonNum:8,  txStart:1133, txEnd:1365 },
            { exonNum:9,  txStart:1366, txEnd:1483 },
            { exonNum:10, txStart:1484, txEnd:1695 },
            { exonNum:11, txStart:1696, txEnd:1800 },
            { exonNum:12, txStart:1801, txEnd:1900 }
        ],
        sampleVariants: [
            { hgvs:'c.1133_1365del', type:'frameshift', desc:'Deleción exón 8 — Duchenne (out-of-frame)' },
            { hgvs:'c.438_1132del',  type:'in-frame',   desc:'Deleción exones 3-7 — Becker (in-frame)' }
        ]
    },

    HEXA: {
        id: 'NM_000520.6', gene: 'HEXA',
        chromosome: '15', strand: '+',
        genomicStart: 72342500,
        description: 'Hexosaminidase alpha — Tay-Sachs disease',
        cdsStart: 85, cdsEnd: 1674,
        exons: [
            { exonNum:1,  txStart:1,    txEnd:120  },
            { exonNum:2,  txStart:121,  txEnd:210  },
            { exonNum:3,  txStart:211,  txEnd:320  },
            { exonNum:4,  txStart:321,  txEnd:415  },
            { exonNum:5,  txStart:416,  txEnd:508  },
            { exonNum:6,  txStart:509,  txEnd:602  },
            { exonNum:7,  txStart:603,  txEnd:776  },
            { exonNum:8,  txStart:777,  txEnd:890  },
            { exonNum:9,  txStart:891,  txEnd:1010 },
            { exonNum:10, txStart:1011, txEnd:1140 },
            { exonNum:11, txStart:1141, txEnd:1300 },
            { exonNum:12, txStart:1301, txEnd:1420 },
            { exonNum:13, txStart:1421, txEnd:1535 },
            { exonNum:14, txStart:1536, txEnd:2050 }
        ],
        sampleVariants: [
            { hgvs:'c.1274_1277dupTATC', type:'frameshift', desc:'p.Tyr427Ilefs*5 — Tay-Sachs insertion' }
        ]
    }
};

// ─── Synthetic sequence generator ────────────────────────────────────────────
// Creates a deterministic, biologically valid mRNA sequence for a transcript.
// 5'UTR = random ACGT; CDS starts with ATG, uses only sense codons, ends with TAA; 3'UTR = random.

function seededRng(seed) {
    let s = Math.abs(seed) || 1;
    return () => {
        s = (s * 1664525 + 1013904223) & 0xFFFFFFFF;
        return (s >>> 0) / 0xFFFFFFFF;
    };
}

export function buildSequence(tx) {
    const bases = ['A','C','G','T'];
    // Seed from the transcript ID string
    let seedVal = 0;
    for (let i = 0; i < tx.id.length; i++) seedVal += tx.id.charCodeAt(i) * (i + 1);
    const rng = seededRng(seedVal);

    const totalLen = tx.exons[tx.exons.length - 1].txEnd;
    const arr = new Array(totalLen);

    // Fill everything with random bases first
    for (let i = 0; i < totalLen; i++) arr[i] = bases[Math.floor(rng() * 4)];

    // Overwrite CDS with valid ATG + sense codons + TAA stop
    const cdsLen = tx.cdsEnd - tx.cdsStart + 1;  // includes stop codon
    arr[tx.cdsStart - 1] = 'A';
    arr[tx.cdsStart]     = 'T';
    arr[tx.cdsStart + 1] = 'G';

    const numCodons = Math.floor((cdsLen - 6) / 3); // -3 for ATG, -3 for TAA
    for (let c = 0; c < numCodons; c++) {
        const codon = SENSE_CODONS[Math.floor(rng() * SENSE_CODONS.length)];
        const offset = tx.cdsStart + 3 + c * 3; // 0-based
        arr[offset]     = codon[0];
        arr[offset + 1] = codon[1];
        arr[offset + 2] = codon[2];
    }
    // Stop codon TAA
    const stopOffset = tx.cdsEnd - 3; // 0-based
    arr[stopOffset]     = 'T';
    arr[stopOffset + 1] = 'A';
    arr[stopOffset + 2] = 'A';

    return arr.join('');
}

// ─── Exon phase annotation ────────────────────────────────────────────────────
// Adds isCoding, codingLen, entryPhase, exitPhase to each exon

export function annotateExons(tx) {
    const { cdsStart, cdsEnd, exons } = tx;
    let cumCoding = 0;
    return exons.map(ex => {
        const overlapStart = Math.max(ex.txStart, cdsStart);
        const overlapEnd   = Math.min(ex.txEnd,   cdsEnd);
        const isCoding = overlapEnd >= overlapStart;

        if (!isCoding) {
            return { ...ex, isCoding: false, codingStart: null, codingEnd: null,
                     codingLen: 0, entryPhase: null, exitPhase: null };
        }
        const codingLen  = overlapEnd - overlapStart + 1;
        const entryPhase = cumCoding % 3;
        cumCoding += codingLen;
        const exitPhase  = cumCoding % 3;

        return { ...ex, isCoding: true,
                 codingStart: overlapStart, codingEnd: overlapEnd,
                 codingLen, entryPhase, exitPhase };
    });
}

// ─── Public API ──────────────────────────────────────────────────────────────

/**
 * Returns a transcript object ready for simulation.
 * Accepts:
 *   - a gene key like 'BRCA1'
 *   - an NM accession like 'NM_007294.4'
 */
export function getTranscript(query) {
    const q = query.trim().toUpperCase();

    // Match by gene key
    if (TRANSCRIPTS[q]) {
        const tx = { ...TRANSCRIPTS[q] };
        tx.sequence = buildSequence(tx);
        tx.exons    = annotateExons(tx);
        tx.proteinLength = Math.floor((tx.cdsEnd - tx.cdsStart + 1 - 3) / 3); // excl stop
        return tx;
    }

    // Match by NM accession (partial match allowed, e.g. NM_007294)
    for (const key of Object.keys(TRANSCRIPTS)) {
        const t = TRANSCRIPTS[key];
        if (t.id.toUpperCase().startsWith(q) || q.startsWith(t.id.toUpperCase().split('.')[0])) {
            const tx = { ...t };
            tx.sequence = buildSequence(tx);
            tx.exons    = annotateExons(tx);
            tx.proteinLength = Math.floor((tx.cdsEnd - tx.cdsStart + 1 - 3) / 3);
            return tx;
        }
    }

    throw new Error(`Transcripto "${query}" no encontrado. Prueba: BRCA1, CFTR, TP53, HBB, SMN1, DMD, HEXA o su NM_*.`);
}
