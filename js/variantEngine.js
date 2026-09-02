/**
 * Variant Engine
 * Parses HGVS cDNA notation, applies the mutation to the mRNA sequence,
 * translates WT and mutant CDS, and classifies the consequence.
 *
 * Coordinate convention:
 *   - All positions in this file are 1-based mRNA positions.
 *   - cdnaPos 1 = first nucleotide of the CDS = transcript position cdsStart.
 *   - Splicing notation: 123+1 means 1 nt into the intron after cdna pos 123.
 */

import { translateCDS, aa1to3, aaFullName } from './geneticCode.js';

// ─── HGVS Parser ─────────────────────────────────────────────────────────────

/**
 * Parse an HGVS cDNA string into a structured variant object.
 * Handles: substitution, deletion, duplication, insertion, delins, splicing_intronic.
 */
export function parseHGVS(hgvs) {
    let raw = hgvs.trim();
    if (raw.startsWith('c.')) raw = raw.slice(2);

    // Splicing intronic: 315+1G>A  774+1G>A  1585-1G>A
    const spl = raw.match(/^(\d+)([\+\-]\d+)([ACGT])>([ACGT])$/i);
    if (spl) return {
        type: 'splicing',
        anchorCdna: parseInt(spl[1]),
        offset: parseInt(spl[2]),
        ref: spl[3].toUpperCase(), alt: spl[4].toUpperCase(),
        hgvs: `c.${raw}`
    };

    // Substitution: 542C>T
    const sub = raw.match(/^([\-\*]?\d+)([ACGT])>([ACGT])$/i);
    if (sub) return {
        type: 'substitution',
        cdnaPos: parseCdna(sub[1]),
        ref: sub[2].toUpperCase(), alt: sub[3].toUpperCase(),
        hgvs: `c.${raw}`
    };

    // Deletion: 68del  68_69del  68_69delAG
    const del = raw.match(/^([\-\*]?\d+)(?:_([\-\*]?\d+))?del([ACGT]*)$/i);
    if (del) {
        const s = parseCdna(del[1]);
        const e = del[2] ? parseCdna(del[2]) : s;
        return { type: 'deletion', cdnaStart: s, cdnaEnd: e, hgvs: `c.${raw}` };
    }

    // Duplication: 5266dup  1274_1277dupTATC
    const dup = raw.match(/^([\-\*]?\d+)(?:_([\-\*]?\d+))?dup([ACGT]*)$/i);
    if (dup) {
        const s = parseCdna(dup[1]);
        const e = dup[2] ? parseCdna(dup[2]) : s;
        return { type: 'duplication', cdnaStart: s, cdnaEnd: e,
                 dupBases: dup[3] ? dup[3].toUpperCase() : null, hgvs: `c.${raw}` };
    }

    // Insertion: 4358_4359insGCTA
    const ins = raw.match(/^([\-\*]?\d+)_([\-\*]?\d+)ins([ACGT]+)$/i);
    if (ins) return {
        type: 'insertion',
        cdnaStart: parseCdna(ins[1]), cdnaEnd: parseCdna(ins[2]),
        inserted: ins[3].toUpperCase(), hgvs: `c.${raw}`
    };

    // Delins: 112_113delinsTG
    const di = raw.match(/^([\-\*]?\d+)(?:_([\-\*]?\d+))?delins([ACGT]+)$/i);
    if (di) {
        const s = parseCdna(di[1]);
        const e = di[2] ? parseCdna(di[2]) : s;
        return { type: 'delins', cdnaStart: s, cdnaEnd: e,
                 inserted: di[3].toUpperCase(), hgvs: `c.${raw}` };
    }

    throw new Error(`Variante "${hgvs}" no reconocida.\nEjemplos válidos:\n  c.542C>T\n  c.68_69delAG\n  c.5266dupC\n  c.315+1G>A\n  c.1521_1523delCTT`);
}

function parseCdna(s) {
    return parseInt(s.replace('*',''), 10);
}

/** Convert cDNA position to 0-based index in the mRNA sequence */
function cdnaToIdx(cdnaPos, cdsStart) {
    // cDNA pos 1 → index cdsStart-1
    return cdsStart - 1 + cdnaPos - 1;
}

// ─── Main simulation ──────────────────────────────────────────────────────────

export function simulateVariant(tx, hgvsStr) {
    const variant = typeof hgvsStr === 'string' ? parseHGVS(hgvsStr) : hgvsStr;
    const seq = tx.sequence;
    const cdsStart = tx.cdsStart; // 1-based

    // WT translation
    const wtCds = seq.slice(cdsStart - 1, tx.cdsEnd);
    const wtResult = translateCDS(wtCds);
    const wtProtein = wtResult.protein;

    let mutSeq = seq;
    let deltaNt = 0;
    let variantTxPos = null;   // 1-based mRNA position of the variant
    let isSplicing = false;

    if (variant.type === 'splicing') {
        // Pin position to the boundary nt of the exon (anchor cdna pos)
        variantTxPos = cdnaToIdx(variant.anchorCdna, cdsStart) + 1; // 1-based
        isSplicing = true;
        // No sequence change for splicing — consequence is predicted mechanistically

    } else if (variant.type === 'substitution') {
        const idx = cdnaToIdx(variant.cdnaPos, cdsStart);
        variantTxPos = idx + 1;
        mutSeq = seq.slice(0, idx) + variant.alt + seq.slice(idx + 1);
        deltaNt = 0;

    } else if (variant.type === 'deletion') {
        const startIdx = cdnaToIdx(variant.cdnaStart, cdsStart);
        const endIdx   = cdnaToIdx(variant.cdnaEnd,   cdsStart) + 1; // exclusive
        variantTxPos = startIdx + 1;
        deltaNt = -(endIdx - startIdx);
        mutSeq = seq.slice(0, startIdx) + seq.slice(endIdx);

    } else if (variant.type === 'duplication') {
        const startIdx = cdnaToIdx(variant.cdnaStart, cdsStart);
        const endIdx   = cdnaToIdx(variant.cdnaEnd,   cdsStart) + 1; // exclusive
        variantTxPos = endIdx; // insert after the duplicated segment
        const dupSeg = variant.dupBases || seq.slice(startIdx, endIdx);
        deltaNt = dupSeg.length;
        mutSeq = seq.slice(0, endIdx) + dupSeg + seq.slice(endIdx);

    } else if (variant.type === 'insertion') {
        const insertAfterIdx = cdnaToIdx(variant.cdnaStart, cdsStart) + 1;
        variantTxPos = insertAfterIdx + 1;
        deltaNt = variant.inserted.length;
        mutSeq = seq.slice(0, insertAfterIdx) + variant.inserted + seq.slice(insertAfterIdx);

    } else if (variant.type === 'delins') {
        const startIdx = cdnaToIdx(variant.cdnaStart, cdsStart);
        const endIdx   = cdnaToIdx(variant.cdnaEnd,   cdsStart) + 1;
        variantTxPos = startIdx + 1;
        deltaNt = variant.inserted.length - (endIdx - startIdx);
        mutSeq = seq.slice(0, startIdx) + variant.inserted + seq.slice(endIdx);
    }

    // Mutant translation
    const mutCdsStart = cdsStart; // CDS start doesn't move for coding variants
    const mutCds = mutSeq.slice(mutCdsStart - 1);
    const mutResult = translateCDS(mutCds);
    const mutProteinFull = mutResult.protein;

    // Find first stop in mutant
    const firstStop = mutProteinFull.indexOf('*');
    const mutProtein = firstStop >= 0 ? mutProteinFull.slice(0, firstStop) : mutProteinFull;

    // Classify
    const consequence = isSplicing
        ? classifySplicing(variant, tx)
        : classifyProtein(wtProtein, mutProteinFull, deltaNt, variant, cdsStart);

    // NMD prediction
    const nmd = predictNMD(tx, consequence.ptcCodonPos, consequence.isPrematureStop);

    return {
        variant,
        variantTxPos,
        isSplicing,
        wtProtein,
        mutProtein,
        mutProteinFull,
        wtCds,
        mutCds: mutSeq.slice(mutCdsStart - 1),
        mutSeq,
        deltaNt,
        readingFrameShift: ((Math.abs(deltaNt) % 3) + 3) % 3,
        consequence,
        nmd
    };
}

// ─── Consequence classifier ───────────────────────────────────────────────────

function classifySplicing(variant, tx) {
    const isDonor = variant.offset > 0;
    const offsetAbs = Math.abs(variant.offset);
    const isCanonical = isDonor
        ? (offsetAbs === 1 || offsetAbs === 2)
        : (offsetAbs === 1 || offsetAbs === 2);

    const mechanism = isCanonical
        ? (isDonor
            ? `Disrupción del sitio donante (+${offsetAbs}). Causa salto del exón o retención del intrón.`
            : `Disrupción del sitio aceptor (-${offsetAbs}). Causa salto del exón siguiente o retención del intrón.`)
        : `Variante intrónica no canónica (offset ${variant.offset}). Posible efecto en splicing críptico.`;

    return {
        category: 'Splicing (Canónico / Intrónico)',
        hgvsProtein: 'p.?',
        isFrameshift: false,
        isPrematureStop: false,
        ptcCodonPos: null,
        description: mechanism,
        badge: 'badge-warning'
    };
}

function classifyProtein(wtProt, mutProtFull, deltaNt, variant, cdsStart) {
    // Start loss
    if (!mutProtFull.startsWith('M')) {
        return {
            category: 'Pérdida de Inicio (Start-Loss)',
            hgvsProtein: 'p.Met1?',
            isFrameshift: false,
            isPrematureStop: false,
            ptcCodonPos: null,
            description: 'La variante elimina o altera el codón de inicio ATG (Met1). No se puede predecir la proteína resultante.',
            badge: 'badge-danger'
        };
    }

    // Find first difference
    const minLen = Math.min(wtProt.length, mutProtFull.length);
    let firstDiff = -1;
    for (let i = 0; i < minLen; i++) {
        if (wtProt[i] !== mutProtFull[i]) { firstDiff = i; break; }
    }
    if (firstDiff === -1) firstDiff = minLen;

    const aaPos = firstDiff + 1; // 1-based
    const wtAa  = wtProt[firstDiff]     || '';
    const mutAa = mutProtFull[firstDiff] || '';
    const isFrameshift = deltaNt !== 0 && (Math.abs(deltaNt) % 3 !== 0);

    // Frameshift
    if (isFrameshift) {
        const stopIdx = mutProtFull.indexOf('*', firstDiff);
        const fsLen   = stopIdx >= 0 ? (stopIdx - firstDiff + 1) : (mutProtFull.length - firstDiff);
        const ptcPos  = stopIdx >= 0 ? stopIdx + 1 : null;
        const hgvs    = `p.${aa1to3(wtAa)}${aaPos}${aa1to3(mutAa)}fs*${fsLen}`;
        return {
            category: 'Frameshift (Cambio de Marco)',
            hgvsProtein: hgvs,
            isFrameshift: true, isPrematureStop: true,
            ptcCodonPos: ptcPos,
            description: `Frameshift desde residuo ${aaPos} (${aa1to3(wtAa)}). Secuencia aberrante de ${fsLen} aa antes del PTC en posición ${ptcPos || 'no encontrado'}.`,
            badge: 'badge-danger'
        };
    }

    // Nonsense
    if (mutAa === '*') {
        return {
            category: 'Nonsense (Ganancia de Parada)',
            hgvsProtein: `p.${aa1to3(wtAa)}${aaPos}*`,
            isFrameshift: false, isPrematureStop: true,
            ptcCodonPos: aaPos,
            description: `Codón de parada prematuro en residuo ${aaPos} — p.${aa1to3(wtAa)}${aaPos}*. La proteína se trunca en ese punto.`,
            badge: 'badge-danger'
        };
    }

    // Silent
    if (firstDiff >= minLen && wtProt.length === mutProtFull.length) {
        return {
            category: 'Sinónima (Silent)',
            hgvsProtein: `p.${aa1to3(wtAa)}${aaPos}=`,
            isFrameshift: false, isPrematureStop: false,
            ptcCodonPos: null,
            description: 'Cambio sinónimo — la secuencia proteica no se altera.',
            badge: 'badge-neutral'
        };
    }

    // In-frame deletion
    if (variant.type === 'deletion' && Math.abs(deltaNt) % 3 === 0) {
        const nDeleted = Math.abs(deltaNt) / 3;
        const hgvs = nDeleted === 1
            ? `p.${aa1to3(wtAa)}${aaPos}del`
            : `p.${aa1to3(wtAa)}${aaPos}_${aa1to3(wtProt[aaPos + nDeleted - 2] || 'X')}${aaPos + nDeleted - 1}del`;
        return {
            category: 'In-Frame Deletion (Deleción en Marco)',
            hgvsProtein: hgvs,
            isFrameshift: false, isPrematureStop: false,
            ptcCodonPos: null,
            description: `Deleción de ${nDeleted} aminoácido(s) en marco. Pauta de lectura conservada.`,
            badge: 'badge-warning'
        };
    }

    // Missense
    if (variant.type === 'substitution' || (deltaNt === 0)) {
        return {
            category: 'Missense (Cambio de Sentido)',
            hgvsProtein: `p.${aa1to3(wtAa)}${aaPos}${aa1to3(mutAa)}`,
            isFrameshift: false, isPrematureStop: false,
            ptcCodonPos: null,
            description: `Sustitución en residuo ${aaPos}: ${aaFullName(wtAa)} (${aa1to3(wtAa)}) → ${aaFullName(mutAa)} (${aa1to3(mutAa)}).`,
            badge: 'badge-warning'
        };
    }

    return {
        category: 'In-Frame Indel',
        hgvsProtein: `p.${aa1to3(wtAa)}${aaPos}indel`,
        isFrameshift: false, isPrematureStop: false,
        ptcCodonPos: null,
        description: 'Indel en marco — pauta de lectura conservada.',
        badge: 'badge-warning'
    };
}

// ─── NMD predictor ────────────────────────────────────────────────────────────

export function predictNMD(tx, ptcCodonPos, isPrematureStop) {
    if (!isPrematureStop || !ptcCodonPos) {
        return {
            applies: false, isEscaping: false,
            status: 'No aplica', badge: 'badge-neutral',
            reason: 'La variante no genera un codón de parada prematuro (PTC).',
            distanceToLastJunction: null, ptcExon: null
        };
    }

    // Convert PTC codon position → mRNA position
    const ptcTxPos = tx.cdsStart + (ptcCodonPos - 1) * 3;

    const codingExons = tx.exons.filter(e => e.isCoding);
    if (codingExons.length <= 1) {
        return {
            applies: true, isEscaping: true,
            status: 'NMD Escape — gen monoexónico',
            badge: 'badge-warning',
            reason: 'Gen con un solo exón codificante. Sin EJCs río abajo: no hay NMD.',
            distanceToLastJunction: null, ptcExon: null
        };
    }

    // Find which exon contains the PTC
    let ptcExon = null;
    for (const ex of codingExons) {
        if (ptcTxPos >= ex.txStart && ptcTxPos <= ex.txEnd) { ptcExon = ex; break; }
    }

    const lastCodingExon = codingExons[codingExons.length - 1];
    const penultCodingExon = codingExons[codingExons.length - 2];

    // PTC in last exon → NMD escape
    if (ptcExon && ptcExon.exonNum === lastCodingExon.exonNum) {
        return {
            applies: true, isEscaping: true,
            status: 'NMD Escape — último exón',
            badge: 'badge-warning',
            reason: `El PTC está en el último exón codificante (Exón ${ptcExon.exonNum}). Sin EJC río abajo → proteína truncada estable.`,
            distanceToLastJunction: 0,
            ptcExon: ptcExon?.exonNum
        };
    }

    // Distance from PTC to last exon-exon junction
    const lastJunctionPos = penultCodingExon.txEnd;
    const dist = lastJunctionPos - ptcTxPos;

    if (dist < 55) {
        return {
            applies: true, isEscaping: true,
            status: `NMD Escape — regla 50-55 nt (${dist} nt)`,
            badge: 'badge-warning',
            reason: `El PTC está a solo ${dist} nt de la última unión exón-exón (<55 nt). El ribosoma desplaza el EJC → escape de NMD.`,
            distanceToLastJunction: dist,
            ptcExon: ptcExon?.exonNum
        };
    }

    return {
        applies: true, isEscaping: false,
        status: 'Trigger NMD',
        badge: 'badge-danger',
        reason: `El PTC está en el Exón ${ptcExon?.exonNum ?? '?'}, a ${dist} nt de la última unión exón-exón (>${55} nt). Los EJCs río abajo reclutan UPF1 → degradación del ARNm (Haploinsuficiencia).`,
        distanceToLastJunction: dist,
        ptcExon: ptcExon?.exonNum
    };
}
