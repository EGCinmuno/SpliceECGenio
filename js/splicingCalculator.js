/**
 * Splicing Consequence and Exon Phase Calculator
 * Models exon skipping, reading frame disruptions, exon compatibility matrices,
 * and splice donor/acceptor mutations.
 */

import { translateCDS } from './geneticCode.js';
import { analyzeNMD } from './variantEngine.js';

/**
 * Calculates exon phases, codon boundary splits, and in-frame vs out-of-frame status for all exons
 * @param {Array<Object>} exons - List of exon objects
 * @param {number} cdsStart - 1-based transcript CDS start
 * @param {number} cdsEnd - 1-based transcript CDS end
 * @returns {Array<Object>} Enriched exons with phase metadata
 */
export function calculateExonPhases(exons, cdsStart, cdsEnd) {
    let cumulativeCodingLength = 0;
    
    return exons.map((exon, index) => {
        const isCoding = exon.end >= cdsStart && exon.start <= cdsEnd;
        let codingStart = null;
        let codingEnd = null;
        let codingLength = 0;
        let entryPhase = null; // 0, 1, 2
        let exitPhase = null;  // 0, 1, 2
        
        if (isCoding) {
            codingStart = Math.max(exon.start, cdsStart);
            codingEnd = Math.min(exon.end, cdsEnd);
            codingLength = (codingEnd - codingStart + 1);
            
            entryPhase = cumulativeCodingLength % 3;
            cumulativeCodingLength += codingLength;
            exitPhase = cumulativeCodingLength % 3;
        }
        
        const isInFrame = (codingLength % 3 === 0);
        const splitBasesAtEnd = exitPhase;
        
        return {
            ...exon,
            isCoding,
            codingStart,
            codingEnd,
            codingLength,
            entryPhase,
            exitPhase,
            isInFrame,
            phaseMod: codingLength % 3,
            phaseDescription: isCoding 
                ? `Phase ${entryPhase} → Phase ${exitPhase} (${codingLength} bp, ${isInFrame ? 'In-frame 3n' : `Shift +${codingLength % 3}`})`
                : 'Non-coding (UTR)'
        };
    });
}

/**
 * Simulates the effect of skipping one or more exons on the spliced transcript and protein
 * @param {Object} transcript - Transcript model
 * @param {Array<number>} skippedExonNumbers - Array of exon numbers to exclude (e.g. [12] or [45, 46, 47, 48, 49, 50])
 * @returns {Object} Outcome of exon skipping
 */
export function simulateExonSkipping(transcript, skippedExonNumbers) {
    const exons = transcript.exons;
    const skippedSet = new Set(skippedExonNumbers);
    
    // Retain only unskipped exons
    const remainingExons = exons.filter(e => !skippedSet.has(e.exonNum));
    const skippedExonsList = exons.filter(e => skippedSet.has(e.exonNum));
    
    let totalSkippedCodingBp = 0;
    skippedExonsList.forEach(e => {
        if (e.isCoding && e.codingLength) {
            totalSkippedCodingBp += e.codingLength;
        }
    });
    
    const isSkippingInFrame = (totalSkippedCodingBp % 3 === 0);
    const readingFrameShift = totalSkippedCodingBp % 3;
    
    // Reconstruct spliced cDNA sequence
    let splicedMutSeq = '';
    remainingExons.forEach(ex => {
        splicedMutSeq += transcript.sequence.substring(ex.start - 1, ex.end);
    });
    
    // Recalculate CDS coordinates in newly spliced transcript
    // Find where CDS starts in the new spliced sequence
    let newCdsStart = 1;
    let accumulated = 0;
    for (const ex of remainingExons) {
        if (transcript.cdsStart >= ex.start && transcript.cdsStart <= ex.end) {
            newCdsStart = accumulated + (transcript.cdsStart - ex.start + 1);
            break;
        }
        accumulated += (ex.end - ex.start + 1);
    }
    
    // Translate mutated spliced mRNA
    const mutCdsSeq = splicedMutSeq.substring(newCdsStart - 1);
    const mutTranslation = translateCDS(mutCdsSeq, 0);
    const fullMutProtein = mutTranslation.protein;
    
    const firstStopIdx = fullMutProtein.indexOf('*');
    const mutProtein = firstStopIdx >= 0 ? fullMutProtein.substring(0, firstStopIdx) : fullMutProtein;
    
    // NMD evaluation for skipped transcript
    let nmd = null;
    if (!isSkippingInFrame && firstStopIdx >= 0) {
        const ptcCodon = firstStopIdx + 1;
        // Mock virtual transcript with remaining exons for NMD calculation
        let virtualExons = [];
        let curr = 1;
        remainingExons.forEach(e => {
            const len = e.end - e.start + 1;
            virtualExons.push({
                exonNum: e.exonNum,
                start: curr,
                end: curr + len - 1,
                isCoding: e.isCoding
            });
            curr += len;
        });
        
        nmd = analyzeNMD({ exons: virtualExons }, ptcCodon, true, newCdsStart, -totalSkippedCodingBp);
    }
    
    return {
        skippedExonNumbers,
        skippedExonsList,
        totalSkippedCodingBp,
        isSkippingInFrame,
        readingFrameShift,
        consequenceSummary: isSkippingInFrame
            ? `In-frame deletion of ${totalSkippedCodingBp / 3} amino acids (${totalSkippedCodingBp} bp). Preserves downstream reading frame.`
            : `Out-of-frame deletion of ${totalSkippedCodingBp} bp (shift of +${readingFrameShift} nt). Causes downstream frameshift and premature termination.`,
        splicedMutSeq,
        mutProtein,
        fullMutProtein,
        nmd,
        ptcPosition: firstStopIdx >= 0 ? firstStopIdx + 1 : null
    };
}

/**
 * Calculates a 2D compatibility matrix between all coding exons
 * Showing whether joining Exon A directly to Exon B preserves reading frame.
 */
export function generateExonCompatibilityMatrix(exons) {
    const codingExons = exons.filter(e => e.isCoding);
    const matrix = [];
    
    for (let i = 0; i < codingExons.length; i++) {
        const row = [];
        const donorExon = codingExons[i];
        
        for (let j = 0; j < codingExons.length; j++) {
            const acceptorExon = codingExons[j];
            
            if (i >= j) {
                row.push({ valid: false, reason: 'Invalid Order' });
                continue;
            }
            
            // Calculate total coding bp skipped between donorExon and acceptorExon
            let skippedBp = 0;
            for (let k = i + 1; k < j; k++) {
                skippedBp += codingExons[k].codingLength || 0;
            }
            
            const inFrame = (skippedBp % 3 === 0);
            row.push({
                donor: donorExon.exonNum,
                acceptor: acceptorExon.exonNum,
                skippedBp,
                inFrame,
                shift: skippedBp % 3,
                label: inFrame ? 'In-Frame' : `Shift +${skippedBp % 3}`
            });
        }
        matrix.push({ donor: donorExon.exonNum, cells: row });
    }
    
    return { codingExons, matrix };
}
