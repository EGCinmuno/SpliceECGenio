# 🧬 SpliceECGenio — Genomic VariantStudio

**SpliceECGenio** es una plataforma bioinformática interactiva y pedagógica para la consulta de estructura génica (exones e intrones en coordenadas GRCh38), visualización a nivel macro y micro de secuencias nucleotídicas, y simulación en tiempo real del impacto de variantes de splicing y exónicas (desplazamiento de marco de lectura / frameshift, codones de parada prematuros PTC y activación de degradación mediada por NMD).

---

## 🚀 Características Principales

1. **Descarga y Mapeo Genómico Automático (Ensembl REST API):**
   - Resolución de variantes en múltiples formatos clínicos (`Chr:Pos Ref>Alt`, `Chr-Pos-Ref-Alt`, `NM_001089.3(ABCA3):c.875A>T`, etc.).
   - Mapeo automático de transcritos canónicos y descarga de estructura exón-intrón real en GRCh38.
   - Cálculo estricto de fases de lectura de entrada y salida (`Fase 0 ➔ 2`, `Fase 1 ➔ 0`, etc.) para cada exón codificante.

2. **Paso 3: Mapa Macro de Exones e Intrones:**
   - Vista balanceada y vista a escala genómica real.
   - Centrado automático interactivo en el locus de la variante.
   - Botones rápidos para saltar al inicio o final de cada exón con indicador de carga.

3. **Paso 4: Visor a Nivel de Bases y Límite de Splicing:**
   - Visualizador de nucleótidos (coloreados por base) con regla genómica, demarcación de límites de empalme y visualización de codones/aminoácidos en chevrons.
   - Soporte completo para variantes puntuales (SNVs), inserciones y deleciones (señalización de todas las bases WT afectadas y representación de alelos mutados con deleciones marcadas).
   - Conmutador entre **Secuencia Salvaje (WT)** y **Variante Mutada**.

4. **Paso 5: Clasificación y Métricas Biológicas:**
   - Deducción automática del tipo de impacto: Splicing Canónico (+1/+2, -1/-2), Splicing Flanqueante, Intrónica Profunda o Exónica.

5. **Paso 6: Simulador de Splicing y Consecuencias de Marco:**
   - **Salto de Exón Completo (*Exon Skipping*):** Comparación 1:1 de la estructura WT vs el transcripto maduro empalmado sin el exón omitido, evaluando si conserva el marco (*In-frame*) o si rompe la pauta (*Frameshift*).
   - **Activación de Sitios Crípticos (*Cryptic Sites*):** Visualización didáctica en dos fases:
     1. *Estructura Primaria (Pre-ARNm):* Secuencia genómica con el fragmento incorporado/eliminado y señalamiento de la variante.
     2. *Empalme Mutado:* ARNm maduro ensamblado y traducido a través de la nueva unión.
   - **Retención de Intrón Completo (*Intron Retention*):** Búsqueda automatizada en secuencia genómica real hasta el primer codón STOP prematuro (`TAA`, `TAG`, `TGA`), cortando la visualización con rigor biológico.
   - **Variantes Exónicas y Frameshift:** Recálculo respetando estrictamente la fase del exón del Paso 3, trazando los aminoácidos aberrantes downstream hasta el STOP prematuro.

---

## 💻 Instalación y Uso Local

No requiere configuración compleja ni librerías externas en el servidor. Funciona nativamente con tecnologías web modernas (HTML5, Vanilla CSS y ES Modules).

### 1. Clonar el repositorio:
```bash
git clone https://github.com/EGCinmuno/SpliceECGenio.git
cd SpliceECGenio
```

### 2. Levantar el servidor local:
Puedes usar Python:
```bash
python -m http.server 8080
```
O Node.js con cualquier servidor estático:
```bash
npx serve .
```

### 3. Abrir en el navegador:
Accede a `http://localhost:8080`.

---

## 🛠️ Tecnologías

- **Frontend:** HTML5 Semántico, Vanilla CSS (Variables CSS, Diseño responsivo, Modo oscuro/claro, Glassmorphism).
- **Lógica & Bioinformática:** JavaScript ES Modules (Vanilla JS nativo).
- **APIs:** [Ensembl REST API](https://rest.ensembl.org/) (Genoma de Referencia Humano GRCh38).

---

## 📄 Licencia

Este proyecto está bajo la Licencia MIT.
Desarrollado para **EGC Inmunogenética**.
