# OBJETIVO: atlas agéntico multiómico de colitis ulcerosa

Construir y dejar funcionando, de principio a fin, el mismo sistema que el atlas de cáncer renal
([agentic-kidney-atlas](https://github.com/rsolerortuno/agentic-kidney-atlas)), aplicado a la **colitis ulcerosa (CU)**:
datos humanos públicos + literatura + APIs externas en un grafo de conocimiento con evidencia trazable, servidores MCP,
un orquestador y especialistas de Claude Code, verificación determinista de cada afirmación y un informe web.

Pregunta: **"¿Qué dianas terapéuticas, vías inmunes y biomarcadores de la colitis ulcerosa tienen la evidencia
convergente más fuerte?"**

Identificadores: MONDO:0005101, EFO_0000729 (Open Targets), DOID:8577. Comprobar en M1 que Open Targets usa EFO_0000729.

## Restricciones (heredadas del atlas renal)
- Windows 11; Python 3.12 con uv; git, **un commit por milestone**. Repositorio nuevo; el atlas renal no se toca.
- LLM: solo Claude Code con la suscripción (sin API de Claude, SDK ni ANTHROPIC_API_KEY). Orquestador: Opus 5.5;
  especialistas: `model: sonnet` (Sonnet 5.5). No lanzar dos `claude -p` a la vez.
- Todo cálculo numérico (expresión diferencial, metaanálisis, convergencia, verificación) es Python determinista.
- Ninguna afirmación sin evidence_id; los números se enlazan a su registro con `numeric_annotation`.
- Solo datos públicos. Nada de datos internos de ninguna empresa. Ningún archivo con abstracts o matrices crudas en git.
- Estado en PROGRESS.md: milestone actual, comprobación superada con cifras, bloqueos.
- Un milestone se cierra solo cuando su comprobación (✔) pasa. Si falla, se arregla o se documenta por qué no.

## Lecciones del atlas renal que se aplican desde el principio
- **Símbolos de genes:** mapear siempre por HGNC con símbolos previos y alias (en el renal, COXFA4L2/NDUFA4L2 perdió su
  replicación por usar el símbolo antiguo).
- **Salida estructurada:** toda respuesta de modelo que lea un programa usa un esquema JSON validado.
- **Paso de cobertura:** los especialistas marcan hallazgos clave (`"key": true`), y `coverage.py` comprueba que llegan
  al informe.
- **Huella de evidencia** (`basis.json`): la verificación final falla si la evidencia citada cambió.
- **Ejecuciones aisladas:** cada `/informe` en su carpeta (`RUN_DIR`), y réplicas independientes.
- **Rutas largas en Windows:** proyecto en `~\Documents\uc-graphrag`.
- **Composición celular antes que "diana":** en CU, mucho de lo que sube en bulk es infiltración inmune (equivalente a la
  pureza tumoral en el renal). Ningún gen se llama diana sin pasar esa comprobación.

---

## Milestones, en orden de ejecución

### Bloque A: motor y datos base

**M0 Entorno y motor genérico.** Repo nuevo con uv. Copiar del atlas renal el motor que no depende de la enfermedad:
`GraphStore` (SQLite + NetworkX), identificadores de evidencia, verificador (existencia, gen, números, dirección,
huella), `coverage.py`, convergencia configurable, andamiaje de servidores MCP, `render`, manifiesto de reproducibilidad
y benchmark de agentes. La enfermedad, los umbrales y las fuentes van en `config.yaml`.
- ✔ `uv run pytest` en verde con datos de juguete; `grep` no encuentra "ccRCC", "renal" ni "kidney" en el motor.

**M1 Verificación y descarga de datos bulk.** Comprobar en GEO, antes de descargar, muestras, plataforma y metadatos
(inflamación, sanos, respuesta, puntuación Mayo) de cada cohorte, y anotarlo en una tabla en `data/README.md`:

| Uso | Conjunto | Lo que se espera (comprobar) |
|---|---|---|
| Descubrimiento | GSE193677 (Mount Sinai, biopsias de adultos con EII y controles, RNA-seq) | Muestras de CU inflamadas y no inflamadas, y controles |
| Replicación | GSE206285 (UNIFI, ustekinumab) | Biopsias de sigma basales de unos 550 pacientes y 18 sanos, micromatriz |
| Replicación | GSE92415 (PURSUIT, golimumab) | 87 basales, 21 sanos, semana 6 |
| Replicación | GSE16879 (infliximab) y GSE73661 (infliximab y vedolizumab) | Basales con respondedores y no respondedores |
| Referencias | HGNC completo, anotación de cada plataforma de micromatriz | |

Manifiesto con URL, fecha, sha256 y dimensiones.
- ✔ Recuentos de cada cohorte iguales a los anotados (o la diferencia explicada); cada muestra con estado (CU
  inflamada / CU no inflamada / sano) y, si existe, respuesta y puntuación endoscópica.

**M2 Expresión diferencial y metaanálisis.** Por cohorte y por separado (las plataformas no se mezclan en crudo):
- CU inflamada frente a sano (modelo lineal sobre expresión log, con edad y sexo si están);
- inflamada frente a no inflamada del mismo paciente, donde exista (análisis emparejado, como tumor/normal en CPTAC);
- metaanálisis por gen (efectos aleatorios), y fracción de cambios que se replican en cada cohorte independiente.

Salida en `data/processed/*.parquet`.
- ✔ Controles positivos significativos (padj < 0,05) y en la misma dirección en al menos 2 cohortes:
  - **suben** S100A8, LCN2, REG4, DUOX2, CXCL8, MMP3, IL1B y OSM;
  - **bajan** AQP8, SLC26A2, HMGCS2 y ABCG2.
- ✔ `tests/test_expression.py` en verde.

**M3 Literatura.** PubTator3 para CU y los 50 genes principales de M2, igual que en el renal; extracción de relaciones
con esquema cerrado (máximo 100 abstracts) y validación literal de cada frase; índice BM25.
- ✔ 0 relaciones sin PMID; las frases que no aparecen literalmente en su abstract se descartan y se cuentan.

**M4 Grafo de conocimiento.** Nodos Gene, Protein, Pathway, Drug, Disease (MONDO:0005101), Dataset, Publication,
Chunk, Evidence, y más adelante CellType, CellCommunication y PathwayActivity. Toda relación pasa por un Evidence.
Enriquecimiento externo: Open Targets (EFO_0000729), ChEMBL, STRING, UniProt, Reactome.
- ✔ Recuentos por tipo; `claims_without_resolvable_evidence = 0`.

**M5 Servidores MCP.** `expression-mcp` (DE por cohorte, metaanálisis, replicación), `kg-mcp`, `external-mcp`, con
caché en disco y salida con evidence_id.
- ✔ Un test por herramienta, más una prueba stdio `list_tools` por servidor.

**M6 Agentes y orquestador.** Especialistas Sonnet (expresión, literatura, externo) y orquestador `/informe` (Opus) con
el flujo completo del renal: integración determinista, `verify.py agents`, narrativa, `coverage.py` (paso 4b),
ensamblado, `verify.py final`. Salida estructurada y hallazgos clave desde el principio.
- ✔ `claude -p "/informe"` genera `findings.json`; `verify.py final` 0 afirmaciones sin soporte y 0 desfasadas;
  `coverage.py` 0 hallazgos clave sin citar ni justificar.

**M7 Informe web.** La plantilla del renal con el menú "Sections" por categorías; volcanos, metaanálisis (forest plot
por gen), tabla de convergencia, red, fármacos, procedencia, limitaciones y apéndice de evidencia.
- ✔ Un único HTML que abre sin servidor; el validador encuentra 0 citas rotas; sin errores de consola ni scroll
  horizontal a 390 px.

**M8 Evaluación.** (a) Dianas clínicas de CU (Open Targets / ChEMBL con precedencia clínica en CU) frente al ranking:
AUROC, recall@100 y comparación con un ranking solo por tractabilidad. Controles esperados: TNF, IL12B, IL23A, IL23R,
ITGA4, ITGB7, MADCAM1, JAK1, TYK2, S1PR1, S1PR5 y TNFSF15 (TL1A). (b) Holdout temporal: dianas cuyo primer ensayo en CU
es posterior a un corte fijado antes de medir. (c) 10 preguntas con respuesta conocida contestadas con `claude -p` +
MCP. (d) Historial de decisiones en `evals/history.md`; ningún peso se cambia después de ver resultados sin registrarlo.
- ✔ `evals/results.md` con (a)–(d).

**M9 Reproducibilidad.** `uv run <cli> run-all` desde cero (descargas con sha256, análisis, literatura congelada,
grafo, convergencia, verificación, informe) con recuentos idénticos; entradas congeladas en `reproducibility/`.
- ✔ Reconstrucción limpia con salida 0 y las mismas cifras; README como caso de estudio.

**M10 Publicación.** Repo en GitHub (privado hasta que el usuario decida), historial sin abstracts ni matrices crudas.
- ✔ `git log --all` sin esos archivos; `gh repo view` correcto.

### Bloque B: capas celulares, espaciales y de descubrimiento

**M11 Células individuales y espacio.**
- Células: Smillie 2019 (SCP259: 366.650 células, 18 CU y 12 sanos; epitelio, estroma e inmune) y un segundo conjunto
  de CU (buscarlo en CELLxGENE). SCP259 puede pedir sesión de Google para descargar; si es así, pedirlo al usuario.
- Pseudobulk por donante; qué tipo celular produce cada gen; replicación de la atribución en el segundo conjunto.
- Espacio: atlas CosMx/Xenium de EII (GSE312415, GSE312420) y Xenium de colitis (GSE282123). Vecindarios celulares de
  los candidatos con permutaciones y semilla fija.
- Nuevo agente `agente-celular` y `cell-mcp`.
- ✔ Atribución esperada: LCN2 y DUOX2 epiteliales, OSM mieloide e IL11 en fibroblastos inflamatorios; replicación de la
  atribución en el segundo conjunto; resultados espaciales idénticos con otra lista de genes.

**M12 Composición celular** (equivalente a la pureza tumoral). Deconvolución de las cohortes bulk con la referencia de
M11. Para cada gen: ¿el cambio persiste ajustando por la composición, o lo explica la entrada de células inflamatorias?
- ✔ Los marcadores de neutrófilos y linfocitos que suben en bulk se marcan como "composición"; el cambio de los genes
  epiteliales de control persiste tras el ajuste.

**M13 Descubrimiento de candidatos nuevos.** Mismo cribado que el renal: suben en varias cohortes, tipo celular
dominante, pocas publicaciones en CU, sin fármaco clínico, tractables, poca expresión en órganos vitales (HPA y GTEx
proteoma), y sobreviven al ajuste de composición. Listas por compartimento (epitelio, estroma, inmune).
- ✔ Recuentos por estado; los controles positivos de M8 no aparecen como "nuevos".

### Bloque C: capas independientes y vías

**M14 Genética e independencia.**
- GWAS de EII y CU vía Open Targets: gen causal más probable (L2G) y colocalización con eQTL; si es posible, con eQTL
  de colon.
- Proteína: aleatorización mendeliana con pQTL de plasma a partir de estadísticos resumen públicos. La proteómica de
  tejido de CU es escasa (estudios pequeños); se dice explícitamente y solo se usa como apoyo.
- Actividad de la enfermedad (sustituye a la supervivencia): correlación con la puntuación endoscópica y la histología
  en las cohortes que la tengan.
- Perturbación (opcional): cribados CRISPR publicados en linfocitos T humanos.
- Convergencia v3 continua, perfiles separados de diana y de biomarcador, y sensibilidad a los pesos con 1.000
  reponderaciones.
- ✔ 0 evidencias sin resolver; la evaluación de M8 se repite con v3 y se registra la diferencia.

**M15 Capa de vías** (lo nuevo frente al renal).
- Actividad de vías (PROGENy) y de factores de transcripción (CollecTRI) por muestra y por tipo celular.
- Comunicación entre células (LIANA): pares ligando-receptor de inflamada frente a sana, con su confirmación espacial.
- Ranking de vías: una vía puntúa cuando varios de sus genes convergen en expresión, genética, tipo celular y respuesta
  a fármacos.
- `pathway-mcp` y nuevo agente `agente-vias`.
- ✔ La vía IL-23/Th17, el TNF y el eje de integrinas de los linfocitos intestinales aparecen entre las vías principales.
  Cada vía del informe cita sus registros de evidencia.

**M16 Respuesta a tratamientos.**
- UNIFI: ustekinumab frente a placebo, con interacción gen × brazo (como IMmotion en el renal); PURSUIT: golimumab;
  infliximab y vedolizumab.
- Genes y vías asociados a la respuesta en cada fármaco, siempre como asociación con su p ajustado.
- Vías que siguen activas en quienes no responden.
- `agente-aplicado` adaptado: tratamientos, composición inmune y seguridad en órganos vitales para anticuerpos
  bloqueantes o que eliminan células.
- ✔ Se reproduce al menos una firma de respuesta publicada (por ejemplo la de infliximab en GSE16879); se dice qué
  genes no separan los brazos.

### Bloque D: robustez y cierre

**M17 Réplicas y benchmark de agentes.** Tres ejecuciones independientes de `/informe` (Jaccard de la lista de dianas)
y el benchmark de agente único frente a multiagente con el mismo protocolo que el renal (calidad, coste y tiempo),
con salida estructurada y paso de cobertura.
- ✔ 0 afirmaciones sin soporte en todas las ejecuciones; Jaccard registrado; tabla de coste y tiempo.

**M18 Informe completo y limpieza.** Pestañas equivalentes a las del renal (Candidates, Expression, Cells & space,
Discovery, Pathways, Treatment, Genetics & validation, Literature & drugs, Network, Agent benchmark, Methods) con
explicación en lenguaje sencillo, botón "?" en cada panel y menú por categorías. Limpieza del repositorio: retirar
intermedios, ajustar `.gitignore` y comprobar que `run-all` sigue funcionando.
- ✔ Informe válido (0 citas rotas), sin errores de consola, sin scroll horizontal en móvil; `uv run pytest` en verde;
  `run-all` limpio.

## Definición de terminado
M0–M18 superados, `uv run pytest` en verde, el informe generado con datos reales, evaluado en `evals/results.md` y con
la verificación final en 0 afirmaciones sin soporte. Si se agota la cuota a mitad, guardar el estado en PROGRESS.md y
continuar desde ese milestone en la siguiente ventana.

## Fuera de alcance
Enfermedad de Crohn (posible ampliación posterior con UNITI, GSE112366), datos internos de cualquier empresa, la
presentación y cualquier afirmación de validación experimental.
