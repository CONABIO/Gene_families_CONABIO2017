
# Instrucciones Día 1

Entra a [http://keri.conabio.gob.mx/dokuwiki/doku.php?id=guia](http://keri.conabio.gob.mx/dokuwiki/doku.php?id=guia).

1) Sigue las instrucciones para conectarte al servidor keri como usuario `tlaloc` en dos pestañas de tu terminal.

La primera pestaña será tu ventana privada, la segunda la compartiremos con el comando:

`tmux at -t G`

2)

Revisa la viggneta [GeneFamilies.Rmd](https://github.com/asishallab/GeneFamilies/blob/master/vignettes/GeneFamilies.Rmd).

# Instrucciones Día 2

Formar grupos para hacer sacrificios humanos para Tezcatlipoca y su madre, que se llama Tezcatli-poca-madre.

Trabajos para cada grupo:

* En un `HOME` de un usuario `gen[0-21]` del grupo clonar el repositorio git `https://github.com/asishallab/GeneFamilies`. Allí crear un git branch llamado como la chamba (seleccionar de abajo) de su grupo y trabajar _solo_ allí. 

## Chambas para repartir entre grupos

### Grupo cargar familias 
* Requisitos: R and command line
* dificultad: básica 

Receta:
* Cargar las familias de genes. Adaptar el script `exec/loadGeneFamilies.R` para formar un `data.frame` que tenga este formato: 
    ```
    id        ath cla cme csa fve
    cluster_1 78  79  77  83  66
    ```
    Guardar _todos_ los objetos creado en el script en una imagen binaria de R.
* Leer, entender u modificar el script `exec/exec/loadInterProData.R` con el final de cargar las anotaciónes de InterPro. Crear un objeto `all.ipr` conteniendo todos `ath.ipr`, ..., `fva.ipr` concatonados. Formato:
    ```
    V1        V2
    AT1G17600.1 IPR000157
    ```
* Guardor aquellos objetos en formato binario en `data`
* Documentar todo en la vignette

#### Grupo anotación de familias con función
_Depende de los resultados del grupo "cargar familias"_
* requisitos: R
* dificultad: básica

Usar el paquete `AHRD.on.gene.clusters` para anotar las families con descripciones:
* Leer dentro de `R` `help( 'AHRD.on.gene.clusters-package' )`
* Bajar la version mas actuál de `interpro.xml` 
* Crear un script en `exec` que anota las familas y guarda los resultados en una columna adicional de `families.df`. Guardarlo en una imagen binaria de R en `data`

#### Grupo filogenias de familias
_Depende de los resultados del grupo "cargar familias"_
* requisitos: R y un poco de shell
* dificultad: básica

Usar las funciones del paquete `GeneFamilies` para correr para cada familia el "phylogeny pipeline". Es decir leer, entender, y usar el script `exec/reconstructFamilyPhylogenies.R`. Meter cada mil reconstrucciones en un "job" de la lista de espera (SLURM). Alicia y Asis ayudarán con el punto de SLURM. 

### Grupo Ortólogos
* Requisitos: Shell y un poquito de R
* Dificultad: avanzada

Seguir la vignette `GeneFamilies`: 
* Realizar búsquedas bidireccionales entre cada pareja de genomas, con Arabidopsis siendo el de anclage.
* Filtrar y concatonar los resultados 
* Cargar los resultados en R usando un script

#### Grupo Species Tree
_Depende de los resultados del grupo "Ortólogos"_
* Requisitos: R y shell
* Dificultad: avanzada

Receta:
* Alienar las secuencias de amino ácidos de cada grupo de proteínas ortólogas con Mafft. _Ojo_: Reemplazar los nombres de los genes con los nombres de la especia antes de correr Mafft. Siempre tener el mismo orden de las especia en cada alineamiento (MSA). Crean un script en `exec` que genera los input para mafft.
* Concatenar los MSA con uso de `csplit` y `cat`. Para eso lo del orden anterior.
* Correr FastTree
* Calibrar el resultado a que sea un Ultrametric tree. Usar el ejemplo al final de este documento ("Ultrametric tree example code"). Para esto editar el archivo `DESCRIPTION` en el paquete `GeneFamilies` y añadir el paquete `phangorn` como nuevo requisito.
* Documentar todo en la vignette

### Grupo Tandems
_Depende de los resultados del grupo "Ortólogos" y "cargar familias"_
* requisitos: R y un poco de shell
* dificultad: mediana

Receta:
* Identificar las funciones que se requieren para identificar los tandems dentro de cada familia. Buscan `tandemsFromFamiliesAndNeighborhood` y en `R`: `?tandemsFromFamiliesAndNeighborhood`
* Seguir la vignette en como preparar los GFF3 para obtener de allí la información sobre vecindad de genes. 
* No olvidar al final de obtener Tandems de cada familia hacer el `setdiff` con los ortólogos.
* Crear un script en `exec` que realice el trabajo y que guarde un archivo binario con los resultados en `data`.
* Crear tres objetos en R: 
    * `tandems`, formato:
    ```
                 Family   Gene
    tandem_cluster_3264 894950
    tandem_cluster_3264 894951
    tandem_cluster_3265 898687
    tandem_cluster_3265 898691
    tandem_cluster_3266 898693
    tandem_cluster_3266 898689
    ```
    * `tandems.genes` un vector de los nombres de genes que son tandems de todos los genomas
    * `tandems.lst` una `list` con nombres como `tandem_cluster_12` y values vectores de genes perteneciéndo a ese tandem cluster
* Ser chingón

### Grupo expansion y contracción
* requisitos: buen shell y un poco de R
* dificultad: mediana

Receta:
* Crear un script en `exec` que genera los archivos input para `cafe`.
* Correr cafe (puede ser que tengan que hacerlo desde el cafe shell - es muy mamón ese café, pero su shell está estable) Leán abajo la sección "Example CAFE Shell or Script"
* modificar script `exec/parseCafeResult.R` 
* documentar todo en el vignette

### Grupo selección 
* requisitos: R y shell
* dificultad: mediana

Por el momento empezamos nada mas con medir `Ka/Ks` en pares de gene homologos de máxima similitud de secuencia. 
* En la Vignette leer el parágrafo "Detect pairwise Ka/Ks ratios within closest homologs"
* Leer, entender y modificar script `exec/computePairwiseKaKsRatios.R`
* Usarlo en modo batches, cada conjunto de tamaño 300 pares de genes (o algún tamaño conveniente).
* Documentar todo en la vignette

# Ultrametric tree example code

```
require(phangorn)

message("USAGE: Rscript generate_ultrametric_tree.R")

input.args <- "./spe8_ortho_concate_ml_gamma_rescaled.newick"

tree.rooted <- ape::root(read.tree(input.args), outgroup = "A.arabicum", resolve.root = TRUE)

#' We know the divergence times between the following clades based on fossil
#' data and other estimates [see references at bottom of this script]:
aly.ath.mrca <- mrca.phylo(tree.rooted, which(tree.rooted$tip.label %in% c("A.thaliana",
    "A.lyrata")))
chi.aly.ath.cru.mrca <- mrca.phylo(tree.rooted, which(tree.rooted$tip.label %in%
    c("A.lyrata", "A.thaliana", "C.rubella", "C.hirsuta")))
bra.esa.mrca <- mrca.phylo(tree.rooted, which(tree.rooted$tip.label %in% c("B.rapa",
    "E.salsugineum")))
ingroup.mrca <- mrca.phylo(tree.rooted, which(tree.rooted$tip.label %in% c("A.thaliana",
    "E.salsugineum")))
aet.split <- mrca.phylo(tree.rooted, which(tree.rooted$tip.label %in% c("A.arabicum",
    "A.thaliana")))
chr.df <- data.frame(node = c(aly.ath.mrca, chi.aly.ath.cru.mrca, bra.esa.mrca, ingroup.mrca,
    aet.split), age.min = c(13, 35.6, 38.4, 43.2, 45), stringsAsFactors = FALSE)
chr.df$age.max <- chr.df$age.min
chr.df[[which(chr.df$node == aet.split), "age.max"]] <- 60
chr.df$soft.bounds <- FALSE

#' Scale the rooted tree
ctrl <- chronos.control(nb.rate.cat = 1)
tree.rooted.chrono <- chronos(tree.rooted, calibration = chr.df, lambda = 3.2)  # Lambda value obtained from references
#' Round the edge lengths:
tree.rooted.chrono$edge.length <- round(tree.rooted.chrono$edge.length, digits = 1)
#' Save the results
write.tree(tree.rooted.chrono, "./spe8_ortho_concate_ml_gamma_rescaled_ultrametric.newick")

#' Plot the tree:
pdf( "./spe8_ortho_concate_ml_gamma_rescaled_ultrametric.pdf" )
plot.phylo(tree.rooted.chrono, show.node.label=TRUE)
edgelabels(tree.rooted.chrono$edge.length)
axisPhylo()
dev.off()

png( "./spe8_ortho_concate_ml_gamma_rescaled_ultrametric.png" )
plot.phylo(tree.rooted.chrono, show.node.label=TRUE)
edgelabels(tree.rooted.chrono$edge.length)
axisPhylo()
dev.off()

#' REFERENCES:
#' [1] http://www.mobot.org/mobot/research/apweb/orders/brassicalesweb.htm
#' [2] Yang, Ruolin, David J. Jarvis, Hao Chen, Mark Beilstein, Jane Grimwood,
#'     Jerry Jenkins, ShengQiang Shu, et al. “The Reference Genome of the
#'     Halophytic Plant Eutrema Salsugineum.” Plant Genetics and Genomics 4
#'     (2013): 46. doi:10.3389/fpls.2013.00046.
```

# Example CAFE Shell or Script
```
#!cafe
#version
#date

load -i core_gene_family.txt -t 10 -l cafe_logfile.txt -p 1
tree (((((A.thaliana:13,A.lyrata:13)1.000:8,C.rubella:21)1.000:4,C.hirsuta:25)1.000:3,(E.salsugineum:17,(B.rapa:14,S.parvula:14):3)1.000:12)1.000:19,A.arabicum:47)
lambda -l 0.00549081688799
report result_cafe
```

The input file for CAFE _needs_ to look like this, i.e. the columns need to be named `description`, `id`, and then the species abbreviations identical to those in the above tree.
```
description     id      A.thaliana      A.lyrata        C.rubella       C.hirsuta       A.arabicum      B.rapa  E.salsugineum   S.parvula
cluster_1       cluster_1       78      79      77      83      66      149     75      76
cluster_2       cluster_2       70      106     67      110     1       47      46      10
cluster_3       cluster_3       39      40      56      93      32      68      34      32
```
