library(here)
for (i in c("BP", "MF", "CC")) {
  for (j in c("lnHMdisp", "voom")) {
    for (k in c("brca", "coad", "kirc", "lihc", "luad", "lusc", "prad", "thca")) {
      assign(
        paste0("full_results_", i, "_", j, "_", k), 
        read.csv(
          here(
            "Results/GSEA results June 2020", 
            paste0("full_results_", i, "_", j, "_", k, ".csv")
          )
        )
      )
      assign(
        paste0("rrvgo_redundancy_results_", i, "_", j, "_", k), 
        read.csv(
          here(
            "Results/GSEA results June 2020", 
            paste0("rrvgo_redundancy_results_", i, "_", j, "_", k, ".csv")
          )
        )
      )
      assign(
        paste0("redundant_terms_removed_", i, "_", j, "_", k), 
        read.csv(
          here(
            "Results/GSEA results June 2020", 
            paste0("redundant_terms_removed_", i, "_", j, "_", k, ".csv")
          )
        )
      )
    }
  }
}
rm(i,j,k)


## Look at significant terms at 0.05 level, up to 10

# brca
redundant_terms_removed_BP_lnHMdisp_brca[1:10, ]
# Homophilic cell adhesion via plasma membrane molecules
# Double-strand break repair
# ncRNA metabolic process
# Macromolecule methylation
# RNA 3'-end processing
# Gene silencing
# DNA recombination
# Telomere organisation
# RNA localisation
# mRNA processing
redundant_terms_removed_BP_voom_brca[1:10, ]
# Extracellular structure organisation
# Regulation of system process
# Homophilic cell adhesion via plasma membrane adhesion molecules
# Regulation of ion transmembrane transport
# Ameboidal type cell migration
# Response to corticosteroid
# Adenylate cyclase activating G-protein coupled receptor signalling pathway
# G-protein coupled receptor signalling pathway coupled to cyclic nucleotide second messenger
# Sensory organ development
# Skeletal system development
redundant_terms_removed_MF_lnHMdisp_brca[1:10, ]
# Histone binding
# Syntaxin 1 binding
# Nuclease activity
# Single-stranded DNA binding
# Transferase activity transferring one-carbon groups
# Chromatin binding
# Basal transcription machinery binding
# RNA polymerase core enzyme binding
# Ubiquitin-like protein ligase activity
# Telomeric DNA binding
redundant_terms_removed_MF_voom_brca[1:10, ]
# Glycosaminoglycan binding
# Growth factor activity
# G-protein coupled receptor activity
# Integrin binding
# Hormone binding
# Sulphur compound binding
# Inorganic anion transmembrane transporter activity
# Calcium ion transmembrane transporter activity
# Protein heterodimerisation activity
# Retinol dehydrogenase activity
redundant_terms_removed_CC_lnHMdisp_brca[1:10, ]
# Endoplasmic reticulum quality control compartment
# Nuclear chromosome telomeric region
# Ciliary tip
# Replisome
# Nuclear pore
redundant_terms_removed_CC_voom_brca[1:10, ]
# Extracellular matrix
# Immunoglobulin complex
# Endoplasmic reticulum lumen
# Apical plasma membrane
# Protein-DNA complex
# Cell-cell junction
# Nuclear nucleosome
# Blood microparticle
# Synaptic membrane
# Collagen trimer

# coad
redundant_terms_removed_BP_lnHMdisp_coad[1:10, ]
# mRNA processing
# RNA 3'-end processing
# Translational initiation
# ncRNA metabolic process
# RNA polyadenylation
# Ribonucleoprotein complex subunit organisation
# Ribosome biogenesis
# Cotranslational protein targeting to membrane
# Establishment of protein localisation to endoplasmic reticulum
redundant_terms_removed_BP_voom_coad[1:10, ]
# Humoral immune response
# Phagocytosis recognition
# B-cell receptor signalling pathway
# Immunoglobulin production
# Adaptive immune response based on somatic recombination of immune receptors built from immunoglobulin superfamily domains
# Membrane invagination
# Defence response to bacterium
# Phagocytosis
# Drug transport
# Receptor-mediated endocytosis
redundant_terms_removed_MF_lnHMdisp_coad[1:10, ]
# Structural constituent of ribosome
# 3'-5' exonuclease activity
# Ubiquitin-like protein-specific protease activity
# Ribonucleoprotein complex binding
# rRNA binding
# Translation factor activity RNA binding
# Polyubiquitin modification-dependent protein binding
redundant_terms_removed_MF_voom_coad[1:10, ]
# Antigen binding
# Immunoglobulin receptor binding
# Cytokine activity
# Sulphur compound binding
# Cation channel activity
# Organic acid binding
# Extracellular matrix structural constituent
# G-protein coupled receptor binding
# G-protein coupled receptor activity
# Glucuronosyl transferase activity
redundant_terms_removed_CC_lnHMdisp_coad[1:10, ]
# mRNA cleavage factor complex
# Inner mitochondrial membrane protein complex
# Mitochondrial matrix
# U2-type spliceosomal complex
# Ribosome
# DNA packaging complex
# Nuclear chromosome telomeric region
# Protein-DNA complex
redundant_terms_removed_CC_voom_coad[1:10, ]
# Immunoglobulin complex
# External side of plasma membrane
# Extracellular matrix
# Receptor complex
# Plasma membrane protein complex
# Blood microparticle
# Cell body
# Anchored component of membrane
# Collagen trimer
# Postsynaptic membrane

# kirc
redundant_terms_removed_BP_lnHMdisp_kirc[1:10, ]
# ncRNA metabolic process
# Negative regulation of cell cycle G2-M phase transition
# Establishment of protein localisation to endoplasmic reticulum
# Cilium organisation
# Interleukin-1 mediated signalling pathway
# Regulation of cellular amino acid metabolic process
# Ribosome biogenesis
# Hematopoietic stem cell differentiation
# Regulation of stem cell differentiation
# Proteasomal protein catabolic process
redundant_terms_removed_BP_voom_kirc[1:10, ]
# Adaptive immune response
# Humoral immune response
# Production of molecular mediator of immune response
# Phagocytosis recognition
# Membrane invagination
# Leukocyte chemotaxis
# Extracellular structure organisation
# Regulation of body fluid levels
# Positive regulation of T-cell proliferation
# Phospholipase C-activating G-protein coupled receptor signalling pathway
redundant_terms_removed_MF_lnHMdisp_kirc[1:10, ]
# Structural component of ribosome
# Ligase activity forming carbon-oxygen bonds
# Threonine-type peptidase activity
# tRNA binding
# Translation factor activity RNA binding
# Polyubiquitin modification-dependent protein binding
# Ribosome binding
redundant_terms_removed_MF_voom_kirc[1:10, ]
# Immunoglobulin receptor binding
# G-protein coupled receptor activity
# Cytokine receptor activity
# Gated channel activity
# Carbohydrate binding
# Cytokine binding
# Glycosaminoglycan binding
# Receptor regulator activity
# Anion transmembrane transporter activity
# Peptide receptor activity
redundant_terms_removed_CC_lnHMdisp_kirc[1:10, ]
# Cilium
# Ribosome
# Proteasome core complex
# Ubiquitin ligase complex
# Plasma membrane-bounded cell projection cytoplasm
# Clathrin adaptor complex
# Clathrin coat of coated pit
# Sodium channel complex
# Trans-gogli network
# U2-type spliceosomal complex
redundant_terms_removed_CC_voom_kirc[1:10, ]
# Immunoglobulin complex
# External side of plasma membrane
# MHC protein complex
# T-cell receptor complex
# Endoplasmic reticulum lumen
# Basolateral plasma membrane
# Secretory granule membrane
# Apical part of membrane
# Basement membrane
# Blood microparticle

# lihc
redundant_terms_removed_BP_lnHMdisp_lihc[1:10, ]
# Complement activation
# Humoral immune response mediated by circulating immunoglobulin
# Homophilic cell adhesion via plasma membrane adhesion molecules
# Organic acid catabolic process
# Phagocytosis recognition
# Organic acid biosynthetic process
# Protein localisation to endoplasmic reticulum
# Drug catabolic process
# Regulation of hormone levels
# Lipid catabolic process
redundant_terms_removed_BP_voom_lihc[1:10, ]
# Humoral immune response mediated by circulating immunoglobulin
# Complement activation
# Adaptive immune response based on somatic recombination of immune receptors built from immunoglobulin superfamily domains
# Membrane invagination
# Phagocytosis
# Positive regulation of B-cell activation
# Small molecule catabolic process
# Leukocyte migration
# Receptor-mediated endocytosis
# Xenobiotic metabolic process
redundant_terms_removed_MF_lnHMdisp_lihc[1:10, ]
# Immunoglobulin receptor binding
# Antigen binding
# Cation transmembrane transporter activity
# Electron transfer activity
# Oxidoreductase activity acting on paired donors with incorporation or reduction of molecular oxygen-reduced flavin or 
# flavoprotein as one donor and incorporation of one atom of oxygen
# Serine-type endopeptidase inhibitor activity
# Vitamin binding
# Extracellular matrix structural constituent conferring compression resistance
# Carbon-carbon lyase activity
# Transferase activity transferring hexosyl groups
redundant_terms_removed_MF_voom_lihc[1:10, ]
# Antigen binding
# Immunoglobulin receptor binding
# Cofactor binding
# Organic acid binding
# Tetrapyrrole binding
# Vitamin binding
# Carbohydrate binding
# Oxygen binding
# Serine-type endopeptidase inhibitor activity
# Extracellular matrix structural constituent
redundant_terms_removed_CC_lnHMdisp_lihc[1:10, ]
# Immunoglobulin complex circulating
# Mitochondrial matrix
# Extracellular matrix
# Blood microparticle
# Organelle inner membrane
# Apical plasma membrane
# Endoplasmic reticulum lumen
# Large ribosomal subunit
# Basolateral plasma membrane
# High density lipoprotein particle
redundant_terms_removed_CC_voom_lihc[1:10, ]
# Immunoglobulin complex circulating
# Blood microparticle
# External side of plasma membrane
# Platelet alpha granule
# Vesicle lumen
# High density lipoprotein particle
# Collagen trimer
# Microbody
# Kinetochore
# Golgi lumen

# luad
redundant_terms_removed_BP_lnHMdisp_luad[1:10, ]
# Chromatin assembly or disassembly
# Glycosylation
# DNA conformation change
# Protein-DNA complex subunit organisation
# Homophilic cell adhesion via plasma membrane adhesion molecules
# Retrograde vesicle-mediated transport golgi to endoplasmic reticulum
# RNA polyadenylation
# Ethanol metabolic process
# Telomere organisation
# Glycoprotein metabolic process
redundant_terms_removed_BP_voom_luad[1:10, ]
# Phagocytosis recognition
# B-cell mediated immunity
# Lymphocyte mediated immunity
# Defence response to bacterium
# Phagocytosis
# Immunoglobulin production
# Membrane invagination
# FC receptor-mediated stimulatory signalling pathway
# Urogenital system development
# Regulation of body fluid levels
redundant_terms_removed_MF_lnHMdisp_luad[1:10, ]
# Oxidoreductase activity acting on the CH-CH group of donors
# Ligase activity forming carbon oxygen bonds
# Transferase activity transferring pentosyl groups
# Ubiquitin-like protein transferase activity
# Cofactor binding
# Ubiquitin-like protein binding
# Symporter activity
# Transferase activity transferring nitrogenous groups
# Protein heterodimerisation activity
# 3'-5' exonuclease activity
redundant_terms_removed_MF_voom_luad[1:10, ]
# Antigen binding
# Immunoglobulin receptor binding
# Receptor regulator activity
# Extracellular matrix structural constituent
# Oxygen binding
# G-protein coupled receptor activity
# Tetrapyrrole binding
# Iron ion binding
# Carbohydrate binding
# Glycosaminoglycan binding
redundant_terms_removed_CC_lnHMdisp_luad[1:10, ]
# DNA packaging complex
# Protein-DNA complex
# Transport vesicle
# Intrinsic component of endoplasmic reticulum membrane
# Coated vesicle
# Nuclear envelope
# Site of polarised growth
# Nuclear chromosome telomeric region
# Lamellar body
# Presynaptic membrane
redundant_terms_removed_CC_voom_luad[1:10, ]
# Immunoglobulin complex circulating
# Extracellular matrix
# External side of plasma membrane
# Collagen trimer
# DNA packaging complex
# Apical plasma membrane
# Blood microparticle
# Anchored component of membrane
# I band
# Basolateral plasma membrane

# lusc
redundant_terms_removed_BP_lnHMdisp_lusc[1:10, ]
# Synapse assembly
# Regulation of system process
# Dendrite development
# Cell-cell adhesion via plasma membrane adhesion molecules
# Gland morphogenesis
# Glycosylation
# Regulation of ion transmembrane transport
# Multicellular organismal signalling
# Organic hydroxy compound transport
# Regulation of exocytosis
redundant_terms_removed_BP_voom_lusc[1:10, ]
# Keratinisation
# Epidermal cell differentiation
# Extracellular structure organisation
# Regulation of body fluid levels
# G-protein coupled receptor signalling pathway coupled to cyclic nucleotide second messenger
# Chromosome segregation
# Cell cycle DNA replication
# Myofibril assembly
# Positive regulation of cytokine production
# Cell fate commitment
redundant_terms_removed_MF_lnHMdisp_lusc[1:10, ]
# Voltage-gated ion channel activity
# Ligase activity
# Phosphoric ester hydrolase activity
# Potassium ion transmembrane transporter activity
# RAC guanyl nucleotide exchange factor activity
# Anion transmembrane transporter activity
# Transferase activity transferring glycosyl groups
# Frizzled binding
# Cofactor binding
# Phospholipid binding
redundant_terms_removed_MF_voom_lusc[1:10, ]
# Carbohydrate binding
# Glycosaminoglycan binding
# Immunoglobulin binding
# Serine-type endopeptidase inhibitor activity
# Cytokine binding
# DNA helicase activity
# Extracellular matrix structural constituent
# Organic acid binding
# Cytokine receptor activity
# Tetrapyrrole binding
redundant_terms_removed_CC_lnHMdisp_lusc[1:10, ]
# Cell body
# Transport vesicle
# Plasma membrane protein complex
# Sarcolemma
# Presynaptic membrane
# Postsynaptic membrane
# Organelle subcompartment
# Cell-cell junction
# Site of polarised growth
# Intrinsic component of organelle membrane
redundant_terms_removed_CC_voom_lusc[1:10, ]
# Extracellular matrix
# Cornified envelope
# Cell-cell junction
# Apical plasma membrane
# Intermediate filament cytoskeleton
# Secretory granule membrane
# MHC protein complex
# I band
# Membrane region
# Actin cytoskeleton

# prad
redundant_terms_removed_BP_lnHMdisp_prad[1:10, ]
# None with FDR < 0.05 
redundant_terms_removed_BP_voom_prad[1:10, ]
# Anion transport
# Regulation of transmembrane transport
# Muscle filament sliding
# Regulation of muscle contraction
# Regulation of membrane potential
# Cognition
# Cyclic nucleotide-mediated signalling
# Multicellular organismal signalling
# Cellular component assembly involved in morphogenesis
# Cell fate commitment
redundant_terms_removed_MF_lnHMdisp_prad[1:10, ]
# None with FDR < 0.05
redundant_terms_removed_MF_voom_prad[1:10, ]
# Actin binding
# Voltage-gated ion channel activity
# Potassium ion transmembrane transporter activity
# Growth factor activity
# Glycosaminoglycan binding
# Integrin binding
# Secondary active transmembrane transporter activity
# Neuropeptide receptor activity
# Actinin binding
# Structural constituent of cytoskeleton
redundant_terms_removed_CC_lnHMdisp_prad[1:10, ]
# None with FDR < 0.05
redundant_terms_removed_CC_voom_prad[1:10, ]
# I band
# Basolateral plasma membrane
# Cell body
# Sarcolemma
# Actomyosin
# Excitatory synapse
# Neuron projection terminus
# Plasma membrane protein complex
# Apical part of cell
# Endoplasmic reticulum lumen

# thca
redundant_terms_removed_BP_lnHMdisp_thca[1:10, ]
# Homophilic cell adhesion via plasma membrane adhesion molecules
# Skeletal system development
# Cartilage development
# Sensory perception
# Striated muscle cell differentiation
# Sensory system development
# Renal system process
# Heart development
# Regulation of membrane potential
# Neurotransmitter transport
redundant_terms_removed_BP_voom_thca[1:10, ]
# Regulation of system process
# Extracellular structure organisation
# Anion transport
# Regulation of neuron projection development
# Bone growth
# Odontogenesis
# Positive regulation of synaptic transmission
# Neurotransmitter transport
# Sensory organ development
# Skin development
redundant_terms_removed_MF_lnHMdisp_thca[1:10, ]
# Extracellular matrix structural constituent
# Sulphur compound binding
# Glycosaminoglycan binding
# Transmembrane receptor protein tyrosine kinase activity
# Voltage-gated ion channel activity
# Cofactor binding
# Extracellular matrix binding
# Organic anion transmembrane transporter activity
# Integrin binding
# Cell adhesion mediator activity
redundant_terms_removed_MF_voom_thca[1:10, ]
# Receptor regulator activity
# Extracellular matrix structural constituent
# Glycosaminoglycan binding
# Metal ion transmembrane transporter activity
# Voltage-gated ion channel activity
# Actin binding
# Hydrolase activity acting on acid phosphorus-nitrogen bonds
# Sulphur compound binding
# Integrin binding
# Growth factor binding
redundant_terms_removed_CC_lnHMdisp_thca[1:10, ]
# Collagen trimer
# Basement membrane
# Postsynaptic membrane
# Presynaptic membrane
# Anchored component of membrane
# I band
# Endoplasmic reticulum lumen
# Dendritic tree
# Basolateral plasma membrane
# Golgi lumen
redundant_terms_removed_CC_voom_thca[1:10, ]
# Postsynaptic membrane
# Anchored component of membrane
# Collagen trimer
# Basement membrane
# Receptor complex
# Golgi lumen
# Endoplasmic reticulum lumen
# Cell-cell junction
# Main axon
# Excitatory synapse


## Assess GO term results with reference to ability of DD (differential dispersion) to 
## identify cancer-related genes

# brca
# DD identifies roughly same number of cancer-related genes as DE
# Strong association for DE, none for DD
# DD terms related to transcription, translation, gene silencing
# DE terms related to ECM, signalling

# coad
# DD identifies roughly same number of cancer-related genes as DE
# Strong association for DE, none for DD
# DD terms related to translation, protein localisation
# DE terms related to immune system, ECM, signalling

# kirc
# DE identifies clearly more cancer-related genes than DD
# Strong association for DE, none for DD
# DD terms related to cell cycle, translation
# DE terms related to immune system

# lihc
# DD identifies clearly more cancer-related genes than DE
# Strong association for DE, none for DD
# DD terms related to immune system, ECM, protein localisation
# DE terms related to immune system

# luad
# DD identifies clearly more cancer-related genes than DE
# Strong association for DE and DD
# DD terms related to translation, chromosome conformation
# DE terms related to immune system

# lusc
# DD identifies roughly same number of cancer-related genes as DE
# Strong association for DE and DD
# DD terms related to signalling, neural function
# DE terms related to ECM, cell cycle, cytoskeleton

# prad
# DE identifies far more cancer-related genes than DD
# Strong association for DE, none for DD
# No terms related to DD
# DE terms related to muscle function, morphogenesis, signalling

# thca
# DE identifies clearly more cancer-related genes than DD
# Strong association for DE and DD
# DD terms related to tissue development, ECM, signalling
# DE terms related to ECM, neural function, organ development

