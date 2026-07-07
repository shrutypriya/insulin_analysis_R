#STEP 1:Loading the tool (CRITICAL: This must run every single time)
library(rentrez)
library(bio3d)
#STEP 2:Set the IDs
human_insulin_id <- "P01308" 
mouse_insulin_id <- "P01325"

#STEP 3:Fetch the datas from NCBI
human_fasta <- entrez_fetch(db="protein", id=human_insulin_id, rettype="fasta")
mouse_fasta <- entrez_fetch(db="protein", id=mouse_insulin_id, rettype="fasta")

#STEP 4:Create the folder and save the files
dir.create("data", showWarnings = FALSE) 
write(human_fasta, file="data/human_insulin.fasta")
write(mouse_fasta, file="data/mouse_insulin.fasta")

#STEP 5:Confirm it worked
cat("Success! NCBI sequence data downloaded and saved to your 'data' folder.")
# STEP 6:Set the IDs
human_insulin_id <- "P01308" 
mouse_insulin_id <- "P01325"

#STEP 7:Fetch the data from NCBI
human_fasta <- entrez_fetch(db="protein", id=human_insulin_id, rettype="fasta")
mouse_fasta <- entrez_fetch(db="protein", id=mouse_insulin_id, rettype="fasta")

#STEP 8:Create the folder and save the files
dir.create("data", showWarnings = FALSE) 
write(human_fasta, file="data/human_insulin.fasta")
write(mouse_fasta, file="data/mouse_insulin.fasta")

#STEP 9:Confirm it worked
cat("Success! NCBI sequence data downloaded and saved to your 'data' folder.")

#Load the structural library explicitly
library(bio3d)

#Download the 3D structure file (.pdb) for Human Insulin into your data folder
pdb_file <- get.pdb("1MSO", path="data", overwrite=TRUE)

#Read the downloaded structure into R so we can work with it
insulin_structure <- read.pdb("data/1MSO.pdb")

#Print a text summary of the protein structure to your console
print(insulin_structure)

#Select only the "Alpha Carbon" atoms (the central backbone of the protein)
calpha_atoms <- atom.select(insulin_structure, elety="CA")

#Generate a professional line plot showing the flexibility (B-factor) of each residue
plot.bio3d(insulin_structure$atom$b[calpha_atoms$atom], 
           sse=insulin_structure, 
           typ="l", 
           ylab="B-Factor (Flexibility)", 
           xlab="Amino Acid Residue Number", 
           main="Insulin Structural Flexibility Analysis")

#Load the visualization library
library(NGLVieweR)

#Render the interactive 3D model of your insulin structure
NGLVieweR("1MSO") %>%
  stageParameters(backgroundColor = "white") %>%
  addRepresentation("cartoon", param = list(name = "cartoon", color = "residueindex"))

#Install BiocManager (the app store for biology packages in R)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

#Install the 'msa' package using BiocManager
BiocManager::install("msa", update = FALSE)

#Load your libraries
library(rentrez)
library(Biostrings) # Automatically installed with msa to handle sequences
library(msa)

#Define the NCBI Accession IDs for our 4 species
species_ids <- c(
  Human = "P01308",
  Mouse = "P01325",
  Chicken = "P01340",
  Zebrafish = "P68082"
)

#Fetch the sequences in FASTA format from NCBI
raw_sequences <- entrez_fetch(db="protein", id=species_ids, rettype="fasta")

#Save this multi-sequence file into your 'data' folder
write(raw_sequences, file="data/insulin_4_species.fasta")

cat("Successfully downloaded all 4 species to data/insulin_4_species.fasta!")

#Read the FASTA file into R
my_protein_seqs <- readAAStringSet("data/insulin_4_species.fasta")

#Give the sequences clean names so your chart looks professional
names(my_protein_seqs) <- c("Human", "Mouse", "Chicken", "Zebrafish")

#Run the alignment using the famous 'ClustalW' algorithm
insulin_alignment <- msa(my_protein_seqs, method="ClustalW")

#Print the alignment to see it in console!
print(insulin_alignment, show=c("alignment"))

#Save the printed alignment matrix directly to the data folder
sink("data/insulin_alignment_results.txt")
print(insulin_alignment, show=c("alignment"))
sink()

cat("Alignment results saved to data/insulin_alignment_results.txt!")


