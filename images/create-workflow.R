# Create the manuscript workflow figure from MetaWeave's main functions.
# Run from the project root after installing DiagrammeR and DiagrammeRsvg.
library(DiagrammeR)
library(DiagrammeRsvg)

workflow <- grViz('
digraph workflow {
  graph [rankdir = LR, bgcolor = white, nodesep = 0.35, ranksep = 0.55,
         fontname = Helvetica, pad = 0.25]
  node [shape = box, style = "rounded,filled", color = "#47616b",
        fillcolor = "#edf4f5", fontname = Helvetica, fontsize = 13,
        margin = "0.18,0.12"]
  edge [color = "#47616b", arrowsize = 0.75]

  ranges [label = "Species distribution rasters", fillcolor = "#f2f2ed"]
  sites [label = "Local species lists", fillcolor = "#f2f2ed"]
  assembly [label = <<FONT FACE="Courier New"><B>assemble_communities()</B></FONT><BR/><FONT POINT-SIZE="11">Rasters to local assemblages</FONT>>]
  direct [label = <<FONT FACE="Courier New"><B>new_assemblage()</B></FONT><BR/><FONT POINT-SIZE="11">Create one assemblage</FONT>>]
  inputs [label = "Regional interactions, blocks,\nor a user prediction function",
          fillcolor = "#f2f2ed"]
  models [label = <<FONT FACE="Courier New"><B>probability_matrix_model()</B></FONT><BR/><FONT FACE="Courier New"><B>block_model()</B></FONT><BR/><FONT FACE="Courier New"><B>new_inference_model()</B></FONT><BR/><FONT POINT-SIZE="11">Define an interaction model</FONT>>]
  infer [label = <<FONT FACE="Courier New"><B>infer_network()</B></FONT><BR/><FONT FACE="Courier New"><B>infer_networks()</B></FONT><BR/><FONT POINT-SIZE="11">Apply a model locally</FONT>>]
  matrices [label = "Local interaction matrices\nSpecies pairs and predicted values",
            fillcolor = "#e3f0e7"]
  summary [label = "Network summaries and maps\nor external analyses",
           fillcolor = "#e3f0e7"]
  spatial [label = <<FONT FACE="Courier New"><B>run_spatial_inference()</B></FONT><BR/><FONT POINT-SIZE="11">Combines the spatial steps</FONT>>,
           fillcolor = "#fff3d6", style = "rounded,filled,dashed"]

  ranges -> assembly -> infer
  sites -> direct -> infer
  inputs -> models -> infer
  infer -> matrices -> summary
  ranges -> spatial [style = dashed]
  models -> spatial [style = dashed]
  spatial -> summary [style = dashed]
}
')

# Replace the manuscript SVG with a portable vector figure.
writeLines(export_svg(workflow), "images/metaweave-workflow.svg", useBytes = TRUE)
