library(shiny)
library(RColorBrewer)
library(gplots)
library(data.table)
library(DT)

pageWithSidebar(
	titlePanel("Yu Lab - NEPC/AR Signature Score"),

	sidebarPanel(
		fileInput("user_fpkm", h4("Upload RNA-Seq data in RPKM/FPKM")),
		h3("Download Data:"),
		downloadButton("downloadDbgapHeatmap", HTML("Download<br/>DBGAP Heatmap")),
		downloadButton("downloadDbgapScores", HTML("Download DBGAP<br/>NEPC/AR Scores Plot")),
		downloadButton("downloadUserHeatmap", HTML("Download Heatmap<br/>for Your Dataset")),
		downloadButton("downloadUserScores", HTML("Download Your<br/>NEPC/AR Scores Plot")),
		width=3
	),
		
	mainPanel(
		tabsetPanel(
			# Instructions Tab
			tabPanel(
			"Instructions",
			h3(strong("NEPC/AR Signature Instructions:"), align="center"),
			htmlOutput("text1"),
			htmlOutput("text2"),
			div(style="display: inline-block;vertical-align:middle; width: 4%;", ""),
				div(style="display: inline-block;vertical-align:top; width: 35%;",	
					img(src='ExampleDataset.png', height = '100%', width='100%')),
				div(style="display: inline-block;vertical-align:top; width: 10%;",	
					img(src='ColorKey.png', height = '100%', width='100%')),
				div(style="display: inline-block;vertical-align:top; width: 45%;",
					htmlOutput("text3"))
			),
			
			# Main panel with all data
			tabPanel(
				"Results", 
				h3(strong("Comparison between human tissue (DBGAP) and your data"), align = "center"),
				fluidRow(splitLayout(cellWidths=c("50%", "50%"), h3("DBGAP Dataset", align="center"), h3("Your Uploaded Data", align="center"))),
				fluidRow(splitLayout(cellWidths=c("50%","50%"), plotOutput("dbgap_heatmap"), plotOutput("user_heatmap"))),
				fluidRow(splitLayout(cellWidths=c("50%","50%"), plotOutput("dbgap_scores"), plotOutput("user_scores")))
			)
			
			# Get working plots then do output
		)
	)
)