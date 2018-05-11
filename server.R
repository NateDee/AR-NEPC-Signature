# Define server logic


server <- function(input, output) {
	
	# Instructions panel:
	output$text1 <- renderText(
		{ paste(
			'<b>',
				'<ul>',
					'<li style="margin:10px 0";>This app is for analyzing how NEPC and AR signaling pathways are altered within your experiments.</li>',
					'<li style="margin:10px 0";>This app requires your RNA-seq data in FPKM/RPKM format.</li>',
					'<li style="margin:10px 0";>The gene lists for analyzing NEPC and AR signaling are adapted from Beltran et al. 2016. (https://www.ncbi.nlm.nih.gov/pubmed/26855148)</li>',
				'</ul>',	
			'</b>'
		)}
	)
	
	output$text2 <- renderText(
		{ paste(
			'<b>',
				'<ol>',
					'<li style="margin:20px 0";>Format your expression data according to specifications below.</li>',
				'</ol>',	
			'</b>'
		)}
	)
	
	output$text3 <- renderText(
		{ paste( 
		'<p>',
			'<ul>',
				'<li style="margin:10px 0";>The first row must contain headers (See Row 1 to the left).</li>',
				'<li style="margin:10px 0";>The first column (Column A) must contain gene names (gene symbols).</li>',
				'<b><li style="margin:10px 0";>The second row must contain a group number, see the color key to choose a group color label.</li></b>',
			'</ul>',	
		'</p>'
		)}
	)

	# User RPKM file
	userdatasetInput <- reactive({
		# validate( )
		inFile <- input$user_fpkm
		if (is.null(inFile))
			return(NULL)
		tbl <- as.data.frame(fread(inFile$datapath, sep="\t", header=TRUE, stringsAsFactors=FALSE))
		return(tbl)
	})
	
	# Output dbgap heatmap, set function
	output$dbgap_heatmap <- renderPlot({
		drawDbgapHeat()
	})
	
	# drawDbgapHeat function
	drawDbgapHeat <- function() {
		# dbgap data = dbgap_fpkm
		# dbgap_scores = dbgap_scores
		columns_colors = as.character(dbgap_scores[4,2:261])
		dbgap.mat = as.matrix(dbgap_fpkm[,4:263])
		row.names(dbgap.mat) = dbgap_fpkm$Gene_name
		row_group = c(rep("red", 30), rep("blue", 68)) # 30 AR genes, 68 NEPC genes
		palette = colorRampPalette(brewer.pal(9, "YlGnBu"))(n=299)
		# Run heatmap
		heatmap.2(dbgap.mat, trace="none", dendrogram="none", col = palette, Colv=FALSE, Rowv=FALSE, scale="none", ColSideColors = columns_colors, RowSideColors = row_group, breaks = seq(0,8,length=300), labCol=FALSE)
		text(x=rep(0.1352252,2), y=c(0.06615201, 0.51807700), srt=90, xpd=TRUE, adj=0, labels=c("NEPC Signature", "AR Signature"), font=2)
		legend("topright", title="Column Labels", legend=c("Benign", "Local PCa", "CRPC", "NEPC"), fill=c("green", "orange", "red", "darkred"), text.font=2, xpd=TRUE, inset=c(0,-0.2))
	}

	output$user_heatmap <- renderPlot({
		drawUserHeat()
	})
	
	drawUserHeat <- function() {
		# After sorting/filtering, need to convert to log2
		user = userdatasetInput()
		colnames(user)[1] = "Gene_name" # to match dbgap for join
		columns_colors = as.integer(user[1, 2:ncol(user)])
		# Recode column colors:
			columns_colors[columns_colors == -2] = "green"
			columns_colors[columns_colors == -1] = "orange"
			columns_colors[columns_colors == 0] = "black"
			columns_colors[columns_colors == 1] = "red"
			columns_colors[columns_colors == 2] = "darkred"
		row_group = c(rep("red", 30), rep("blue", 68)) # 30 AR genes, 68 NEPC genes
		palette = colorRampPalette(brewer.pal(9, "YlGnBu"))(n=299)
		
		user_filt = join(dbgap_fpkm[1:2], y=user, type="left")
		user_filt$Ensembl_ID = NULL
		# Replace duplicate genes with correct values, AR-, NKX3-1-, KLK3-,
			user_filt[user_filt$Gene_name == "AR-",][2:ncol(user_filt)] = user_filt[user_filt$Gene_name == "AR", ][2:ncol(user_filt)]
			user_filt[user_filt$Gene_name == "NKX3-1-",][2:ncol(user_filt)] = user_filt[user_filt$Gene_name == "NKX3-1", ][2:ncol(user_filt)]
			user_filt[user_filt$Gene_name == "KLK3-",][2:ncol(user_filt)] = user_filt[user_filt$Gene_name == "KLK3", ][2:ncol(user_filt)]
		user.mat = as.matrix(user_filt[,2:ncol(user_filt)])
		rownames(user.mat) = user_filt$Gene_name
		# Convert to log2(FPKM + 1)
		user.mat = log2(user.mat + 1)
		heatmap.2(user.mat, trace="none", dendrogram="none", col = palette, Colv=FALSE, Rowv=FALSE, scale="none", ColSideColors = columns_colors, RowSideColors = row_group, breaks = seq(0,8,length=300), labCol=FALSE)
		text(x=rep(0.1352252,2), y=c(0.06615201, 0.51807700), srt=90, xpd=TRUE, adj=0, labels=c("NEPC Signature", "AR Signature"), font=2)
	}
	
	output$dbgap_scores <- renderPlot({
		drawDbgapScores()
	})
	
	# Draw DBGAP NEPC/AR Scores line plot
	
	drawDbgapScores <- function() {
		columns_colors = as.character(dbgap_scores[4,2:261])
		plot(as.numeric(dbgap_scores[2,2:261]), type="l", col="black", ylab="Score", lwd=2, ylim=c(-0.2,1), xlim=c(10,251), cex.lab=1.5, font=2, xlab="")
		lines(as.numeric(dbgap_scores[1,2:261]), type="l", col="grey", lwd=2)
		legend(150, 0.5, legend=c("NEPC Score", "AR Score"), col=c("black", "grey"), lwd=2, cex=1.25)
		points(x=seq(1,260), y=rep(1,260), type="p", pch=15, col=columns_colors)
	}
	
	output$user_scores <- renderPlot({
		drawUserScores()
	})

	drawUserScores <- function() {
		# After sorting/filtering, need to convert to log2
		user = userdatasetInput()
		colnames(user)[1] = "Gene_name" # to match dbgap for join
		columns_colors = as.integer(user[1, 2:ncol(user)])
		# Recode column colors:
			columns_colors[columns_colors == -2] = "green"
			columns_colors[columns_colors == -1] = "orange"
			columns_colors[columns_colors == 0] = "black"
			columns_colors[columns_colors == 1] = "red"
			columns_colors[columns_colors == 2] = "darkred"
		
		user_filt = join(dbgap_fpkm[1:2], y=user, type="left")
		user_filt$Ensembl_ID = NULL
		# Replace duplicate genes with correct values, AR-, NKX3-1-, KLK3-,
			user_filt[user_filt$Gene_name == "AR-",][2:ncol(user_filt)] = user_filt[user_filt$Gene_name == "AR", ][2:ncol(user_filt)]
			user_filt[user_filt$Gene_name == "NKX3-1-",][2:ncol(user_filt)] = user_filt[user_filt$Gene_name == "NKX3-1", ][2:ncol(user_filt)]
			user_filt[user_filt$Gene_name == "KLK3-",][2:ncol(user_filt)] = user_filt[user_filt$Gene_name == "KLK3", ][2:ncol(user_filt)]
		user.mat = as.matrix(user_filt[,2:ncol(user_filt)])
		rownames(user.mat) = user_filt$Gene_name
		# Convert to log2(FPKM + 1)
		user.mat = log2(user.mat + 1)
		ar_cor_user = cor(as.numeric(ar_ref[,3]), user.mat[1:30,], method="pearson")
		nepc_cor_user = cor(as.numeric(nepc_ref[,3]), user.mat[31:98,], method="pearson")
		# Run plots of data
		plot(as.numeric(nepc_cor_user), type="o", col="black", ylab="Score", ylim=c(-0.2,1), lwd=2, font=2, cex=1.5, xlab="")
		lines(as.numeric(ar_cor_user), type="o", col="grey", lwd=2)
		legend("bottomleft", legend=c("NEPC Score", "AR Score"), col=c("black","grey"), lwd=2, xpd=TRUE, inset=c(0,-.3))
		points(x=seq(1:length(nepc_cor_user)), y=rep(1, length(nepc_cor_user)), type="p", pch=15, col = columns_colors)
	}
	
	# Create download boxes
	output$downloadDbgapHeatmap <- downloadHandler(
		filename <- function() {
			paste0("DBGAP_AR-NEPC_Signatures_heatmap.png")
		},
		content <- function(file) {
			png(file, width = 8, height = 10, unit="in", res=300)
			# drawDbgapHeat()
				# dbgap data = dbgap_fpkm
			# dbgap_scores = dbgap_scores
			columns_colors = as.character(dbgap_scores[4,2:261])
			dbgap.mat = as.matrix(dbgap_fpkm[,4:263])
			row.names(dbgap.mat) = dbgap_fpkm$Gene_name
			row_group = c(rep("red", 30), rep("blue", 68)) # 30 AR genes, 68 NEPC genes
			palette = colorRampPalette(brewer.pal(9, "YlGnBu"))(n=299)
			# Run heatmap
			heatmap.2(dbgap.mat, trace="none", dendrogram="none", col = palette, Colv=FALSE, Rowv=FALSE, scale="none", ColSideColors = columns_colors, RowSideColors = row_group, breaks = seq(0,8,length=300), labCol=FALSE)
			text(x=rep(0.1352252,2), y=c(0.06615201, 0.51807700), srt=90, xpd=TRUE, adj=0, labels=c("NEPC Signature", "AR Signature"), font=2)
			legend("topright", title="Column Labels", legend=c("Benign", "Local PCa", "CRPC", "NEPC"), fill=c("green", "orange", "red", "darkred"), text.font=2, xpd=TRUE, inset=c(0,-0.05))
			dev.off()
		}
	)
	
	output$downloadDbgapScores <- downloadHandler(
		filename <- function() {
			paste0("DBGAP_AR-NEPC_Signatures_Scores-Plot.png")
		},
		content <- function(file) {
			png(file, width = 10, height = 4, unit="in", res=300)
			drawDbgapScores()
			dev.off()
		}
	)
	
	output$downloadUserHeatmap <- downloadHandler(
		filename <- function() {
			paste0(basename(unlist(strsplit(input$user_fpkm$name, split='.txt', fixed=TRUE))[1]),
			"_AR-NEPC_Signatures_heatmap.png")
		},
		content <- function(file) {
			png(file, width = 8, height = 10, unit="in", res=300)
			drawUserHeat()
			dev.off()
		}
	)
	
	output$downloadUserScores <- downloadHandler(
		filename <- function() {
			paste0(basename(unlist(strsplit(input$user_fpkm$name, split='.txt', fixed=TRUE))[1]),
			"_AR-NEPC_Signatures_Scores-Plot.png")
		},
		content <- function(file) {
			png(file, width = 10, height = 4, unit="in", res=300)
			drawUserScores()
			dev.off()
		}
	)
}