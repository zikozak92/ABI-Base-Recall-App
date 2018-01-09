library("shiny")
library("shinyBS")
library("Biostrings")
library("sangerseqR")
source("mainFunctions.R")
# ui.R

shinyUI(fluidPage(theme = "bootstrap.css",
  titlePanel(" "),

  sidebarLayout(
    sidebarPanel(
      conditionalPanel(
        condition = "!$('html').hasClass('shiny-busy')",
        tags$img(src="logoABI2.png", alt="ABI Base Recall", width="90%", class="logo")),
      
      conditionalPanel(
        condition = "$('html').hasClass('shiny-busy')",
        tags$img(src="ajax-loader.gif", alt="Processing in progress", width="43%", class="logo")),
        #tags$img(src="ssbmgif.gif", alt="Processing in progress", width="47%", class="logo")),
        #tags$img(src="giphy.gif", alt="Processing in progress", width="76.5%", class="logo")),
        #tags$img(src="tumblr.gif", alt="Processing in progress", width="58.7%", class="logo")),
      hr(),
      
      h4(tags$span(style="color: blue","File Input")),
      
      conditionalPanel(condition = "!input.example",
                       helpText(tags$small("The file uploaded must be in Abi Format (*.ab1)")),
                       fileInput('files', 'Select Abi file',
                                 multiple = FALSE,
                                 accept = c('text/csv', '.ab1'))
      ),
      
      checkboxInput('example', 'Load Example Data', FALSE),
      
      bsAlert(anchorId = "seq_corrup"),
      
      h4(tags$span(style="color: blue",'Processing')),
      
      conditionalPanel(condition = "output.fileUploaded",
                       
                       h5('End Clipping'),
                       conditionalPanel(condition = "!input.autoTrim",
                                        
                                        tags$div(class="float-left",numericInput("trim5", label = "5'",  value = "0", min = 0)),
                                        tags$div(class="float-right",numericInput("trim3", label = "3'",  value = "0", min = 0)),
                                        tags$br(style="clear:both"),
                                        bsAlert(anchorId = "long_trim"),
                                        bsAlert(anchorId = "na_trim"),
                                        helpText(tags$small("Specify the number of bases to remove from sequence's ends"))
                       ),
                       checkboxInput('autoTrim', 'Auto Trim', FALSE),
                       helpText(tags$small("Check to perform an automatic trimming")),
                       checkboxInput('recall', 'Recall N peaks', FALSE),
                       helpText(tags$small('Check for automatic correction of the N bases'))
      ),
        
        h4(tags$span(style="color: blue",'Chromatogram Options')),
      
        conditionalPanel(condition = "output.fileUploaded",
                         checkboxInput('showRmBases', 'Show Removed Bases', FALSE),
                         helpText(tags$small("Check to show the bases trimmed from the ends")),
                         conditionalPanel(condition = "input.recall",
                                          checkboxInput('showOrgBases', 'Show Original Bases', FALSE),
                                          helpText(tags$small("check to show the original sequence along with the modified one"))),
                         tags$div(class="float-left",numericInput("x", label = "Bases per Line", value = "100", min = 10, max = 400, step = 10)),
                         tags$br(style="clear:both"),
                         bsAlert(anchorId = "na_width"),
                         helpText(tags$small("Enter an appriximative number of bases per row to display in the chromatogram"))
        ),
      
     h4(tags$span(style="color: blue","Alignment")),
      
      conditionalPanel(
        condition = "output.fileUploaded",
        
        h5('Reference Sequence:'),
        
        conditionalPanel(condition = "!input.example",
                         
                        tags$textarea(id="refseq", rows=8, cols=50, ""),
                        fileInput('refFile', 'Select Fasta file',
                                  multiple = FALSE,
                                  accept = c(
                                    'text/csv',
                                    '.fasta')
                        )),
        
        conditionalPanel(condition = "input.example",
                         tags$textarea(id="refseqsample", rows=8, cols=50,
                                       "TGGGGAATCCTGCTGAACCAAGCCTTATGATCGACGGAATTCTGTGGGAAGGTTTCGGCGGAGATCCTTGCGATCCTTGC
ACCACTTGGTGTGACGCTATCAGCATGCGTATGGGTTACTATGGTGACTTTGTTTTCGACCGTGTTTTGAAAACAGATGT
GAATAAAGAATTCCAAATGGGTGACAAGCCTACAAGTACTACAGGCAATGCTACAGCTCCAACCACTCTTACAGCAAGAG
AGAATCCTGCTTACGGCCGACATATGCAGGATGCTGAGATGTTTACAAATGCCGCTTGCATGGCATTGAATATTTGGGAT
CGCTTTGATGTATTCTGTACACTAGGAGCCTCTAGCGGATACCTTAAAGGAAACTCTGCTTCTTTCAATTTAGTTGGATT
GTTTGGAGATAATGAAAATCAAAGCACGGTCAAAACGAATTCTGTACCAAATATGAGCTTAGATCAATCTGTTGTTGAAC
TTTACACAGATACTGCCTTCTCTTGGAGCGTGGGCGCTCGAGCAGCTTTGTGGGAGTGCGGATGTGCGACTTTAGGGGCT
TCTTTCCAATACGCTCAATCTAAACCTAAAGTCGAAGAATTAAACGTTCTCTGTAACGCAGCTGAGTTTACTATCAATAA
GCCTAAAGGATATGTAGGGCAAGAATTCCCTCTTGCACTCATAGCAGGAACTGATGCAGCGACGGGCACTAAAGATGCCT
CTATTGATTACCATGAGTGGCAAGCAAGTTTAGCTCTCTCTTACAGATTGAATATGTTCACTCCCTACATTGGAGTTAAA
TGGTCTCGAGCAAGTTTTGATGCCG")
        ),

        helpText(tags$small("Sequence should include at least the region covered by the sequence results. Non-DNA characters will automatically be removed.")),
        
        bsCollapse(multiple = FALSE, id = "collapse1",
                   bsCollapsePanel("Alignment Type",
                                   radioButtons("align_type", "",
                                                choices=c("Global" = "global",
                                                          "Local" = "local",
                                                          "Overlap" = "overlap",
                                                          "Global-Local" = "global-local",
                                                          "Local-Global" = "local-global"),
                                                selected = "global" ))
        ),
        
        helpText(tags$small("Choose the alignment type to be performed, set on Global by default.")),

        selectInput("blastProg", "BLAST Programs", 
                    choices = c("Nucleotide, BLASTn", 
                                "Translated, BLASTx", 
                                "Translated, TBLASTx", 
                                "MegaBLAST")),
        
        uiOutput("blastButt"),

        helpText(tags$small("Click on Blast to load the sequence in the selected blast program."))
      )
    ),
    
    mainPanel(
      
      tabsetPanel(id="maintabset", selected="Instructions",
                  
                  tabPanel("Chromatogram",
                           conditionalPanel(
                             condition = "!output.fileUploaded", 
                             tags$p("A Chromatogram plot will be shown here when a valid input has been uploaded.")
                           ), 
                           conditionalPanel(
                             condition = "output.fileUploaded",
                             tags$p(""),
                             uiOutput("exportPdf"),
                             hr(),
                             plotOutput("chrom.plot", height="auto")
                           )
                  ),
                  
                  tabPanel("Sequence Info",
                           conditionalPanel(
                             condition = "!output.fileUploaded", 
                             tags$p("Sequence set of information will be shown here when a valid input has been uploaded.")
                           ), 
                           conditionalPanel(
                             condition = "output.fileUploaded",
                             h5('Sequence:'),
                             verbatimTextOutput("seqShow"),
                             uiOutput("exportFasta"),
                             hr(),
                             h5('Quality Histogram :'),
                             plotOutput("qualityHist"),
                             h5('Cumulative Quality Histogram:'),
                             plotOutput("cumQualityHist"),
                             hr(),
                             h5('Details:'),
                             uiOutput("detailsWidget")
                           )
                  ),
                  
                  tabPanel("Alignment",
                           conditionalPanel(
                             condition = "!output.refSeqUploaded", 
                             tags$p("pairwise alignment will be shown here when a valid reference sequence has been selected.")
                           ), 
                           conditionalPanel(
                             condition = "output.refSeqUploaded",
                             checkboxInput('reverseSeq', 'Reverse the main sequence.', FALSE),
                             helpText(tags$small("check to reverse the sequence for the alignment.")),
                             hr(),
                             h5('Alignment Header:'),
                             verbatimTextOutput("AlignHead"),
                             h5('Alignment:'),
                             verbatimTextOutput("AlignDispl"),
                             h5('Alignment Program:'),
                             verbatimTextOutput("AlignProg")
                           )
                  ),
                  
                  tabPanel("Report",
                           conditionalPanel(
                             condition = "!output.fileUploaded", 
                             tags$p("Options to generate a report will be shown here when a valid input has been uploaded.")
                           ), 
                           conditionalPanel(
                             condition = "output.fileUploaded",
                             uiOutput("cbReportSeq"),
                             conditionalPanel(
                               condition = "output.refSeqUploaded",
                               uiOutput("cbReportAlign")
                             ),
                             radioButtons("reportFormat", "Report Format:",
                                          choices=c('PDF', 'HTML', 'Word')),
                             hr(),
                             uiOutput("ReportFileName"),
                             downloadButton('dlreport', 'Download')
                           ) 
                  ),
                  
                  tabPanel("Multiple Correction",
                           
                           fileInput('multifiles', '',
                                     multiple = TRUE,
                                     accept = c('text/csv', '.ab1')),
                           
                           checkboxInput('multirecall', 'Recall N peaks', FALSE),
                           
                           conditionalPanel(condition = "!input.multiautoTrim",
                           numericInput("multitrim5", label = "5'",  value = "0", min = 0),
                           numericInput("multitrim3", label = "3'",  value = "0", min = 0)
                           ),
                           
                           checkboxInput('multiautoTrim', 'Auto Trim', FALSE),
                           
                           hr(),
                           
                           bsAlert(anchorId = "unprocessed"),
                           
                           uiOutput("multiexportFasta")
                           
                  ),
                  
                  tabPanel("Instructions",
                           includeHTML("instructions.html")
                  )
                  
      )
      
    )
  )
))