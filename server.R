# server.R
library("shiny")
library("shinyBS")
library("Biostrings")
library("sangerseqR")
library("rmarkdown")
library("knitr")
source("mainFunctions.R")

options(shiny.maxRequestSize = 50*1024^2)
sample<-"Sample/MAR035_P1_2008-07-24.ab1"

shinyServer(function(input,output,session) {
  
  #--------Reactives----------#
  
  sampleName <- reactive({ 
    if(input$example) {
      return(sub("([^.]+)\\.[[:alnum:]]+$", "\\1",basename(sample)))
    } else if(!is.null(input$files)) {
      return(sub("([^.]+)\\.[[:alnum:]]+$", "\\1",unlist(input$files$name)))
    } else return(NULL)
  })
  
  objABIF <- reactive({ 
    if(input$example) {
      obj <- read.abif(sample)
      return(obj) 
    }
    else if(!is.null(input$files)) {
      obj <- read.abif(input$files$datapath)
      return(obj)
    }
    else
      return(NULL)
  })
 
  inputdata <- reactive({
    
    closeAlert(session, alertId="alert")
    if(!is.null(objABIF())){
      if (is.null(objABIF()@data$DATA.9) ||
            is.null(objABIF()@data$DATA.10) ||
            is.null(objABIF()@data$DATA.11) ||
            is.null(objABIF()@data$DATA.12) ||
            is.null(objABIF()@data$PLOC.2) ||
            is.null(objABIF()@data$PBAS.2) ||
            is.null(objABIF()@data$PCON.2)){
        createAlert(session, anchorId = "seq_corrup",
                    alertId="alert",
                    content = "The file uploaded is corrupted or doesn't contain some important informations.",
                    style = "danger",
                    append=FALSE,
                    #block = FALSE,
                    dismiss = FALSE
        )
        
        return(NULL)
      }
      else
        
        return(baseCall_fun(objABIF(), sampleName(), input$recall))
    }
    else
      return(NULL)
  })
  
  baseTrim <- reactive({
    closeAlert(session, alertId="alert2")
    closeAlert(session, alertId="alert3")
    if (input$autoTrim)
      t5t3 <- calcul.t5t3(inputdata()$Q)
    else
      
      if (is.na(input$trim5) | is.na(input$trim3)){
        createAlert(session, anchorId = "na_trim",
                    alertId="alert3",
                    content = "Can't be NA! Specify valid numbers (default : 0). Or check Auto Trim for an atomatic trimming.",
                    style = "warning",
                    append=FALSE,
                    #block = FALSE,
                    dismiss = FALSE
        )
        
        t5t3 <- list(t5=0, t3=0)
      }

    else if( length(inputdata()$seq) <= (input$trim5+input$trim3) ){
      createAlert(session, anchorId = "long_trim",
                  alertId="alert2",
                  content = paste0("Number of bases to trim is greater than the length of the sequence which equals to ",
                                   length(inputdata()$seq), ". (default : 0)"),
                  style = "warning",
                  append=FALSE,
                  #block = FALSE,
                  dismiss = FALSE
      )
      
      t5t3 <- list(t5=0, t3=0)
    }
    else
      
      t5t3 <- list(t5=input$trim5, t3=input$trim3)
    
    return(t5t3)
  })
  
  chromWidth <- reactive({
    closeAlert(session, alertId="alert4")
      if (is.na(input$x)){
        createAlert(session, anchorId = "na_width",
                    alertId="alert4",
                    content = "Can't be NA! Enter a valid number (default : 100).",
                    style = "warning",
                    append=FALSE,
                    #block = FALSE,
                    dismiss = FALSE
        )
        wid <- 100
      }
    else
      wid <- input$x
    return(wid)
  })
  
  trimmedSequence <- reactive({
    
    if (sum( baseTrim()$t5 , baseTrim()$t3 ) > length(inputdata()$seq))
      return ("")
    else
      return (inputdata()$seq[(baseTrim()$t5+1):(length(inputdata()$seq)-baseTrim()$t3)])

  })
  
  h <- reactive({
    if (!is.null(input$files) | input$example) {
      chromHeight(inputdata(), baseTrim()$t5, baseTrim()$t3, width=chromWidth(), showtrim=input$showRmBases)
    }
  })
  
  refSequence <- reactive({
    
    if(input$example){
      refseq <- cleanrefseq(input$refseqsample)
      return (list(name = "Reference",
                   width = nchar(refseq),
                   seq = refseq))
    }
    
    if (!is.null(input$refFile)){
      refseq <- readDNAStringSet(input$refFile$datapath)
      return (list(name = refseq@ranges@NAMES,
                   width = refseq@ranges@width,
                   seq = cleanrefseq(toString(refseq))))
    }
    
    if(input$refseq != ""){
      refseq <- cleanrefseq(input$refseq)
      return (list(name = "Reference",
                   width = nchar(refseq),
                   seq = refseq))
    }
    
    else
      return(NULL) 
  })
  
  alignementPw <- reactive({
    if(!is.null(trimmedSequence()) & !is.null(refSequence())){
      mainSeq<-ifelse(input$reverseSeq, char2string(rev(trimmedSequence())), char2string(trimmedSequence()))
      return(pairwiseAligner(mainSeq,refSequence()$seq,input$align_type))
    }
  })
  
  progBLASTquery <- reactive({
    if (is.null(trimmedSequence()))
      return(NULL)
    switch(input$blastProg,
           "Nucleotide, BLASTn"= blaster(trimmedSequence(),sampleName(),"BLASTn"),
           "Translated, BLASTx"= blaster(trimmedSequence(),sampleName(),"BLASTx"),
           "Translated, TBLASTx"= blaster(trimmedSequence(),sampleName(),"TBLASTx"),
           "MegaBLAST"= blaster(trimmedSequence(),sampleName(),"MegaBLAST"))
  })
  
  #-------------UIs------------#
  
  output$exportFasta <- renderUI({
    bsCollapse(multiple = FALSE, id = "collapse2",
               bsCollapsePanel("Export as Fasta",
                               textInput("seq_header", label = "Header:", value = sampleName()),
                               numericInput("seq_bpl", label = " Number of Bases per Line:", value = "60", min = 10, max = 100, step = 10),
                               textInput("seq_fname", label = "File Name:", value = sampleName()),
                               downloadButton('dlfasta', 'Download'))
    )
  })
  
  output$multiexportFasta <- renderUI({
    bsCollapse(multiple = FALSE, id = "collapse4",
               bsCollapsePanel("Export as Fasta",
                               numericInput("multiseq_bpl", label = " Number of Bases per Line:", value = "60", min = 10, max = 100, step = 10),
                               textInput("multiseq_fname", label = "File Name:", value = "Multi_fasta"),
                               downloadButton('multidlfasta', 'Download'))
    )
  })
  
  output$exportPdf <- renderUI({
    bsCollapse(multiple = FALSE, id = "collapse3",
               bsCollapsePanel("Export as PDF",
                               textInput("chrom_fname", label = "File Name:", value = sampleName()),
                               downloadButton('dlchrom', 'Download'))
    )
  })

  output$blastButt <- renderUI({
    a(actionButton("action", label = "Blast"), href=progBLASTquery(), target="_blank" )
  })
  
  output$cbReportSeq<- renderUI({
    checkboxGroupInput("variable1", "plots and informations to include in the report:",
                       c("Sequence" = "seq",
                         "Quality Histogram" = "qualHist",
                         "Cumulative Quality Histogram" = "cumQualHist"))
  })
  
  output$cbReportAlign<- renderUI({
    checkboxGroupInput("variable2", "Alignment informations to include in the report:",
                       c("Reference Sequence" = "refSeq",
                         "Alignment" = "align"))
  })
  
  output$ReportFileName<- renderUI({
    textInput("report_fname", label = "File Name :", value = sampleName())
  })
  
  output$detailsWidget<- renderUI({
    selectInput('in3', '', dumbExtractor(detailSeq(objABIF())), multiple=TRUE, selectize=FALSE,width='100%')
  })

  #--------ReactiveOutputs-------#

  output$h <- reactive(h())
  
  output$fileUploaded <- reactive({return(!is.null(inputdata()))})
  outputOptions(output, "fileUploaded", suspendWhenHidden = FALSE)
  
  output$refSeqUploaded <- reactive({return(!is.null(refSequence()))})
  outputOptions(output, "refSeqUploaded", suspendWhenHidden = FALSE)
  
  #------------Prints------------#
  
  output$seqShow <- renderPrint({ 
    if (!is.null(trimmedSequence()))
      cat(char2string(trimmedSequence()))
     })
  
  output$AlignHead <- renderPrint({ 
    if (!is.null(alignementPw()))
      cat(alignementPw()$header, sep='\n')
    })
  
  output$AlignDispl <- renderPrint({ 
    if (!is.null(alignementPw()))
      cat(alignementPw()$align, sep='\n')
    })
  
  output$AlignProg <- renderPrint({ 
    if (!is.null(alignementPw()))
      cat(alignementPw()$prog)
    })
  
  #------------Plots-------------#
  
  output$qualityHist <- renderPlot({ 
    if(!is.null(inputdata()))
      plotHist.fun(inputdata()$Q,
                   paste(sampleName(),": Quality Histogram",sep=" "),
                   "Phred Qualities",
                   baseTrim()$t5,
                   length(inputdata()$seq)-baseTrim()$t3)
  })

  output$cumQualityHist <- renderPlot({
    if(!is.null(inputdata()))
      plotHist.fun(cumsum(inputdata()$Q),
                   paste(sampleName(),": Cumulative Quality Histogram",sep=" "),
                   "Cumulative Phred Qualities",
                   baseTrim()$t5,
                   length(inputdata()$seq)-baseTrim()$t3)
  })
  
  output$chrom.plot <- renderPlot(
    
    chromatogram.plot(inputdata(), baseTrim()$t5, baseTrim()$t3,
                      showOrigin = input$showOrgBases,
                      width=chromWidth(), height=1, cex.mtext=1, cex.base=1, ylim=3, 
                      filename=NULL, showtrim=input$showRmBases, showhets=TRUE), height=reactive(h())
    )
  
  #---------DLHandlers---------#
  
  output$dlfasta<- downloadHandler(
    filename = function() {
      paste0(input$seq_fname, '.fasta')
    },
    content = function(file) {
      writeFasta(trimmedSequence(), input$seq_header, file, "w", nbchar = input$seq_bpl)
    }
  )
  
  output$multidlfasta<- downloadHandler(
    filename = function() {
      paste0(input$multiseq_fname, '.fasta')
    },
    content = function(file) {
      writeFasta(kkk(), names(kkk()), file, "w", nbchar = input$multiseq_bpl)
    }
  )
  
  output$dlchrom<- downloadHandler(
    filename = function() {
      paste0(input$chrom_fname, '.pdf')
    },
    content = function(file) {
      chromatogram.plot(inputdata(), baseTrim()$t5, baseTrim()$t3,
                        showOrigin = input$showOrgBases,
                        width=chromWidth(), height=1, filename=file,
                        showtrim=input$showRmBases, showhets=TRUE)
    }
  )
  
  output$dlreport <- downloadHandler(
    filename = function() {
      paste(input$report_fname, sep = '.', switch(
        input$reportFormat, PDF = 'pdf', HTML = 'html', Word = 'docx'
      ))
    },
    content = function(file) {
      src <- normalizePath('report.Rmd')
      
      # temporarily switch to the temp dir, in case you do not have write
      # permission to the current working directory
      owd <- setwd(tempdir())
      on.exit(setwd(owd))
      file.copy(src, 'report.Rmd')

      out <- render('report.Rmd', switch(
        input$reportFormat,
        PDF = pdf_document(), HTML = html_document(), Word = word_document()
      ))
      file.rename(out, file)
    }
  )
  
  #-----------Observe/Update--------------#
  
  observe({
    
    if(!is.null(inputdata()))
      updateTabsetPanel(session, "maintabset", selected = "Chromatogram")
    else
      updateTabsetPanel(session, "maintabset", selected = "Instructions") 
  })
  
  observe({
    
    if(!input$recall)
      updateCheckboxInput(session,'showOrgBases', 'Show Original Bases', FALSE)
  })
  
  observe({
    
    if(is.null(inputdata())){
      
      updateNumericInput(session,"trim5",  value = "0")
      
      updateNumericInput(session,"trim3",  value = "0")
      
      updateCheckboxInput(session,'autoTrim', 'Auto Trim', FALSE)
      
      updateCheckboxInput(session,'recall', 'Recall N peaks', FALSE)
      
      updateCheckboxInput(session,'showRmBases', 'Show Removed Bases',FALSE)
      
      updateCheckboxInput(session,'showOrgBases', 'Show Original Bases', FALSE)
      
      updateNumericInput(session,"x",  value = "100")
      
      updateRadioButtons(session,"align_type", selected = "global")
      
      updateSelectInput(session,"blastProg", selected = "Nucleotide, BLASTn")
    }
  })
  
  #-----------Multiple--------------#
  
  kkk <- reactive({ 
    if(!is.null(input$multifiles)) {
      kk<-multiProcess(input$multifiles,input$multirecall,
                       input$multiautoTrim,input$multitrim5,input$multitrim3)
      
      closeAlert(session, alertId="alert_kk")
      if(!is.null(kk$unprocessed)){
          createAlert(session, anchorId = "unprocessed",
                      alertId="alert_kk",
                      content = paste("the following files couldn't be analyzed:",
                                      paste(kk$unprocessed,collapse = ', '),sep="\t"),
                      style = "warning",
                      append=FALSE,
                      #block = FALSE,
                      dismiss = FALSE
          )
    }
    return(kk$processed)
    }
  })
  
  
  
})



