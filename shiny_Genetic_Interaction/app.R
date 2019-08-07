library(shiny)
library(tidyverse)
library(stringr)
library(grid)

load("./appdata/VST_normalized_filtered_shiny_new.RData")

#replace Ho by control
VST_normalized_filtered$LibraryName[VST_normalized_filtered$LibraryName == "HO"] <- "HO(control)"
#VST_normalized_filtered<-rename(VST_normalized_filtered, `Query Mutant`=LibraryName)
colnames(VST_normalized_filtered)[colnames(VST_normalized_filtered)=="LibraryName"] <- "Query Mutant"
VST_normalized_filtered <- VST_normalized_filtered %>%
  mutate(GrowStage = factor(ifelse(GrowStage == "G", "proliferation", "quiescence")))

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  if (numPlots==1) {
    print(plots[[1]])
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

#assign color to each library
Query <- c("HO(control)","RIM15","TOR1","PHO85")
col <- c("#E5C494","#E78AC3","#66C2A5","#8DA0CB")

col_tb<-tbl_df(cbind(Query,col))

plotgoi <- function( GeneID="HRD1 HRD3 TOR1", 
                     Library="RIM15 PHO85 TOR1 HO(control)",
                     Growth = "proliferation quiescence"){
  GeneID <- str_split(GeneID,"\\s+")[[1]]
  GeneID <- GeneID[grepl("\\S",GeneID)]
  GeneID <- toupper(GeneID)
  
  Library <- str_split(Library,"\\s+")[[1]]
  Library <- Library[grepl("\\S",Library)]
  Library <- toupper(Library)
  
  Growth <- str_split(Growth,"\\s+")[[1]]
  Growth <- Growth[grepl("\\S",Growth)]
  Growth <- tolower(Growth)

lm_plot<-VST_normalized_filtered %>%
  filter(`Query Mutant` %in% c("HO(control)",Library[1]) & Common == GeneID[1] & GrowStage == Growth[1]) %>%
  filter(!(Common == GeneID[1] & `Query Mutant` == GeneID[1])) %>%
  ggplot(aes(x=SamplingTime, y=Rela2WT, group=`Query Mutant`, color=`Query Mutant`)) +
  geom_point(aes(shape = `Query Mutant`))+
  geom_smooth(method = "lm", aes(fill = `Query Mutant`))+
  ggtitle(paste0(GeneID[1], " in ", Growth[1]))+  
  xlab("Hours post inoculation")+
  ylab("Relative frequency to query mutant")+
  scale_color_manual(values=col_tb$col[col_tb$Query %in% c("HO(control)", Library[1])])+
  scale_fill_manual(values=col_tb$col[col_tb$Query %in% c("HO(control)", Library[1])])+
  facet_grid(~nutrient, scales = "free")+
  theme_classic()
  
  #caculate the confidence interval for ATG7 slopes
  dt <- data.frame()
  for (n in sort(unique(VST_normalized_filtered$nutrient))){
    for (i in unique(VST_normalized_filtered$`Query Mutant`)){
      for (m in unique(VST_normalized_filtered$GrowStage)){
        tmp <- VST_normalized_filtered %>%
          filter(`Query Mutant` == i, Common == GeneID[1] & nutrient == n & GrowStage == m) %>%
          group_by(nutrient, `Query Mutant`, Strain, GrowStage) %>%
          lm(Rela2WT ~ SamplingTime, .) 
      
      coef <- summary(tmp)
      #caculate the 95% confidence interval
      coef <- coef$coefficients[2,1] 
      low_ci <- confint(tmp, level=0.95)[2,1]
      high_ci <- confint(tmp, level=0.95)[2,2]
      
      df <- data.frame(n, i, m, coef, low_ci, high_ci)
      dt <- rbind(dt,df)
    }
  }
}
  
  colnames(dt)<-c("nutrient","Query Mutant","GrowStage","fitness","low_ci","high_ci")
  library(ggsignif)
  ##dot plot of slopes with 95% confidenc interval
  dt$`Query Mutant` <- factor(dt$`Query Mutant`, levels = c("HO(control)","TOR1","RIM15","PHO85"))
  
  dot_ci<-ggplot(dt %>% dplyr::filter(GrowStage == Growth[1] & `Query Mutant` %in% c("HO(control)",Library[1])), aes(x=`Query Mutant`, y=fitness, fill=`Query Mutant`)) +
    geom_hline(yintercept=0, col="grey", linetype = "dashed") +
    geom_point(aes(col = `Query Mutant`, size=1.2))+
    geom_errorbar(aes(ymin = low_ci, ymax = high_ci, col = `Query Mutant`), width = 0.2)+
    scale_color_manual(values=col_tb$col[col_tb$Query %in% c("HO(control)", Library[1])])+
    scale_fill_manual(values=col_tb$col[col_tb$Query %in% c("HO(control)", Library[1])])+
    facet_grid(~nutrient)+
    xlab("Query mutant")+
    ylab("Relative fitness/survival rate")+
    theme_classic()
  
#  return(lm_plot)
  multiplot(lm_plot, dot_ci)
}

ui <- fluidPage(
  titlePanel("Fitness/Survival of single and double mutant for highthrouput genetic interaction"),
# change to multi-page navbar, so can put experimental details in here
  fluidRow(p("Supplementary ShinyApp to plot data associated with the experiment 
              from Siyu, Anastasia, Brandt, and Gresham 2019."),
           p("Technical issues, please email: dg107@nyu.edu.")),
  hr(),
  sidebarLayout(
    sidebarPanel(
      #helpText("Enter in the gene features you're interested in plotting."),
      #hr(),
      textInput(inputId="GeneID", value = "hrd1",
        label=h4("Array mutant (common name)"),
        placeholder="hrd1 hrd3 msn2 GUA1 IMD2 IMD3 IMD4"
        ),
      radioButtons(inputId="query_gene",
        label=h4("Signaling kinases (query mutant)"),
        choices = list("RIM15"=1, "TOR1"= 2, "PHO85" = 3),
        #default select
        selected = 1),
      radioButtons(inputId="Grow_Stage",
                   label=h4("Look at genetic interactions under:"),
                   choices = list("proliferation"=1, "quiescence"= 2),
                   #default select
                   selected = 1),
      width=3
      ),
    mainPanel(
      fluidRow(plotOutput("plot", width  = "700px",height = "500px"),
      fluidRow(textOutput("genez_report"))
      )
      )
    )
  )

server <- function(input, output) {
  output$plot <- renderPlot({
    plotgoi(input$GeneID,
      Library=switch(input$query_gene,`1`="RIM15",`2`="TOR1",`3`="PHO85"),
      Growth=switch(input$Grow_Stage, `1`="proliferation", `2`="quiescence")
      )
    })
}

shinyApp(ui = ui, server = server)

#you will be able to reach the html with pasting the link into your browser: 
#http://shiny.bio.nyu.edu/users/ss6025/shiny_Genetic_Interaction/