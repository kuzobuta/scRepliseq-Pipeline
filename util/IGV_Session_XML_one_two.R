####IGV session XML####

##Variables##
#files: vector of each file path
#name: vector of each file name
#############

##blank file##
#blankfile="mm9.blank.bed"

#For binary data
#max=1
#min=-1

#For log2 med data
#max=0.75
#min=-0.75

##############################################################################################
IGV_XML_export = function(files,out_name,out_dir_m,out_dir,genome,max,min){
library("XML")
XML_path_m=paste(out_dir_m,"/",out_name,sep="")
XML_path=paste(out_dir,"/",out_name,sep="")
#0.Session
Session=xmlNode("Session",attrs=c(genome=genome,
                            hasGeneTrack="true",
                            hasSequenceTrack="true",
                            locus="All",path=XML_path,
                            version="8"))

#1.Resources variable: file (abs path), name (name only)
Resources=xmlNode("Resources")
##Generate it for each files##
for (i in 1:length(files)){

file=files[i]

R_file=addAttributes(xmlNode("Resource"),path=file)
Resources=addChildren(Resources,R_file)
}
##############################

#2. Panel1 (data file), variable: file, name
Panel1=addAttributes(xmlNode("Panel"),height="509",name="DataPanel",width="1133")
##Generate it for each files##
for (i in 1:length(files)){

file=files[i]
name=basename(file)
#name=names[i]

Track1=addAttributes(xmlNode("Track"),
                    altColor="204,204,0",
                    autoScale="false",
                    clazz="org.broad.igv.track.DataSourceTrack",
                    color="0,0,178",
                    colorScale=paste0("ContinuousColorScale;0.0;",min,";0.0;",max,";204,204,0;255,255,255;0,0,178"),
                    displayMode="COLLAPSED",featureVisibilityWindow="-1",fontSize="10",height="5",
                    id=file,
                    name=name,
                    renderer="HEATMAP",
                    sortable="true",
                    visible="true")
DataRange1=addAttributes(xmlNode("DataRange"),
                        baseline="0.0",
                        drawBaseline="true",
                        flipAxis="false",
                        maximum=max,
                        minimum=min,
                        type="LINEAR")
Track1=addChildren(Track1,DataRange1)
Panel1=addChildren(Panel1,Track1)
}

##############################

#3. Panel2 (refseq file)
Panel2=addAttributes(xmlNode("Panel"),height="72",name="FeaturePanel",width="1133")
Track2=addAttributes(xmlNode("Track"),
                     altColor="0,0,178",
                     autoScale="false",
                     clazz="org.broad.igv.track.SequenceTrack",
                     color="0,0,178",
                     displayMode="COLLAPSED",
                     featureVisibilityWindow="-1",
                     fontSize="10",
                     id="Reference sequence",
                     name="Reference sequence",
                     sortable="false",
                     visible="true")
Track3=addAttributes(xmlNode("Track"),
                     altColor="0,0,178",
                     autoScale="false",
                     clazz="org.broad.igv.track.FeatureTrack",
                     color="0,0,178",colorScale="ContinuousColorScale;0.0;226.0;255,255,255;0,0,178",
                     displayMode="COLLAPSED",featureVisibilityWindow="-1",fontSize="10",
                     height="35",id=paste0(genome,"_genes"),name="Refseq genes",renderer="BASIC_FEATURE",sortable="false",
                     visible="true",windowFunction="count")
DataRange1=addAttributes(xmlNode("DataRange"),
                         baseline="0.0",
                         drawBaseline="true",
                         flipAxis="false",
                         maximum="226.0",
                         minimum="0.0",type="LINEAR")
Track3=addChildren(Track3,DataRange1)
Panel2=addChildren(Panel2,Track2,Track3)

#4. PanelLayout
PanelLayout=addAttributes(xmlNode("PanelLayout"),dividerFractions="0.8690476190476191")

#Merge the Children
IGV=addChildren(Session,Resources,Panel1,Panel2,PanelLayout)

#Save
saveXML(IGV,file=XML_path_m,prefix='<?xml version="1.0" encoding="UTF-8" standalone="no"?>\n')

#Print the log
print(paste("Save as ",XML_path,sep=""))
print(paste(length(files)," files",sep=""))
}

##############################################################################################
IGV_XML_export_two_set = function(files1,files2,out_name,out_dir_m,out_dir,genome,blankfile,max1,min1,max2,min2){
  library("XML")
  XML_path_m=paste(out_dir_m,"/",out_name,sep="")
  XML_path=paste(out_dir,"/",out_name,sep="")
  #0.Session
  Session=xmlNode("Session",attrs=c(genome=genome,
                                    hasGeneTrack="true",
                                    hasSequenceTrack="true",
                                    locus="All",path=XML_path,
                                    version="8"))
  
  #1.Resources variable: file (abs path), name (name only)
  Resources=xmlNode("Resources")
  ##Generate it for each files (1st set)##
  for (i in 1:length(files1)){
    
    file=files1[i]
    
    R_file=addAttributes(xmlNode("Resource"),path=file)
    Resources=addChildren(Resources,R_file)
  }
  ##Generate it for each files (2nd set)##
  for (i in 1:length(files2)){
    
    file=files2[i]
    
    R_file=addAttributes(xmlNode("Resource"),path=file)
    Resources=addChildren(Resources,R_file)
  }
  ##For blank file
  R_file=addAttributes(xmlNode("Resource"),path=blankfile)
  Resources=addChildren(Resources,R_file)
  
  ##############################
  
  #2. Panel1 (data file), variable: file, name
  Panel1=addAttributes(xmlNode("Panel"),height="509",name="DataPanel",width="1133")
  ##Generate it for each files(1st set)##
  for (i in 1:length(files1)){
    
    file=files1[i]
    name=basename(file)
    #name=names[i]
    
    Track1=addAttributes(xmlNode("Track"),
                         altColor="204,204,0",
                         autoScale="false",
                         clazz="org.broad.igv.track.DataSourceTrack",
                         color="0,0,178",
                         colorScale=paste0("ContinuousColorScale;0.0;",min1,";0.0;",max1,";204,204,0;255,255,255;0,0,178"),
                         displayMode="COLLAPSED",featureVisibilityWindow="-1",fontSize="10",height="5",
                         id=file,
                         name=name,
                         renderer="HEATMAP",
                         sortable="true",
                         visible="true")
    DataRange1=addAttributes(xmlNode("DataRange"),
                             baseline="0.0",
                             drawBaseline="true",
                             flipAxis="false",
                             maximum=max1,
                             minimum=min1,
                             type="LINEAR")
    Track1=addChildren(Track1,DataRange1)
    Panel1=addChildren(Panel1,Track1)
  }
  ##add blank file##
  file=blankfile
  name=basename(file)
  Track1=addAttributes(xmlNode("Track"),
                       altColor="0,0,178",
                       autoScale="false",
                       clazz="org.broad.igv.track.FeatureTrack",
                       color="204,0,204",
                       colorScale="ContinuousColorScale;0.0;2.0;255,255,255;204,0,204",
                       displayMode="COLLAPSED",
                       featureVisibilityWindow="-1",
                       fontSize="10",height="20",
                       id=file,
                       name=name,
                       renderer="BASIC_FEATURE",
                       sortable="false",
                       visible="true",
                       windowFunction="count")
  DataRange1=addAttributes(xmlNode("DataRange"),
                           baseline="0.0",
                           drawBaseline="true",
                           flipAxis="false",
                           maximum="2.0",
                           minimum="0.0",
                           type="LINEAR")
  Track1=addChildren(Track1,DataRange1)
  Panel1=addChildren(Panel1,Track1)
  ##Generate it for each files(2nd set)##
  for (i in 1:length(files2)){
    
    file=files2[i]
    name=basename(file)
    #name=names[i]
    
    Track1=addAttributes(xmlNode("Track"),
                         altColor="204,204,0",
                         autoScale="false",
                         clazz="org.broad.igv.track.DataSourceTrack",
                         color="0,0,178",
                         colorScale=paste0("ContinuousColorScale;0.0;",min2,";0.0;",max2,";204,204,0;255,255,255;0,0,178"),
                         displayMode="COLLAPSED",featureVisibilityWindow="-1",fontSize="10",height="5",
                         id=file,
                         name=name,
                         renderer="HEATMAP",
                         sortable="true",
                         visible="true")
    DataRange1=addAttributes(xmlNode("DataRange"),
                             baseline="0.0",
                             drawBaseline="true",
                             flipAxis="false",
                             maximum=max2,
                             minimum=min2,
                             type="LINEAR")
    Track1=addChildren(Track1,DataRange1)
    Panel1=addChildren(Panel1,Track1)
  }
  ##############################
  
  #3. Panel2 (refseq file)
  Panel2=addAttributes(xmlNode("Panel"),height="72",name="FeaturePanel",width="1133")
  Track2=addAttributes(xmlNode("Track"),
                       altColor="0,0,178",
                       autoScale="false",
                       clazz="org.broad.igv.track.SequenceTrack",
                       color="0,0,178",
                       displayMode="COLLAPSED",
                       featureVisibilityWindow="-1",
                       fontSize="10",
                       id="Reference sequence",
                       name="Reference sequence",
                       sortable="false",
                       visible="true")
  Track3=addAttributes(xmlNode("Track"),
                       altColor="0,0,178",
                       autoScale="false",
                       clazz="org.broad.igv.track.FeatureTrack",
                       color="0,0,178",colorScale="ContinuousColorScale;0.0;226.0;255,255,255;0,0,178",
                       displayMode="COLLAPSED",featureVisibilityWindow="-1",fontSize="10",
                       height="35",id=paste0(genome,"_genes"),name="Refseq genes",renderer="BASIC_FEATURE",sortable="false",
                       visible="true",windowFunction="count")
  DataRange1=addAttributes(xmlNode("DataRange"),
                           baseline="0.0",
                           drawBaseline="true",
                           flipAxis="false",
                           maximum="226.0",
                           minimum="0.0",type="LINEAR")
  Track3=addChildren(Track3,DataRange1)
  Panel2=addChildren(Panel2,Track2,Track3)
  
  #4. PanelLayout
  PanelLayout=addAttributes(xmlNode("PanelLayout"),dividerFractions="0.8690476190476191")
  
  #Merge the Children
  IGV=addChildren(Session,Resources,Panel1,Panel2,PanelLayout)
  
  #Save
  saveXML(IGV,file=XML_path_m,prefix='<?xml version="1.0" encoding="UTF-8" standalone="no"?>\n')
  
  #Print the log
  print(paste("Save as ",XML_path,sep=""))
  print(paste(length(files1)," files (1st set)",sep=""))
  print(paste(length(files2)," files (2nd set)",sep=""))
  
}

