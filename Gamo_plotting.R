###Basic plotting of Gamo Results for SAA Poster

#Necessary Functions

plottR=function(x,y,
                percent=TRUE,
                spacing=5,
                min_LYCO=0,
                min_PHYT=0,
                min_POL=0,
                scale_exclude=FALSE,
                stacked=FALSE,
                set_windows=FALSE,
                window_size=100,
                stem_labels=NA,
                exclude_columns=NA,
                exclude_size=25,
                #exclude=NA, #Let's lose the "exclude" command, include "independent" plotting?
                plot_boxes=TRUE,
                box_width=5,
                pl_colors=0,
                sample_lines=FALSE,
                by_depth=TRUE,
                plot_title=NA,
                point_limit=1,
                log_concentrations=FALSE,
                manual_range=FALSE,
                min_y=0,
                max_y=100,
                label_buffer=75,
                label_anchor=5,
                plot_zones=FALSE,
                clust=NA,
                n_groups=3,
                zone_labels=NA
                #aspct_ratio=
)
{ 
  if(class(stem_labels)=="character"){
    taxa=stem_labels
  } else {
    taxa=colnames(x)
  }
  
  
  if(percent==TRUE){ #Path for calculating percents rather than raw values ###CALCULATING PERCENTS IS STUPID
    
    #Scan for values to exclude
    LYCO_n=y$LYCO_n
    x=x[LYCO_n>min_LYCO,]
    
    if(by_depth==FALSE){
      all_sample_depths=y$cal_BP
      cut_sample_depths=all_sample_depths[LYCO_n<=min_LYCO]
      use_sample_depths=all_sample_depths[LYCO_n>min_LYCO]
    } else {
      all_sample_depths=y$depth
      cut_sample_depths=all_sample_depths[LYCO_n<=min_LYCO]
      use_sample_depths=all_sample_depths[LYCO_n>min_LYCO]
    }
    
    
    #Need maximum widths of each type
    #taxa=colnames(x) - testing something here
    x[is.na(x)] = 0
    x_max=apply(x,2,max)
    
    if(scale_exclude==TRUE){ #We're dropping in spaces for the appended data
      x_max[(length(x_max)-(length(exclude_columns)-1)):length(x_max)]=exclude_size
    }
    
    windows=ceiling(x_max) #Each morphotype gets a "window" set just above its max value.
    
    anchors=vector("numeric",length=length(windows))
    for(i in 1:length(anchors)){
      if(i==1){
        anchors[i]=0
      } else {
        anchors[i]=sum(windows[1:(i-1)])+(spacing*(i-1))
      }
    }
    
    plot_width=sum(windows)+(spacing*(length(windows)-1)) #Total width for the plot in percent.
    
    if(by_depth==TRUE){ #Fork here to plot by date
      y_var=y$depth[LYCO_n>min_LYCO]
    } else {
      y_var=y$cal_BP[LYCO_n>min_LYCO] #First crack at adding dates... 3.6.2024-CAK
    } #this is where you would add dates if applicable
    
    #Break here to add if statement to manage y axis plotting.
    if(manual_range==FALSE){
      plot(0,0,xlim=c(0-(plot_width*0.05),plot_width),ylim=c(max(y_var)*-1,label_buffer),pch=NA,ann=FALSE,axes=FALSE) 
      #Do we need divergent paths for plotting by date or depth?
    } else {
      plot(0,0,xlim=c(0-(plot_width*0.05),plot_width),ylim=c(max_y*-1,(min_y+label_buffer)),pch=NA,ann=FALSE,axes=FALSE)
    }
    ### Plotting uses manual manipulation of the upper part of the y axis to create space for labels. This needs to be managed as a %. 
    
    if(is.character(plot_title)==TRUE){
      title(main=plot_title)
    }
    
    #if(plot_zones==TRUE){ #adjusting anchors and taxa list for zone title
    #  anchors=c(0-(plot_width*0.049),anchors)
    #  taxa = c("zones",taxa)
    #}
    
    text(anchors,label_anchor,labels=taxa,srt=45,adj=c(0,0),cex=0.75)
    #This is where we break it for mannual_range
    if(manual_range==FALSE){
      axis(2,seq(max(y_var)*-1,0,5),las=1,labels=seq(max(y_var),0,-5),cex.axis=0.75)
    } else {
      axis(2,seq(max_y*-1,min_y*-1,max_y*0.05),las=1,labels=seq(max_y,min_y,max_y*-0.05),cex.axis=0.75)
    }
    
    
    for(i in 1:ncol(x)){ #Loop for plotting data by columns, note that everything that's submitted gets plotted, "exclude" only refers to % calculation
      top_pts=c(rep(anchors[i],2))
      btm_pts=c(0,max(y_var)*-1)
      #print(c(top_pts,btm_pts)) #Testing if this was working, if there's problems, try turning this on.
      #Forking here to manage drawing stems
      if(manual_range==TRUE){
        btm_pts=c(min_y,max_y)*-1
      }
      lines(top_pts,btm_pts)
      
      if(plot_boxes==TRUE){ #Splitting path here to focus on plotting boxes, can incorporate polygons later...
        if(scale_exclude==TRUE){
          if(i<ncol(x)-(length(exclude_columns)-1)){
            for(ii in 1:nrow(x)){
              if(x[ii,i]>point_limit){
                box_x=c(0,0,x[ii,i],x[ii,i])+anchors[i]
                box_y=c(y_var[ii],y_var[ii]+box_width,y_var[ii]+box_width,y_var[ii])*-1
                #polygon(box_x,box_y,col="black")
                polygon(box_x,box_y,col=pl_colors[i])
              } else {
                pch_val=21
                pch_val[x[ii,i]==0]=NA
                points(x[ii,i]+anchors[i],(y_var[ii]+(box_width*0.5))*-1,pch=pch_val,bg="black",cex=0.5)
              }
              
            }
          } else {
            box_mx=max(x[,i],na.rm=TRUE)
            box_pr=x[,i]/box_mx
            box_points=box_pr*exclude_size
            for(ii in 1:nrow(x)){
              box_x=c(0,0,box_points[ii],box_points[ii])+anchors[i]
              box_y=c(y_var[ii],y_var[ii]+box_width,y_var[ii]+box_width,y_var[ii])*-1
              #polygon(box_x,box_y,col="black")
              polygon(box_x,box_y,col=pl_colors[i])
            }
          }
        } else {
          for(ii in 1:nrow(x)){
            if(x[ii,i]>point_limit){
              box_x=c(0,0,x[ii,i],x[ii,i])+anchors[i]
              box_y=c(y_var[ii],y_var[ii]+box_width,y_var[ii]+box_width,y_var[ii])*-1
              #polygon(box_x,box_y,col="black")
              polygon(box_x,box_y,col=pl_colors[i])
            } else {
              pch_val=21
              pch_val[x[ii,i]==0]=NA
              points(x[ii,i]+anchors[i],(y_var[ii]+(box_width*0.5))*-1,pch=pch_val,bg="black",cex=0.5)
            }
            
          }
        }
        
        
        
      } else { #This is where we do polygons
        
        if(scale_exclude==TRUE){
          if(i<ncol(x)-(length(exclude_columns)-1)){
            
            x_points = c(x[,i],rep(0,length(x[,i])))+anchors[i]
            y_points = c(y_var,rev(y_var))*-1
            polygon(x_points,y_points, col=pl_colors[i])
            
            for(ii in 1:nrow(x)){
              if(x[ii,i]>point_limit){
                if(sample_lines==TRUE){
                  lines(c(0,x[ii,i])+anchors[i],c(y_var[ii],y_var[ii])*-1)
                }
                #box_x=c(0,0,x[ii,i],x[ii,i])+anchors[i]
                #box_y=c(y_var[ii],y_var[ii]+box_width,y_var[ii]+box_width,y_var[ii])*-1
                #polygon(box_x,box_y,col="black")
                #polygon(box_x,box_y,col=pl_colors[i])
              } else {
                pch_val=21
                pch_val[x[ii,i]==0]=NA
                points(x[ii,i]+anchors[i],(y_var[ii]+(box_width*0.5))*-1,pch=pch_val,bg="black",cex=0.5)
              }
              
            }
          } else {
            box_mx=max(x[,i],na.rm=TRUE)
            box_pr=x[,i]/box_mx
            box_points=box_pr*exclude_size
            
            x_points = c(box_points,rep(0,length(box_points)))+anchors[i]
            y_points = c(y_var,rev(y_var))*-1
            #x_points=c(box_points,rep(0,length(box_points)))
            #y_points=c(y_var,rev(y_var))
            
            polygon(x_points,y_points,col=pl_colors[i])
            
            for(ii in 1:nrow(x)){
              #box_x=c(0,0,box_points[ii],box_points[ii])+anchors[i]
              #box_y=c(y_var[ii],y_var[ii]+box_width,y_var[ii]+box_width,y_var[ii])*-1
              #polygon(box_x,box_y,col="black")
              #polygon(box_x,box_y,col=pl_colors[i])
            }
          }
        } else {
          
          x_points = c(x[,i],rep(0,length(x[,i])))+anchors[i]
          y_points = c(y_var,rev(y_var))*-1
          polygon(x_points,y_points, col=pl_colors[i])
          
          for(ii in 1:nrow(x)){
            if(x[ii,i]>point_limit){
              if(sample_lines==TRUE){
                lines(c(0,x[ii,i])+anchors[i],c(y_var[ii],y_var[ii])*-1)
              }
              #box_x=c(0,0,x[ii,i],x[ii,i])+anchors[i]
              #box_y=c(y_var[ii],y_var[ii]+box_width,y_var[ii]+box_width,y_var[ii])*-1
              #polygon(box_x,box_y,col="black")
              #polygon(box_x,box_y,col=pl_colors[i])
            } else {
              pch_val=21
              pch_val[x[ii,i]==0]=NA
              points(x[ii,i]+anchors[i],(y_var[ii]+(box_width*0.5))*-1,pch=pch_val,bg="black",cex=0.5)
            }
            
          }
        }
        
        
        
      }
      
      if(i<ncol(x)-(length(exclude_columns)-1)){
        axis(1,at=c(anchors[i],anchors[i]+windows[i]),labels=c(NA,windows[i]),cex.axis=0.75)
      } else {
        axis(1,at=c(anchors[i],anchors[i]+windows[i]),labels=c(NA,ceiling(max(x[,i],na.rm=TRUE))),cex.axis=0.75)
      }
      
      #HOLD POINT FOR ZONES
      
    } #Closing column plotting loop.
    
    if(plot_zones==TRUE){ ###Adding zones to plots, supply a chclust object and n groups.
      clusters = cutree(clust, k = n_groups)
      text_anchors=vector("numeric",length = length(n_groups))
      text(0-(plot_width*0.04),label_anchor,labels="zones",srt=45,adj=c(0,0),cex=0.75) #adding zone label
      if(class(zone_labels)=="character"){
        
      } 
      else {
        zone_labels=as.character(c(1:n_groups))
      }
      for(j in 2:n_groups){ #Now we find lines between the major groups. Can append this to the cluster building.
        xj = max(y_var[clusters == j-1])*-1
        yj = min(y_var[clusters == j])*-1
        b = (xj + (yj - xj)/2)
        lines(c(0-(plot_width*0.05),plot_width), c(b, b), lty = 3)
        text_anchors[j]=b
      }
      for(k in 1:n_groups){
        if(k!=n_groups){
          text(0-(plot_width*0.03),text_anchors[k+1]+((text_anchors[k]-text_anchors[k+1])/2),zone_labels[k],cex=0.6,adj=c(1,0))
        } else {
          text(0-(plot_width*0.03),text_anchors[k]-((text_anchors[k]+max_y)/2),zone_labels[k],cex=0.6,adj=c(1,0))
        }
        
      }
    } else {}
    
    points(-rep(spacing,length(use_sample_depths)),(use_sample_depths*-1)-(box_width*0.5),pch=21,cex=0.5)
    points(-rep(spacing,length(cut_sample_depths)),(cut_sample_depths*-1)-(box_width*0.5),pch=4,cex=0.5)
    
  } else { #Closing percent calculation path, this is where raw value plots would begin.
    
    #if(stacked==TRUE & summarize_taxa==TRUE){ #Building out loop to handle summarizing results by taxa groups.
    
    #We need summarized results for select samples (thresholds for sample size and lycopodium?), then we plot as points or boxes?
    
    #} #Closing out summarizing by group plots.
    
    if(stacked==TRUE & set_windows==TRUE){ #Opening section for plotting stacked results with even windows.
      
      if(log_concentrations==TRUE){
        conc_finder=grep("conc",colnames(x))
        for(i in 1:length(conc_finder)){
          x[,conc_finder[i]]=log(x[,conc_finder[i]])
        }
      }
      
      LYCO_n=y$LYCO_n
      
      all_sample_depths=y$depth
      cut_sample_depths=all_sample_depths[LYCO_n<=min_LYCO]
      use_sample_depths=all_sample_depths[LYCO_n>min_LYCO]
      
      x=x[LYCO_n>min_LYCO,]
      y=y[LYCO_n>min_LYCO,]
      #taxa=colnames(x)
      x_max=apply(x,2,max)
      
      if(scale_exclude==TRUE){ #We're dropping in spaces for the appended data
        x_max[(length(x_max)-(length(exclude_columns)-1)):length(x_max)]=exclude_size
      }
      
      windows=rep(window_size,ncol(x)/2) #Each morphotype gets a "window" set just above its max value.
      
      anchors=vector("numeric",length=length(windows))
      for(i in 1:length(anchors)){
        if(i==1){
          anchors[i]=0
        } else {
          anchors[i]=sum(windows[1:(i-1)])+(spacing*(i-1))
        }
      }
      
      plot_width=sum(windows)+(spacing*(length(windows)-1)) #Total width for the plot in percent.
      
      if(by_depth==TRUE){ #Fork here to plot by date
        y_var=y$depth
      } else {} #this is where you would add dates if applicable
      
      
      #Break here to add if statement to manage y axis plotting.
      if(manual_range==FALSE){
        plot(0,0,xlim=c(0-(plot_width*0.05),plot_width),ylim=c(max(y_var)*-1,label_buffer),pch=NA,ann=FALSE,axes=FALSE) 
        #Do we need divergent paths for plotting by date or depth?
      } else {
        plot(0,0,xlim=c(0-(plot_width*0.05),plot_width),ylim=c(max_y*-1,(min_y+label_buffer)),pch=NA,ann=FALSE,axes=FALSE)
      }
      ### Plotting uses manual manipulation of the upper part of the y axis to create space for labels. This needs to be managed as a %. 
      
      if(is.character(plot_title)==TRUE){
        title(main=plot_title)
      }
      
      #Can we draw boxes straight out?
      
      evens=seq(2,ncol(x),2)
      odds=seq(1,ncol(x),2)
      
      plot_placement=vector("character",length=ncol(x))
      
      plot_placement[evens]="front"
      plot_placement[odds]="back"
      
      #taxa=stem_labels
      
      anchors=rep(anchors,2)
      anchors=anchors[order(anchors)]
      
      text(anchors[evens]+(window_size*0.3),5,labels=taxa[evens],adj=c(0,0),cex=0.75) ###This is broken in stacked logic - FIX
      text(anchors[odds]+(window_size*0.3),15,labels=taxa[odds],adj=c(0,0),cex=0.75)
      #Break here with manual_range
      if(manual_range==FALSE){
        axis(2,seq(max(y_var)*-1,0,5),las=1,labels=seq(max(y_var),0,-5),cex.axis=0.75)
      } else {
        axis(2,seq(max_y*-1,min_y*-1,5),las=1,labels=seq(max_y,min_y,-5),cex.axis=0.75)
      }
      
      
      for(i in 1:ncol(x)){
        #print(c(top_pts,btm_pts)) #Testing if this was working, if there's problems, try turning this on.
        if(plot_boxes==TRUE){ #Splitting path here to focus on plotting boxes, can incorporate polygons later..
          if(plot_placement[i]=="front"){
            box_mx=max(x[,i-1],na.rm=TRUE)
            box_pr=x[,i]/box_mx
            box_points=box_pr*window_size
            for(ii in 1:nrow(x)){
              
              if(x[ii,i]>point_limit){
                
                box_x=c(0,0,box_points[ii],box_points[ii])+anchors[i]
                box_y=c(y_var[ii],y_var[ii]+box_width,y_var[ii]+box_width,y_var[ii])*-1
                polygon(box_x,box_y,col="gray")
                
              } else {
                pch_val=21
                pch_val[x[ii,i]==0]=NA
                points(x[ii,i]+anchors[i],(y_var[ii]+(box_width*0.5))*-1,pch=pch_val,bg="gray",cex=0.5)
                #points(x[ii,i]+anchors[i],(y_var[ii]+(box_width*0.5))*-1,pch=21,bg="gray",cex=0.5)
              }
              
            }
            
          } else { #"BACK" set of plots.
            box_mx=max(x[,i],na.rm=TRUE)
            box_pr=x[,i]/box_mx
            box_points=box_pr*window_size
            for(ii in 1:nrow(x)){
              
              if(x[ii,i]>point_limit){
                
                box_x=c(0,0,box_points[ii],box_points[ii])+anchors[i]
                box_y=c(y_var[ii],y_var[ii]+box_width,y_var[ii]+box_width,y_var[ii])*-1
                polygon(box_x,box_y,col="black")
                
              } else {
                pch_val=21
                pch_val[x[ii,i]==0]=NA
                points(x[ii,i]+anchors[i],(y_var[ii]+(box_width*0.5))*-1,pch=pch_val,bg="black",cex=0.5)
                #points(x[ii,i]+anchors[i],(y_var[ii]+(box_width*0.5))*-1,pch=21,bg="black",cex=0.5)
              }
              
            }
            
            labelbox_x=c(0,0,25,25)+anchors[i]
            labelbox_y1=c(5,5+box_width,5+box_width,5)
            labelbox_y2=c(15,15+box_width,15+box_width,15)
            
            polygon(labelbox_x,labelbox_y1,col="gray")
            polygon(labelbox_x,labelbox_y2,col="black")
            
            top_pts=c(rep(anchors[i],2))
            btm_pts=c(0,max(y_var)*-1)
            #Forking here to manage drawing stems
            if(manual_range==TRUE){
              btm_pts=c(min_y,max_y)*-1
            }
            lines(top_pts,btm_pts)
            axis(1,at=c(anchors[i],anchors[i]+window_size),labels=c(NA,ceiling(max(x[,i],na.rm=TRUE))),cex.axis=0.75)
            
          }
          
          #if(i<ncol(x)-(length(exclude_columns)-1)){
          #  axis(1,at=c(anchors[i],anchors[i]+windows[i]),labels=c(0,windows[i]),cex.axis=0.5)
          #} else {
          #  axis(1,at=c(anchors[i],anchors[i]+windows[i]),labels=c(0,ceiling(max(x[,i],na.rm=TRUE))),cex.axis=0.5)
          #}
          
        } else {} #This is where we would do polygons #Closing if/else for box plotting
        
      } #Closing Column loop
      
      points(-rep(spacing,length(use_sample_depths)),(use_sample_depths*-1)-(box_width*0.5),pch=21,cex=0.5)
      points(-rep(spacing,length(cut_sample_depths)),(cut_sample_depths*-1)-(box_width*0.5),pch=4,cex=0.5)
      
    } else {} #Closing loop for Stacked and Set Window == TRUE/FALSE
    
  } #Closing else loop for percent==TRUE
  
} #Closing out function

#Logical controls on plotting figs

save_figs = TRUE

#Data Read

GAMO_all <- read.csv("GAMO_all.csv", header = TRUE) #Revised files

GAMO_LOI <- read.csv("GAMO_LOI.csv", header = TRUE)

GAMO_AMS <- read.csv("GAMO_AMS.csv", header = TRUE)

#Parse LOI results with names

LOI_cores = as.list(unique(GAMO_LOI$Core.ID))

LOI_cores = LOI_cores[order(unique(GAMO_LOI$Core.ID))]

#Isolate core details

GAMO_details = as.data.frame(t(GAMO_all[GAMO_all$CODE == "00_CORE", 6:ncol(GAMO_all)]))

#Add new headers

colnames(GAMO_details) = GAMO_all[1:6,3]


#Isolate counts

GAMO_counts = as.data.frame(t(GAMO_all[GAMO_all$CODE != "00_CORE" & GAMO_all$CODE != "01_PTERID", 6:ncol(GAMO_all)]))

#Add new headers

colnames(GAMO_counts) = GAMO_all[GAMO_all$CODE != "00_CORE" & GAMO_all$CODE != "01_PTERID", 3]

#Get rows for different cores

cores = list("CHA", "KAO", "CHO")

core_finder = lapply(cores, function(x){
  grep(x, row.names(GAMO_counts))
})

names(core_finder) = cores

#Sums and Concentrations

GAMO_sums = lapply(cores, function(x){
  apply(GAMO_counts[core_finder[[x]],],1,sum)
})

names(GAMO_sums) = cores

GAMO_conc = lapply(cores, function(x){
  (GAMO_sums[[x]]/GAMO_details$LYCO_n[core_finder[[x]]])*((GAMO_details$LYCO_mean[core_finder[[x]]]*2)/(GAMO_details$weight_g[core_finder[[x]]]))
})

names(GAMO_conc) = cores

conc_limit = lapply(cores, function(x){
  ceiling(max(GAMO_conc[[x]]))
})

names(conc_limit) = cores

#Simple plots

if(save_figs == TRUE){
  setEPS()
  tiff("Figure-2_all-cores.tiff", height = 3000, width = 2500, res = 300)
}

spacer = 50

plot(0, 0, axes = FALSE, ann = FALSE, pch = NA, xlim = c(-30,400), ylim = c(-180,0))

axis(2, seq(0,-160,-20))

for(i in 1:length(cores)){
  lines(-1*(GAMO_LOI$l.o.i...TOM.[GAMO_LOI$Core.ID == LOI_cores[[i]]])+((i-1)*(100+spacer)), -1*GAMO_LOI$sample.depth_from[GAMO_LOI$Core.ID == LOI_cores[[i]]], lty = 2)
  axis(3,at = c(-30,0)+((i-1)*(100+spacer)), labels = c(50,0), cex.axis = 0.8)
}

for(i in 1:length(cores)){
  for(j in 1:length(core_finder[[i]])){
    polygon(c(rep(100*(GAMO_conc[[i]][j]/conc_limit[[i]])+((i-1)*(100+spacer)),2),rep(0+((i-1)*(100+spacer)),2)),
            -1*c(GAMO_details$depth[core_finder[[i]]][j]+2,GAMO_details$depth[core_finder[[i]]][j],GAMO_details$depth[core_finder[[i]]][j],GAMO_details$depth[core_finder[[i]]][j]+2),
            col = "darkgrey")
    
  }
  
  for(j in 1:length(core_finder[[i]])){
    polygon(c(rep(((100*(GAMO_conc[[i]][j]/conc_limit[[i]]))*(GAMO_details$IND[core_finder[[i]][j]]/GAMO_sums[[i]][j]))+((i-1)*(100+spacer)),2),rep(0+((i-1)*(100+spacer)),2)),
            -1*c(GAMO_details$depth[core_finder[[i]]][j]+2,GAMO_details$depth[core_finder[[i]]][j],GAMO_details$depth[core_finder[[i]]][j],GAMO_details$depth[core_finder[[i]]][j]+2),
            col = "black")
  }
  axis(1, at = c(0,100)+((i-1)*(100+spacer)), labels = c(0, conc_limit[[i]]), cex.axis = 0.8)
}

for(i in 1:length(cores)){
  points(rep(0+((i-1)*(100+spacer)), length(core_finder[[i]])), -1*(GAMO_details$depth[core_finder[[i]]]+1), pch = 21, bg = "goldenrod")
  text(50+((i-1)*(100+spacer)), 3, paste(cores[[i]]))
}

for(i in 1:length(cores)){
  points(rep(0+((i-1)*(100+spacer)), length(GAMO_AMS$Core[GAMO_AMS$Core==cores[[i]]])), -1*GAMO_AMS$depth.top..cm.[GAMO_AMS$Core==cores[[i]]], pch = 24, bg = "forestgreen")
  text(rep(0+((i-1)*(100+spacer)+25), length(GAMO_AMS$Core[GAMO_AMS$Core==cores[[i]]])), -1*GAMO_AMS$depth.top..cm.[GAMO_AMS$Core==cores[[i]]],
       paste0(GAMO_AMS$cal.BP.2sig[GAMO_AMS$Core==cores[[i]]]), cex = 0.8)
}

legend(310, -150, c("LOI%", "conc/g", "IND", "Pollen", "14C"), lty = c(2, NA, NA, NA, NA), pch = c(NA, 22, 22, 21, 24), pt.bg = c(NA, "darkgrey", "black", "goldenrod", "forestgreen"))

if(save_figs == TRUE){
  dev.off()
}

#Pollen Diagram by summed categories

GAMO_taxa = GAMO_all[GAMO_all$CLASS != "Core",c(1:5)]
GAMO_groups = unique(GAMO_taxa$CODE)
group_names = unique(GAMO_taxa$ECOLOGY)

CHO_counts = GAMO_counts[core_finder[["CHO"]],GAMO_taxa$CODE != "00_CORE"]
CHO_detail = GAMO_details[core_finder[["CHO"]],GAMO_all$CODE == "00_CORE"]

CHO_groups = matrix(nrow = nrow(CHO_counts), ncol = length(GAMO_groups))
row.names(CHO_groups) = row.names(CHO_counts)
colnames(CHO_groups) = GAMO_groups


for(i in 1:length(GAMO_groups)){
  n = CHO_counts[,GAMO_taxa$CODE == GAMO_groups[i]]
  if(is.vector(n) == TRUE){
    CHO_groups[,i] = n
  } else {
    CHO_groups[,i] = apply(n, 1, sum)
  }
}

#Percents

CHO_pct = (CHO_groups/GAMO_sums$CHO)*100

#Set colors

Ecol_colors = c("darkgreen", "lightblue", "goldenrod", "darkorange", "forestgreen", "gold", "black")

#PlottR

if(save_figs == TRUE){
  setEPS()
  tiff("Figure-4_CHO-.tiff", height = 2000, width = 3000, res = 300)
}

par(mar = c(5, 4, 4, 5) + 0.1)

plottR(CHO_pct, CHO_detail, point_limit = 3, pl_colors = Ecol_colors, stem_labels = group_names, plot_title = "CHO Pollen %")

if(save_figs == TRUE){
  dev.off()
}

par(mar = c(5, 4, 4, 2) + 0.1)




