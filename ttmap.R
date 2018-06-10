get_clusters_by_level<-
function(ttmap_part1_hda,m1,select = row.names(ttmap_part1_hda$Dc.Dmat),ddd, e,filename = "TEST", n = 3, ad = 0, 
         bd = 0, piq = 1, dd = generate_mismatch_distance(ttmap_part1_hda = ttmap_part1_hda,select = select),
         mean_value_m1 = "N", ni = 2){
  # just writing for user
  write.table(dd, file = "distance_matrix.txt",
              sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)
  beg <- ls()
  # annotations: to write names in the graph, if the clusters are all females, I'll display it
  annot <- function(q, n = n) {
                q <- lapply(seq_len(length(q)),
                            function(i){ddd[q[[i]], n]})
                return(q)}
  qmin <- min(m1)
  q25 <- quantile(m1, 0.25)
  q5 <- quantile(m1, 0.5)
  q75 <- quantile(m1, 0.75)
  qmax <- max(m1)
  images = list(all = m1,
                low = m1[m1[] < q25],
                mid1 = m1[m1[] < q5 & m1[] >= q25],
                mid2 = m1[m1[] >= q5 & m1[] < q75],
                high = m1[m1[] >= q75]
                )
  single_linkage_algorithm<-function(image, epsilon){
      if(length(image) == 0){
        image_map <- as.matrix(0)}
      else{
        # single linkage algorithm: low_map$1 is a list of the elem of first cluster
        image_map <- mapper1(as.matrix(dd[names(image),names(image)]), 
                             e = epsilon)}
        rownames(image_map)<-names(image)
        colnames(image_map)<-names(image)
      return(image_map)}
  clusters_by_level = lapply(images, single_linkage_algorithm, epsilon=e)
  names(clusters_by_level) <- c("all_map", "low_map", "mid1_map", "mid2_map", "high_map")
  if(bd != 0){  # remove singletons (param = piq)
    # check the syntax
    clusters_by_level<-mapply(cutoff_low,
                              clusters_by_level,
                              para=rep(piq,length(clusters_by_level)),
                              text=c( "all", "low, mid1", "mid2", "high"))}
  return(clusters_by_level)}
    
calculate_size<-function(clusters_in_one_level, dd){
  if(dim(clusters_in_one_level)[1] == 1){ # do the same thing, for the global clusters
    if(clusters_in_one_level != 0){
      q <- as.list(colnames(clusters_in_one_level))
      # names(q) <- colnames(clusters_in_one_level)
    }    
    else{
      q <- c()
    }
  }
  else{
    q <- apply(clusters_in_one_level, 1, grep, pattern = 1)
    q <- as.list(q)
    q <- lapply(seq_len(length(q)), 
                    function(i){ colnames(as.matrix(dd))[q[[i]]] } )
  }
  # counting how many samples there are in each cluster to set the size of the sphere
  q1 <- q[ lapply(q, length) > 0]
  size <- lapply(q1, length)
  return(c(q1,size))
}

draw<-function(q1_and_sizes_by_level,
               ttmap_part1_hda,
               m1,
               select = row.names(ttmap_part1_hda$Dc.Dmat),
               ddd, 
               e,
               filename = "TEST", 
               n = 3, 
               ad = 0, 
               bd = 0, 
               piq = 1,
               dd = generate_mismatch_distance(ttmap_part1_hda = ttmap_part1_hda,
                                               select = select),
               mean_value_m1 = "N", 
               ni = 2,
               qmin,
               qmax,
               level = 1){
  # here I start drawing. Open3d() has be called only once, before the call to this function
  squize <- function(q1_all, m1, size_all, n = n){
    r_all <- create_colors(q1_all, m1, size_all)
    sort_r <- sort(r_all$average, index.return = TRUE)
    size_all <- size_all[sort_r$ix]
    f_all <- create_places(size_all)
    r_all$col <- r_all$col[sort_r$ix]
    q1_all <- q1_all[sort_r$ix]
    q1_all_a <- annot(q1_all, n = n)
    out1 <- list(r = r_all, s = size_all,
                 f = f_all, q1_all = q1_all, q1_all_a = q1_all_a)
    return(out1)
  }
  p <- squize(q1_all, m1, size_all, n = n)
  lans <- lapply(seq_len(length(size_all)), function(i){
    spheres3d(p$f[i], 0, 0, radius = p$s[i],
              color=rgb((p$r)$col[[i]][1], (p$r)$col[[i]][2],
                        (p$r)$col[[i]][3], alpha = 1))
    text3d(p$f[i], max(as.matrix(unlist(p$s))) + 5, 0, i)
    if(ad == 0){text3d(p$f[i], 
                       -max(as.matrix(unlist(p$s))) - 5, 0, 
                       paste(unique(p$q1_all_a[[i]]), collapse=""))
    }
  })
  m <- 0
  if(dim(low_map)[1] == 1 && low_map == 0){
    du <- max(as.matrix(unlist(size_all))) + minespace
    d_low <- du
    du <- du + minespace
    size_low <- c()
    p_low <- list()
    p_low$f <- 0
  }
  else{
    p_low <- squize(q1_low, m1, size_low, n=n)
    du <- max(as.matrix(unlist(size_all))) + 
      max(as.matrix(unlist(size_low))) + minespace
    d_low <- du
    lans <- lapply(seq_len(length(p_low$s)), function(i){
      spheres3d(p_low$f[i], du, 
                0, radius = p_low$s[i], 
                color = rgb((p_low$r)$col[[i]][1],
                            (p_low$r)$col[[i]][2], (p_low$r)$col[[i]][3],alpha = 1))
      text3d(p_low$f[i], 
             du + max(as.matrix(unlist(size_low))) + 5, 0,
             length(size_all) + i)
      if(ad == 0){text3d(p_low$f[i], 
                         (du - max(as.matrix(unlist(size_low))) - 5), 0,
                         paste(unique(p_low$q1_all_a[[i]]), collapse = ""))
      }
    })
    create_links(q = p$q1_all, q1 = p_low$q1_all, f = p$f,
                 f1 = p_low$f, m, du, 0)
    du <- du + max(as.matrix(unlist(size_low))) + minespace
  }
  if(dim(mid1_map)[1] == 1 && mid1_map == 0){
    d_mid1 <- du
    du <- du + minespace
    size_mid1 <- c()
    p_mid1 <- list()
    p_mid1$f <- 0}
  else{
    m <- max(p_low$f)
    p_mid1 <- squize(q1_mid1, m1, size_mid1, n = n)
    d_mid1 <- du + max(as.matrix(unlist(size_mid1)))
    du <- d_mid1
    l <- 10
    lans <- lapply(seq_len(length(size_mid1)), function(i){
      spheres3d(m + p_mid1$f[i], du, l,
                radius = p_mid1$s[i], color = rgb((p_mid1$r)$col[[i]][1],
                                                  (p_mid1$r)$col[[i]][2], (p_mid1$r)$col[[i]][3],alpha = 1))
      text3d(m + p_mid1$f[i], du + max(as.matrix(
        unlist(size_mid1))) + 5, l,
        length(size_all) + length(size_low) + i)
      if(ad == 0){text3d(m + p_mid1$f[i],
                         (du - max(as.matrix(unlist(size_mid1))) -5), l,
                         paste(unique(p_mid1$q1_all_a[[i]]), collapse = ""))}
    })
    create_links(p$q1_all, p_mid1$q1_all, p$f, p_mid1$f, m, du, l)
    du <- du + max(as.matrix(unlist(size_mid1)))+minespace
  }
  m <- max(p_mid1$f) + m
  if(dim(mid2_map)[1] == 1 && mid2_map == 0){
    d_mid2 <- du
    du <- du + minespace
    p_mid2 <- list()
    p_mid2$f <- 0
    size_mid2 <- c()}
  else{
    du <- du + max(as.matrix(unlist(size_mid2)))
    d_mid2 <- du
    p_mid2 <- squize(q1_mid2, m1, size_mid2, n = n)
    l <- 20
    lans <- lapply(seq_len(length(size_mid2)),function(i){
      spheres3d(m + p_mid2$f[i],
                du, l, radius = p_mid2$s[i], 
                color = rgb((p_mid2$r)$col[[i]][1],
                            (p_mid2$r)$col[[i]][2], (p_mid2$r)$col[[i]][3],alpha = 1))
      text3d(m + p_mid2$f[i],
             du + max(as.matrix(unlist(size_mid2))) + 5,
             l,length(size_all) + length(size_low) +
               length(size_mid1) + i)
      if(ad == 0){text3d(m + p_mid2$f[i],
                         (du - max(as.matrix(unlist(size_mid2))) -5), l,
                         paste(unique(p_mid2$q1_all_a[[i]]), 
                               collapse = ""))}
    })
    create_links(p$q1_all, p_mid2$q1_all, p$f, p_mid2$f, m, du, l)
    du <- du + minespace +  max(as.matrix(unlist(size_mid2)))
  }
  m <- max(p_mid2$f) + m   
  if(dim(high_map)[1] == 1 && high_map == 0){
    d_high <- du
    p_high <- list()
    p_high$f <- 0
    size_high <- c()}
  else{
    p_high <- squize(q1_high, m1, size_high, n=n)
    du <- du + max(as.matrix(unlist(size_high)))
    d_high <- du
    l <- 30
    lans <- lapply(seq_len(length(size_high)), function(i){
      spheres3d(m + p_high$f[i], du, l,
                radius = p_high$s[i], color = rgb((p_high$r)$col[[i]][1],
                                                  (p_high$r)$col[[i]][2], (p_high$r)$col[[i]][3], alpha=1))
      text3d(m + p_high$f[i],
             du + max(as.matrix(unlist(size_high))) + 5, l,
             length(size_all) + 
               length(size_low) + length(size_mid1) +
               length(size_mid2) + i)
      if(ad == 0){text3d(m + p_high$f[i],
                         (du - max(as.matrix(unlist(size_high))) -5), l,
                         paste(unique(p_high$q1_all_a[[i]]), 
                               collapse = ""))}
    })
    create_links(p$q1_all, p_high$q1_all, p$f, p_high$f, m, du, l)
  }
  annot_right(p_high, p, l, d_high, d_mid2, d_mid1, d_low, m)
  u <- min(- minespace,-max(as.matrix(unlist(size_all))) - 10)
  lans <- lapply(seq_len(100), function(i){
    segments3d(c((((i - 1) / 100)*(max(m + max(p_high$f),
                                       max(p$f)))), ((i / 100)*(max(m + max(p_high$f),
                                                                    max(p$f))))), c(u, u), c(l, l),
               col= matlab.like2(100)[i], lwd = 20)
  })
  text3d(-1-100, u, l, "Mean Deviation")
  if(mean_value_m1 == "N"){
    text3d(-40, u, l, round(qmin / length(select),2))
    text3d(max(m + max(p_high$f),
               max(p$f)) + 40, u, l, 
           round(qmax / length(select),2))
  }
  else{
    text3d(-40, u, l, round(qmin, 2))
    text3d(max(m + max(p_high$f), max(p$f)) + 40, 
           u, l, round(qmax, 2))
  }
  return(end_out)
}
############################################################################################################
  #### descriptions files .txt
  # from now on is just an annotation of the pictures. Separate in another function.
description_files<-function(ars){  
  e <- list()
  if(length(p$q1_all) != 0){
    for(i in seq_len(length(size_all))){
          e_1 <- paste(p$q1_all[[i]][1],
                       paste(ddd[p$q1_all[[i]][1], ni], ")"), 
                       sep=" (")
          if(length(p$q1_all[[i]]) > 1){
            for(j in 2:length(p$q1_all[[i]])){
              e_1 <- paste(paste(e_1, 
                                 p$q1_all[[i]][j], 
                                 sep=","),
                           paste(ddd[p$q1_all[[i]][j], ni],
                                 ")"), 
                           sep=" (")
            }
          }
          e[[i]] <- paste(i, e_1, sep = ":" )
    }
  }
  d <- length(size_all)
  e[[d + 1]]<- "The first quartile"
  if(length(p_low$q1_all) != 0){
    for(i in seq_len(length(size_low))){
      e_1 <- paste(p_low$q1_all[[i]][1],
      paste(ddd[p_low$q1_all[[i]][1], ni], ")"),sep="(")
      if(length(p_low$q1_all[[i]]) > 1){
        for(j in 2:length(p_low$q1_all[[i]])){
          e_1 <- paste(paste(e_1,
                             p_low$q1_all[[i]][j], 
                             sep=","),
                       paste(ddd[p_low$q1_all[[i]][j], ni], ")"), 
                       sep="(")
        }
      }
      e[[i + d + 1]] <- paste(i + d, e_1, sep = ":" )
    }
  }   
  d <- length(size_all) + length(size_low) + 1
  e[[d+1]]<- "The second quartile"
  if(length(p_mid1$q1_all) != 0){
    for(i in seq_len(length(size_mid1))){
      e_1 <- paste(p_mid1$q1_all[[i]][1],
                   paste(ddd[p_mid1$q1_all[[i]][1], ni],
                         ")"),
                   sep="(")
      if(length(p_mid1$q1_all[[i]]) > 1){
        for(j in 2:length(p_mid1$q1_all[[i]])){
          e_1 <-paste(paste(e_1,
                            p_mid1$q1_all[[i]][j], 
                            sep = ","),
                      paste(ddd[p_mid1$q1_all[[i]][j], ni],")"),
                      sep="(")
        }
      }
      e[[i + d + 1]]<- paste(i + d, e_1, sep = ":" )
    }
  }  
  d <- length(size_all) + length(size_low) + length(size_mid1) + 2
  e[[d + 1]]<- "The third quartile"
  if(length(p_mid2$q1_all) != 0){
    for(i in seq_len(length(size_mid2))){
      e_1 <- paste(p_mid2$q1_all[[i]][1],
                   paste(ddd[p_mid2$q1_all[[i]][1], ni], 
                         ")"),
                    sep = "(")
      if(length(p_mid2$q1_all[[i]]) > 1){
        for(j in 2:length(p_mid2$q1_all[[i]])){
          e_1 <- paste(paste(e_1, 
                             p_mid2$q1_all[[i]][j], 
                             sep = ","),
                       paste(ddd[p_mid2$q1_all[[i]][j], ni], 
                             ")"), 
                       sep = "(")
        }
      }
      e[[i + d + 1]]<-paste(i + d, e_1, sep = ":" )
    }
  }
  d <- length(size_all) + 
       length(size_low) +
       length(size_mid1) + 
       length(size_mid2) + 
       3
  e[[d + 1]]<- "The fourth quartile"
  if(length(p_high$q1_all) != 0){
    for(i in seq_len(length(size_high))){
      e_1 <- paste(p_high$q1_all[[i]][1],
                   paste(ddd[p_high$q1_all[[i]][1], ni], 
                         ")"), 
                   sep = "(")
      if(length(p_high$q1_all[[i]]) > 1){
        for(j in 2:length(p_high$q1_all[[i]])){
          e_1 <- paste(paste(e_1, 
                             p_high$q1_all[[i]][j], 
                             sep = ","),
                       paste(ddd[p_high$q1_all[[i]][j], ni], 
                             ")"), 
                       sep = "(")
        }
      }
      e[[i + d + 1]]<-paste(i + d, e_1, sep = ":" )
    }
  }
  # here I write
  write.table(e, 
              file=paste(filename,
                         "description.txt", 
                         sep = "_"),
              quote = FALSE, 
              sep = "\n", 
              row.names = FALSE,
              col.names = FALSE)
  end_out <- list(low_map = p_low$q1_all,
                  mid1_map=p_mid1$q1_all,
                  mid2_map=p_mid2$q1_all, 
                  high_map = p_high$q1_all,
                  all = p$q1_all)
  vari <- setdiff(ls(), list(beg)[[1]])
  vari <- setdiff(list(vari)[[1]], "end_out")
  rm(list = vari)
  return(end_out)
}
