library(airway)
library(colorRamps)
data(airway)
dw <- rowSums(assay(airway))>80
dw <- names(dw[dw==TRUE])
airway <- airway[dw,]
assay(airway) <- log(assay(airway)+1,2)
experiment <- TTMap::make_matrices(airway,seq_len(4),seq_len(4)+4,
                                     NAME = rownames(airway), CLID =rownames(airway))
E=1
Pw=1.1
Bw=0
TTMAP_part1prime <-TTMap::control_adjustment(normal.pcl = experiment$CTRL,
                                              ALPHA,
                                              tumor.pcl = experiment$TEST,
                                              normalname = "The_healthy_controls", dataname = "The_effect_of_cancer",
                                              org.directory = getwd(), e=E,P=Pw,B=Bw);
TTMAP_part1_hda <- TTMap::hyperrectangle_deviation_assessment(x = TTMAP_part1prime,k=dim(TTMAP_part1prime$Normal.mat)[2],
                                                              dataname = "The_effect_of_cancer", normalname = "The_healthy_controls");
head(TTMAP_part1_hda$Dc.Dmat)
library(rgl)
ALPHA <- 1
annot <- c(paste(colnames(experiment$TEST[,-seq_len(3)]),"Dis",sep="."),
           paste(colnames(experiment$CTRL[,-seq_len(3)]),"Dis",sep="."))
annot <- cbind(annot,annot)
rownames(annot)<-annot[,1]
dd5_sgn_only<-TTMap::generate_mismatch_distance(TTMAP_part1_hda,
                                                select=rownames(TTMAP_part1_hda$Dc.Dmat),alpha = ALPHA)
TTMAP_part1_hda$m['SRR1039516.Dis']=4
TTMAP_part1_hda$m['SRR1039517.Dis']=3
TTMAP_part1_hda$m['SRR1039520.Dis']=2
TTMAP_part1_hda$m['SRR1039521.Dis']=100
clusters_by_level<-get_clusters_by_level(TTMAP_part1_hda,
                                        TTMAP_part1_hda$m,
                                        select=rownames(TTMAP_part1_hda$Dc.Dmat),
                                        ddd=annot,
                                        e=TTMap::calcul_e(dd5_sgn_only,0.95,TTMAP_part1prime,1),
                                        filename="first_comparison",
                                        n=1,
                                        dd=dd5_sgn_only)
q1_and_sizes_by_level <- lapply(clusters_by_level, calculate_size, dd=dd5_sgn_only)
names(q1_and_sizes_by_level) <- c("all", "low", "mid1", "mid2", "high")
open3d()
TTMAP_part2_gtlmap <- draw(q1_and_sizes_by_level,
                           TTMAP_part1_hda,
                           TTMAP_part1_hda$m,
                           select=rownames(TTMAP_part1_hda$Dc.Dmat),
                           ddd=annot,
                           e=TTMap::calcul_e(dd5_sgn_only,0.95,TTMAP_part1prime,1),
                           filename="first_comparison",
                           n=1,
                           dd=dd5_sgn_only,
                           qmin=min(m1),
                           qmax=max(m1),
                           level=0)
TTMAP_part1_hda$m['SRR1039516.Dis']=4
TTMAP_part1_hda$m['SRR1039517.Dis']=3
TTMAP_part1_hda$m['SRR1039520.Dis']=2
TTMAP_part1_hda$m['SRR1039521.Dis']=1
clusters_by_level<-get_clusters_by_level(TTMAP_part1_hda,
                                         TTMAP_part1_hda$m,
                                         select=rownames(TTMAP_part1_hda$Dc.Dmat),
                                         ddd=annot,
                                         e=TTMap::calcul_e(dd5_sgn_only,0.95,TTMAP_part1prime,1),
                                         filename="first_comparison",
                                         n=1,
                                         dd=dd5_sgn_only)
q1_and_sizes_by_level <- lapply(clusters_by_level, calculate_size, dd=dd5_sgn_only)
names(q1_and_sizes_by_level) <- c("all", "low", "mid1", "mid2", "high")
TTMAP_part2.1_gtlmap <- draw(q1_and_sizes_by_level,
                           TTMAP_part1_hda,
                           TTMAP_part1_hda$m*200,
                           select=rownames(TTMAP_part1_hda$Dc.Dmat),
                           ddd=annot,
                           e=TTMap::calcul_e(dd5_sgn_only,0.95,TTMAP_part1prime,1),
                           filename="first_comparison",
                           n=1,
                           dd=dd5_sgn_only,
                           qmin=min(m1),
                           qmax=max(m1),
                           level=1)
rgl.postscript("first_output.pdf","pdf")
TTMap::ttmap_sgn_genes(TTMAP_part2_gtlmap,
                       TTMAP_part1_hda, TTMAP_part1prime,
                       annot, n = 2, a = ALPHA,
                       filename = "first_trial", annot = TTMAP_part1prime$tag.pcl, col = "NAME",
                       path = getwd(), Relaxed = 0)
