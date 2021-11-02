#library(rEDM)
#data(sardine_anchovy_sst)
find_E<-function(ce_data){
  sim_1<- simplex(ce_data[,1],lib=c(1,0.6*nrow(ce_data)),pred=c(0.6*nrow(ce_data),nrow(ce_data)),E=c(2:11))
  E_1 <-sim_1[which.max(sim_1$rho),"E"][1]
  sim_2<- simplex(ce_data[,2],lib=c(1,0.6*nrow(ce_data)),pred=c(0.6*nrow(ce_data),nrow(ce_data)),E=c(2:11))
  E_2 <-sim_2[which.max(sim_2$rho),"E"][1]
  return(c(E_1,E_2))
}

sur_results<-function(ce_data,E_1,E_2){
  sur<- make_surrogate_data(ce_data[,1],method ="ebisuzaki",num_surr = 100)
  sur_rho<- do.call(cbind, lapply(seq_len(ncol(sur)),
                                  function(i){
                                    a <- cbind(sur[, i], ce_data[,2])
                                    ccm_infer(a,E_1,E_2)[[2]]
                                  }))
  return(sur_rho)}
ccm_infer<-function(ce_data,E_1,E_2){
  my_lib_sizes <- c(seq(5, 20, by = 5), seq(21, nrow(ce_data), by = 20))
  a_xmap_b<- ccm(ce_data, E = E_1, lib_column = 1,
                 target_column = 2, lib_sizes = my_lib_sizes, num_samples = 30,
                 random_libs = TRUE, replace = TRUE,silent = TRUE)
  b_xmap_a<- ccm(ce_data, E = E_2, lib_column = 2, target_column = 1,
                 lib_sizes = my_lib_sizes, num_samples = 30, random_libs = TRUE,
                 replace = TRUE,silent = TRUE)
  a_xmap_b_means <- ccm_means(a_xmap_b)
  b_xmap_a_means <- ccm_means(b_xmap_a)
  par(mar = c(4, 4, 1, 1), mgp = c(2.5, 1, 0),mfrow = c(1, 1))
  y1 <- pmax(0, a_xmap_b_means$rho)
  y2 <- pmax(0, b_xmap_a_means$rho)
  list(y1,y2,my_lib_sizes)
}

plot_ccm<-function(ce_data){
  E<-find_E(ce_data)
  results<-ccm_infer(ce_data,E[1],E[2])
  sur<-sur_results(ce_data,E[1],E[2])
  mrow<-length(results[[3]])
  ind<-order(sur[mrow,],decreasing=TRUE)[0.05*ncol(sur)+1]
  sig_data<-sur[,ind]
  plot(results[[3]],results[[1]], type = "l", col = "red", xlab = "Library Size",
       ylab = "Cross Map Skill (rho)", ylim = c(0, 0.5))
  lines(results[[3]], results[[2]], col = "blue")
  legend(x = "topleft", legend = c(paste(colnames(ce_data)[1],'xmap',colnames(ce_data)[2]), paste(colnames(ce_data)[2],'xmap',colnames(ce_data)[1])), col = c("red","blue"), lwd = 1, bty = "n", inset = 0.02, cex = 0.8)
  cord.x <-c(5,results[[3]],max(results[[3]]))
  cord.y <- c(0,sig_data,0)
  polygon(cord.x,cord.y,border = NA,density = 10)
}
#plot_ccm(cbind(sardine_anchovy_sst$np_sst,sardine_anchovy_sst$anchovy))




