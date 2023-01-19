# estimate cell probabilities as requested by Yeunhwan
# return a matrix S (one for Riley and Miller each) with S_11 = P(9<R<9.02,1<M<1.01), S_12 = P(9<R<9.02,1.01<M<1.02), ...
# use the kernel density estimate to approximate these probabilities. for example, S_11 ~ fhat(9.01,1.005) x 0.02 x 0.01 etc
# create a grid of eval.points E (the mid-points of each recatngle) to feed into kde
# Ev = [9.01   1.005
#       9.01   1.015
#       .       .
#       .       .
#       9.01   1.995
#       9.03   1.005
#       .       .
#       .       .
#       9.03   1.995
#       .       .
#       .       .]

# Update on Jul 08 --> automate code to handle arbitrary values of delta_M and delta_R


# -------------------------------------- #
# load libraries
library(ks)
library(ggplot2)
library(reshape2)

# stuff related to R grid
#delta_R = 0.02
delta_R = 0.05
Rgrid = seq(9,16,by=delta_R)
#Rgrid_midpts = seq(9.01,15.99,by=delta_R)
lower_midpt_R = 9 + delta_R/2
upper_midpt_R = 16 - delta_R/2
Rgrid_midpts = seq(lower_midpt_R,upper_midpt_R,by=delta_R)
n_R = length(Rgrid_midpts)

# stuff related to M grid
#delta_M = 0.01
delta_M = 0.02
Mgrid = seq(1,2,by=delta_M)
lower_midpt_M = 1 + delta_M/2
upper_midpt_M = 2 - delta_M/2
#Mgrid_midpts = seq(1.005,1.995,by=delta_M)
Mgrid_midpts = seq(lower_midpt_M,upper_midpt_M,by=delta_M)
n_M = length(Mgrid_midpts)

# make E matrix
n_eval = n_R*n_M                           # total number of intervals (also, # points where kde is evaluated)
Ev = matrix(0,nrow=n_eval,ncol=2)
cter = 0

for (rr in 1:n_R)
{
  for (mm in 1:n_M)
  {
    cter = cter+1
    Ev[cter,] = c(Rgrid_midpts[rr],Mgrid_midpts[mm])
  }
}


# -- Create S matrix for Riley data -- #
RM_Riley_all <- read.delim("riley/ST_PST/run1/run1_nlive1000_eff0.3_noCONST_noMM_noIS_tol-1post_equal_weights.dat", sep=" ", header = FALSE)
RM_Riley_cols <-RM_Riley_all[,c(13,9)]
colnames(RM_Riley_cols) <-c("R","M")        ## column 1 is radius, column 2 is mass
RM_Riley_cols1= na.omit(RM_Riley_cols)      ## omit NA, in data.frame format
RM_Riley = as.matrix(RM_Riley_cols1)        ## convert to matrix

H_Riley = Hscv(RM_Riley)                    ## smoothed cv bandwidth matrix (lscv gives warning due to duplicated values)
kRp = kde(RM_Riley, H = H_Riley)
plot(kRp,cont = c(5,10,25,50,68,90,95),main="Riley",xlab="Radius",ylab="Mass",xlim=c(9,16),ylim=c(1,2))           ## contour plot of fitted density

#### extract 90% contour coordinates
#contour_90_Riley <- with(kRp, contourLines(x=eval.points[[1]], y=eval.points[[2]],
#                                            z=estimate, levels=cont["10%"])[[1]])
#RM90_Riley= as.matrix(cbind(contour_90_Riley$x,contour_90_Riley$y))
#write.csv(RM90_Riley,"RM90_Riley.csv")

# contour_95_Riley <- with(kRp, contourLines(x=eval.points[[1]], y=eval.points[[2]],
#                                            z=estimate, levels=cont["5%"])[[1]])
# RM95_Riley= as.matrix(cbind(contour_95_Riley$x,contour_95_Riley$y))
# write.csv(RM95_Riley,"RM95_Riley.csv")

contour_68_Riley <- with(kRp, contourLines(x=eval.points[[1]], y=eval.points[[2]],
                                           z=estimate, levels=cont["32%"])[[1]])
RM68_Riley= as.matrix(cbind(contour_68_Riley$x,contour_68_Riley$y))
write.csv(RM68_Riley,"RM68_Riley.csv")

df_contour_Riley <- data.frame(contour_90_Riley)
d_Riley <-data.frame(R=RM_Riley[,1], M=RM_Riley[,2])
p <- ggplot(data=d_Riley, aes(R, M)) +
  geom_point(color="grey") +
  geom_path(aes(x,y), data=df_contour_Riley) +
  theme_bw()

##### kde etc

kR = kde(RM_Riley, H = H_Riley, eval.points = Ev)
S_Riley_vect = (delta_R*delta_M)*kR$estimate
S_Riley_vect = S_Riley_vect/sum(S_Riley_vect)
S_Riley = matrix(S_Riley_vect,nrow=n_R,ncol=n_M,byrow = TRUE)     #n_R x n_M matrix of output cell probs
S_Riley[S_Riley<=0] = 0
#write.csv(S_Riley,"S_Riley.csv")
write.csv(S_Riley,"S_Riley_coarsegrid.csv")

# some plotting
R = Rgrid_midpts
M = Mgrid_midpts
data <- expand.grid(X=R,Y=M)
data$S = as.vector(c(S_Riley))
ggplot(data, aes(X, Y, fill= S)) + geom_tile()


# -- create S matrix for Miller data -- #
RM_Miller <- read.table("miller/J0030_3spot_RM.txt",header = FALSE)
colnames(RM_Miller) <-c("R","M")

H_Miller = 4*Hscv(RM_Miller)                    ## smoothed cv bandwidth matrix (lscv gives warning due to duplicated values)
kMp = kde(RM_Miller, H = H_Miller)
#kMp = kde(RM_Miller)
plot(kMp,cont = c(5,10,25,50,68,90,95),main="Miller",xlab="Radius",ylab="Mass",xlim=c(9,16),ylim=c(1,2))           ## contour plot of fitted density

#### extract 90% contour coordinates
# contour_90_Miller <- with(kMp, contourLines(x=eval.points[[1]], y=eval.points[[2]],
#                                            z=estimate, levels=cont["10%"])[[1]])
# RM90_Miller= as.matrix(cbind(contour_90_Miller$x,contour_90_Miller$y))
# write.csv(RM90_Miller,"RM90_Miller.csv")

contour_95_Miller <- with(kMp, contourLines(x=eval.points[[1]], y=eval.points[[2]],
                                            z=estimate, levels=cont["5%"])[[1]])
RM95_Miller= as.matrix(cbind(contour_95_Miller$x,contour_95_Miller$y))
#write.csv(RM95_Miller,"RM95_Miller.csv")

contour_68_Miller <- with(kMp, contourLines(x=eval.points[[1]], y=eval.points[[2]],
                                            z=estimate, levels=cont["32%"])[[1]])
RM68_Miller= as.matrix(cbind(contour_68_Miller$x,contour_68_Miller$y))
#write.csv(RM68_Miller,"RM68_Miller.csv")

df_contour_95_Miller <- data.frame(contour_95_Miller)
df_contour_68_Miller <- data.frame(contour_68_Miller)

d_Miller <-data.frame(R=RM_Miller[,1], M=RM_Miller[,2])
p <- ggplot(data=d_Miller, aes(R, M)) +
  geom_point(color="grey") +
  geom_path(aes(x,y), data=df_contour_68_Miller) +
  geom_path(aes(x,y), data=df_contour_95_Miller) +
  theme_bw()


#### kde etc

kM = kde(RM_Miller, H = H_Miller, eval.points = Ev)
S_Miller_vect = (delta_R*delta_M)*kM$estimate
S_Miller_vect = S_Miller_vect/sum(S_Miller_vect)
S_Miller = matrix(S_Miller_vect,nrow=n_R,ncol=n_M,byrow = TRUE)                     # n_R x n_M matrix of output cell probs
S_Miller[S_Miller<=0] = 0
#write.csv(S_Miller,"S_Miller.csv")
write.csv(S_Miller,"S_Miller_coarsegrid.csv")

# some plotting
R = Rgrid_midpts
M = Mgrid_midpts
data <- expand.grid(X=R,Y=M)
data$S = as.vector(c(S_Miller))
ggplot(data, aes(X, Y, fill= S)) + geom_tile()



### Pooled data likelihood
RM_Miller2<- RM_Miller[sample(1:10000,12242,replace=TRUE),]
pool_miller_riley <- rbind(RM_Miller2,RM_Riley)
H_pool= Hscv(pool_miller_riley)                    ## smoothed cv bandwidth matrix (lscv gives warning due to duplicated values)
kRp = kde(pool_miller_riley, H = 3*H_pool)
plot(kRp,cont = c(5,10,25,50,68,90,95),main="pool",xlab="Radius",ylab="Mass",xlim=c(9,16),ylim=c(1,2))           ## contour plot of fitted density

##### pooled kde etc

kR = kde(pool_miller_riley, H =3*H_pool, eval.points = Ev)
S_pool_vect = (delta_R*delta_M)*kR$estimate
S_pool_vect = S_pool_vect/sum(S_pool_vect)
S_pool = matrix(S_pool_vect,nrow=n_R,ncol=n_M,byrow = TRUE)     #n_R x n_M matrix of output cell probs
S_pool[S_pool<=0] = 0
#write.csv(S_Riley,"S_Riley.csv")
write.csv(S_pool,"S_pool_coarsegrid_eq_sam.csv")

##### kde just using miller data

H_miller = Hscv(RM_Miller)

kR = kde(RM_Miller, H =3*H_miller, eval.points = Ev)
S_miller_vect = (delta_R*delta_M)*kR$estimate
S_miller_vect = S_pool_vect/sum(S_pool_vect)
S_miller = matrix(S_pool_vect,nrow=n_R,ncol=n_M,byrow = TRUE)     #n_R x n_M matrix of output cell probs
S_miller[S_miller<=0] = 0
#write.csv(S_Riley,"S_Riley.csv")
write.csv(S_pool,"S_Miller_coarsegrid.csv")
