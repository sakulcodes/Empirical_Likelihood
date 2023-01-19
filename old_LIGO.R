# estimate cell probabilities for (M,L) as requested by Yeunhwan
# return a matrix S (one for each LIGO star) with S_11 = P(l_1 < L < l_2, m_1 < M < m_2), S_12 = P((l_1 < L < l_2, m_2 < M < m_3), ...

# create a grid of eval.points E (the mid-points of each recatngle) to feed into kde just like for the mass-radius case
# grid for L (tidal deformability) goes from 0 to 1500 with a spacing of 15



# -------------------------------------- #
# load libraries
library(ks)
library(ggplot2)

# stuff related to L grid
delta_L = 15
Lgrid = seq(0,1500,by=delta_L)
lower_midpt_L = 0 + delta_L/2
upper_midpt_L = 1500 - delta_L/2
Lgrid_midpts = seq(lower_midpt_L,upper_midpt_L,by=delta_L)
n_L = length(Lgrid_midpts)

# stuff related to M grid
delta_M = 0.02
Mgrid = seq(1,2,by=delta_M)
lower_midpt_M = 1 + delta_M/2
upper_midpt_M = 2 - delta_M/2
Mgrid_midpts = seq(lower_midpt_M,upper_midpt_M,by=delta_M)
n_M = length(Mgrid_midpts)

# make E matrix
n_eval = n_L*n_M                           # total number of intervals (also, # points where kde is evaluated)
Ev = matrix(0,nrow=n_eval,ncol=2)
cter = 0

## vectorized 2D grid for (L1, M1) OR (L2, M2)

for (ll in 1:n_L)
{
  for (mm in 1:n_M)
  {
    cter = cter+1
    Ev[cter,] = c(Lgrid_midpts[ll],Mgrid_midpts[mm])
  }
}

## vectorized 4D grid for (L1, M1) AND (L2, M2)
cter = 0
Ev4 = matrix(0,nrow=n_eval^2,ncol=4)
for (l1 in 1:n_L)
{
  for (m1 in 1:n_M)
  {
    for (l2 in 1:n_L)
    {
      for (m2 in 1:n_M)
      {
        cter = cter+1
        Ev4[cter,] = c(Lgrid_midpts[l1],Mgrid_midpts[m1],Lgrid_midpts[l2],Mgrid_midpts[m2])
      }
    }
  }
}

# -- read LIGO data -- #
ligo <- read.table("EoS-insensitive_posterior_samples.dat.txt",header = TRUE)
colnames(ligo) <-c("M1","M2", "L1","L2", "R1", "R2")


# -- Create S matrix for first LIGO data -- #
LM_ligo_1 = cbind(ligo$L1,ligo$M1)

H_1 = Hscv(LM_ligo_1)                    ## smoothed cv bandwidth matrix (lscv gives warning due to duplicated values)
kRp = kde(LM_ligo_1, H = H_1)
plot(kRp,cont = c(5,10,25,50,68,90,95),main="LIGO1",xlab="Tidal deformability",ylab="Mass",xlim=c(0,1500),ylim=c(1,2))           ## contour plot of fitted density

## Extract 90% contour coordinates
# contour_90_ligo_1 <- with(kRp, contourLines(x=eval.points[[1]], y=eval.points[[2]],
#                                     z=estimate, levels=cont["10%"])[[1]])
# LM90_ligo_1= as.matrix(cbind(contour_90_ligo_1$x,contour_90_ligo_1$y))
# write.csv(LM90_ligo_1,"LM90_ligo_1.csv")

# contour_95_ligo_1 <- with(kRp, contourLines(x=eval.points[[1]], y=eval.points[[2]],
#                                             z=estimate, levels=cont["5%"])[[1]])
# LM95_ligo_1= as.matrix(cbind(contour_95_ligo_1$x,contour_95_ligo_1$y))
# write.csv(LM95_ligo_1,"LM95_ligo_1.csv")

contour_68_ligo_1 <- with(kRp, contourLines(x=eval.points[[1]], y=eval.points[[2]],
                                            z=estimate, levels=cont["32%"])[[1]])
LM68_ligo_1= as.matrix(cbind(contour_68_ligo_1$x,contour_68_ligo_1$y))
write.csv(LM68_ligo_1,"LM68_ligo_1.csv")

df_contour_ligo_1 <- data.frame(contour_90_ligo_1)
d_ligo_1 <-data.frame(L=ligo$L1, M=ligo$M1)
p <- ggplot(data=d_ligo_1, aes(L, M)) +
  geom_point() +
  geom_path(aes(x,y), data=df_contour_ligo_1) +
  theme_bw()

### now kde

kR = kde(LM_ligo_1, H = H_1, eval.points = Ev)
S_ligo_1_vect = (delta_L*delta_M)*kR$estimate
S_ligo_1_vect = S_ligo_1_vect/sum(S_ligo_1_vect)
S_ligo_1 = matrix(S_ligo_1_vect,nrow=n_L,ncol=n_M,byrow = TRUE)     #n_L x n_M matrix of output cell probs
S_ligo_1[S_ligo_1<=0] = 0
write.csv(S_ligo_1,"S_ligo_1.csv")

# some plotting
L = Lgrid_midpts
M = Mgrid_midpts
data <- expand.grid(X=L,Y=M)
data$S = as.vector(c(S_ligo_1))
ggplot(data, aes(X, Y, fill= S)) + geom_tile()


# -- Create S matrix for second LIGO data -- #
LM_ligo_2 = cbind(ligo$L2,ligo$M2)

H_2 = Hscv(LM_ligo_2)                    ## smoothed cv bandwidth matrix (lscv gives warning due to duplicated values)
kRp = kde(LM_ligo_2, H = H_2)
plot(kRp,cont = c(5,10,25,50,68,90,95),main="LIGO2",xlab="Tidal deformability",ylab="Mass",xlim=c(0,1500),ylim=c(1,2))           ## contour plot of fitted density

## Extract 90% contour coordinates
# contour_90_ligo_2 <- with(kRp, contourLines(x=eval.points[[1]], y=eval.points[[2]],
#                                             z=estimate, levels=cont["10%"])[[1]])
# LM90_ligo_2 = as.matrix(cbind(contour_90_ligo_2$x,contour_90_ligo_2$y))
# write.csv(LM90_ligo_2,"LM90_ligo_2.csv")

# contour_95_ligo_2 <- with(kRp, contourLines(x=eval.points[[1]], y=eval.points[[2]],
#                                             z=estimate, levels=cont["5%"])[[1]])
# LM95_ligo_2 = as.matrix(cbind(contour_95_ligo_2$x,contour_95_ligo_2$y))
# write.csv(LM95_ligo_2,"LM95_ligo_2.csv")

contour_68_ligo_2 <- with(kRp, contourLines(x=eval.points[[1]], y=eval.points[[2]],
                                            z=estimate, levels=cont["32%"])[[1]])
LM68_ligo_2 = as.matrix(cbind(contour_68_ligo_2$x,contour_68_ligo_2$y))
write.csv(LM68_ligo_2,"LM68_ligo_2.csv")




df_contour_ligo_2 <- data.frame(contour_90_ligo_2)
d_ligo_2 <-data.frame(L=ligo$L2, M=ligo$M2)
p <- ggplot(data=d_ligo_2, aes(L, M)) +
  geom_point() +
  geom_path(aes(x,y), data=df_contour_ligo_2) +
  theme_bw()



kR = kde(LM_ligo_2, H = H_2, eval.points = Ev)
S_ligo_2_vect = (delta_L*delta_M)*kR$estimate
S_ligo_2_vect = S_ligo_2_vect/sum(S_ligo_2_vect)
S_ligo_2 = matrix(S_ligo_2_vect,nrow=n_L,ncol=n_M,byrow = TRUE)     #n_L x n_M matrix of output cell probs
S_ligo_2[S_ligo_2<=0] = 0
write.csv(S_ligo_2,"S_ligo_2.csv")

# some plotting
L = Lgrid_midpts
M = Mgrid_midpts
data <- expand.grid(X=L,Y=M)
data$S = as.vector(c(S_ligo_2))
ggplot(data, aes(X, Y, fill= S)) + geom_tile()

#### 4-variate KDE for L1, M1, L2, M2 ####

LM_ligo = cbind(LM_ligo_1,LM_ligo_2)
H = Hscv(LM_ligo)
kR = kde(LM_ligo, H = H, eval.points = Ev4)
S_ligo_vect = ((delta_L*delta_M)^2)*kR$estimate
S_ligo_vect = S_ligo_vect/sum(S_ligo_vect)
S_ligo = array(S_ligo_vect,c(n_L,n_M,n_L,n_M))     #(n_L x n_M)^2 array of output cell probs
S_ligo[S_ligo<=0] = 0



### Arrange the 4D array using reshape
S_ligo_melt<-melt(S_ligo)
names(S_ligo_melt) <-c("L1", "M1", "L2","M2", "prob")
write.csv(S_ligo_melt,"S_ligo.csv",row.names = FALSE)
