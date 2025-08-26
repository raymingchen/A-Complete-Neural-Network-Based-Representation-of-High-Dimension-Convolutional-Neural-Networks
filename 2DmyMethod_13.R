####codes for my manuscript - 2D image implementation via linearization method
rm(list=ls());
############ import related packages #####################
library(recolorize);  ##image to matrix
library(raster); ##matrix to image
library(imager); ##to show image to image
library(OpenImageR); ## to rotate an image

############ key functions ###############################
##To begin with, let us implement the bijective function in Lemma 4 and Corollary 7
TensorIndex_to_Num=function(Nvec,ivec)  ##Nvec is the domain;ivec is a vector in prod. set of Nvec
{
M=length(ivec);
  for(k in 2:M)
 {
 ab=(ivec[k]-1)*prod(Nvec[1:(k-1)]);
 }
label_to=ivec[1]+ab;
return(label_to);
}

##trial 1 & 2
Nvec=c(14,27);
ivec=c(12,16);
tensorIndex_to_Num=TensorIndex_to_Num(Nvec,ivec);tensorIndex_to_Num;

Nvec=c(148,297);
ivec=c(112,194);
tensorIndex_to_Num=TensorIndex_to_Num(Nvec,ivec);tensorIndex_to_Num;

## calculate the inverse function Num_to_TensorIndex of TensorIndex_to_Num
## it suffices to take the tensor as matrix in our case
Num_to_MatrixIndex=function(Nvec,n)  ##Nvec is the domain;ivec is a vector in prod. set of Nvec
{
rc=rep(0,2);  ##rc=(row,column) index
m1=Nvec[1];m2=Nvec[2];
quo=n%/%m1;
rem=n%%m1;
  if(rem==0){rc=c(m1,quo)}else{rc=c(rem,quo+1)};
return(rc);
}

##trial 1 & 2
Nvec=c(14,27);
n=222;
num_to_MatrixIndex=Num_to_MatrixIndex(Nvec,n); num_to_MatrixIndex;

Nvec=c(148,297);
n=28676;
num_to_MatrixIndex=Num_to_MatrixIndex(Nvec,n); num_to_MatrixIndex;


###################### main body ##########################################
##image to matrix
Puppy=readImage(file.choose());dim(Puppy);imageShow(Puppy); ##dim=408x531x3
Puppy_R=Puppy[,,1];Puppy_R;dim(Puppy_R); ##dim_Red=408x531;
Puppy_R=resizeImage(Puppy_R, width = 200, height = 200, method = 'bilinear', normalize_pixels = TRUE)
imageShow(Puppy_R);
##Puppy_G=Puppy[,,2];Puppy_G;dim(Puppy_G);      ##dim_Green=408x531;
##Puppy_B=Puppy[,,3];Puppy_B;dim(Puppy_B);       ##dim_Blue=408x531;

##input tensor
InputTensor=Puppy_R;imageShow(InputTensor);
TenDim=dim(Puppy_R);TenDim;

##(radnom kernels)
F1=matrix(runif(36,-100,100),6,6);F;imageShow(F1);
f1=dim(F1);f1;
F2=matrix(runif(16,-8,8),4,4);imageShow(F2);
f2=c(4,4);
F3=diag(c(1,1,1,1,1));imageShow(F3);
f3=c(5,5);
F4=matrix(0,4,4);F4[1,]=F4[4,]=1;F4[2,]=F4[3,]=0;imageShow(F4);
f4=c(4,4);
##stride vector;
sig=c(3,5);

##power tensor and convolution: with kernel/filter F1
pf1_1=floor((TenDim[1]-f1[1])/sig[1])+1;pf1_1;
pf1_2=floor((TenDim[2]-f1[2])/sig[2])+1;pf1_2;
P_f1=as.list(numeric(pf1_1*pf1_2));dim(P_f1)=c(pf1_1,pf1_2);
Ind_f1=as.list(numeric(pf1_1*pf1_2));dim(Ind_f1)=c(pf1_1,pf1_2);
Index_f1=as.list(numeric(pf1_1*pf1_2));dim(Index_f1)=c(pf1_1,pf1_2);
encodedIndex_f1=as.list(numeric(pf1_1*pf1_2));dim(encodedIndex_f1)=c(pf1_1,pf1_2);
conv1=matrix(0,pf1_1,pf1_2);conv1;dim(conv1);
ma1=matrix(0,2,f1[1]);
for(i in 1:pf1_1)
{
v1=(1+sig[1]*(i-1)):(1+sig[1]*(i-1)+f1[1]-1);
   for(j in 1:pf1_2)
   {
   w1=(1+sig[2]*(j-1)):(1+sig[2]*(j-1)+f1[2]-1);
   P_f1[[i,j]]=Puppy_R[v1,w1];
   ma1[1,]=v1;ma1[2,]=w1;Ind_f1[[i,j]]=ma1;
   Index_f1[[i,j]]=as.matrix(expand.grid(v1,w1));
   mma1=matrix(0,dim(Index_f1[[i,j]])[1],1);
        for(k in 1:dim(Index_f1[[i,j]])[1])
        {
        mma1[k]=TensorIndex_to_Num(TenDim,Index_f1[[i,j]][k,])
        }
   encodedIndex_f1[[i,j]]=mma1;
   conv1[i,j]=sum(P_f1[[i,j]]*F1);
   }
}  
imageShow(conv1);

##power tensor and convolution: with kernel/filter F2
pf2_1=floor((TenDim[1]-f2[1])/sig[1])+1;pf2_1;
pf2_2=floor((TenDim[2]-f2[2])/sig[2])+1;pf2_2;
P_f2=as.list(numeric(pf2_1*pf2_2));dim(P_f2)=c(pf2_1,pf2_2);
Ind_f2=as.list(numeric(pf2_1*pf2_2));dim(Ind_f2)=c(pf2_1,pf2_2);
Index_f2=as.list(numeric(pf2_1*pf2_2));dim(Index_f2)=c(pf2_1,pf2_2);
encodedIndex_f2=as.list(numeric(pf2_1*pf2_2));dim(encodedIndex_f2)=c(pf2_1,pf2_2);
conv2=matrix(0,pf2_1,pf2_2);conv2;dim(conv2);
ma2=matrix(0,2,f2[1]);
for(i in 1:pf2_1)
{
v2=(1+sig[1]*(i-1)):(1+sig[1]*(i-1)+f2[1]-1);
   for(j in 1:pf2_2)
   {
   w2=(1+sig[2]*(j-1)):(1+sig[2]*(j-1)+f2[2]-1);
   P_f2[[i,j]]=Puppy_R[v2,w2];
   ma2[1,]=v2;ma2[2,]=w2;Ind_f2[[i,j]]=ma2;
   Index_f2[[i,j]]=as.matrix(expand.grid(v2,w2));
   mma2=matrix(0,dim(Index_f2[[i,j]])[1],1);
        for(k in 1:dim(Index_f2[[i,j]])[1])
        {
        mma2[k]=TensorIndex_to_Num(TenDim,Index_f2[[i,j]][k,])
        }
   encodedIndex_f2[[i,j]]=mma2;
   conv2[i,j]=sum(P_f2[[i,j]]*F2);
   }
}  
imageShow(conv2);

##power tensor and convolution: with kernel/filter F3
pf3_1=floor((TenDim[1]-f3[1])/sig[1])+1;pf3_1;
pf3_2=floor((TenDim[2]-f3[2])/sig[2])+1;pf3_2;
P_f3=as.list(numeric(pf3_1*pf3_2));dim(P_f3)=c(pf3_1,pf3_2);
Ind_f3=as.list(numeric(pf3_1*pf3_2));dim(Ind_f3)=c(pf3_1,pf3_2);
Index_f3=as.list(numeric(pf3_1*pf3_2));dim(Index_f3)=c(pf3_1,pf3_2);
encodedIndex_f3=as.list(numeric(pf3_1*pf3_2));dim(encodedIndex_f3)=c(pf3_1,pf3_2);
conv3=matrix(0,pf3_1,pf3_2);conv3;dim(conv3);
ma3=matrix(0,2,f3[1]);
for(i in 1:pf3_1)
{
v3=(1+sig[1]*(i-1)):(1+sig[1]*(i-1)+f3[1]-1);
   for(j in 1:pf3_2)
   {
   w3=(1+sig[2]*(j-1)):(1+sig[2]*(j-1)+f3[2]-1);
   P_f3[[i,j]]=Puppy_R[v3,w3];
   ma3[1,]=v3;ma3[2,]=w3;Ind_f3[[i,j]]=ma3;
   Index_f3[[i,j]]=as.matrix(expand.grid(v3,w3));
   mma3=matrix(0,dim(Index_f3[[i,j]])[1],1);
        for(k in 1:dim(Index_f3[[i,j]])[1])
        {
        mma3[k]=TensorIndex_to_Num(TenDim,Index_f3[[i,j]][k,])
        }
   encodedIndex_f3[[i,j]]=mma3;
   conv3[i,j]=sum(P_f3[[i,j]]*F3);
   }
}  
imageShow(conv3);

##power tensor and convolution: with kernel/filter F4
pf4_1=floor((TenDim[1]-f3[1])/sig[1])+1;pf4_1;
pf4_2=floor((TenDim[2]-f3[2])/sig[2])+1;pf4_2;
P_f4=as.list(numeric(pf4_1*pf4_2));dim(P_f4)=c(pf4_1,pf4_2);
Ind_f4=as.list(numeric(pf4_1*pf4_2));dim(Ind_f4)=c(pf4_1,pf4_2);
Index_f4=as.list(numeric(pf4_1*pf4_2));dim(Index_f4)=c(pf4_1,pf4_2);
encodedIndex_f4=as.list(numeric(pf4_1*pf4_2));dim(encodedIndex_f4)=c(pf4_1,pf4_2);
conv4=matrix(0,pf4_1,pf4_2);conv4;dim(conv4);
ma4=matrix(0,2,f4[1]);
for(i in 1:pf4_1)
{
v4=(1+sig[1]*(i-1)):(1+sig[1]*(i-1)+f4[1]-1);
   for(j in 1:pf4_2)
   {
   w4=(1+sig[2]*(j-1)):(1+sig[2]*(j-1)+f4[2]-1);
   P_f4[[i,j]]=Puppy_R[v4,w4];
   ma4[1,]=v4;ma4[2,]=w4;Ind_f4[[i,j]]=ma4;
   Index_f4[[i,j]]=as.matrix(expand.grid(v4,w4));
   mma4=matrix(0,dim(Index_f4[[i,j]])[1],1);
        for(k in 1:dim(Index_f4[[i,j]])[1])
        {
        mma4[k]=TensorIndex_to_Num(TenDim,Index_f4[[i,j]][k,])
        }
   encodedIndex_f4[[i,j]]=mma4;
   conv4[i,j]=sum(P_f4[[i,j]]*F4);
   }
}  
imageShow(conv4);

### our linearization method

##(step1) Linearize the input tensor to form the nodes/content in the 1st input layer
vecInputTensor=as.vector(InputTensor);vecInputTensor;  
nodes1stLayer=length(vecInputTensor);nodes1stLayer;  ##200*200=40,000.

##(step2) Calculate the number of nodes in the second Layer (sop layer)
nodes2ndLayer=pf1_1*pf1_2+pf2_1*pf2_2+pf3_1*pf3_2+pf4_1*pf4_2;nodes2ndLayer;#4 kernels, in total (65*39)*4=10,140.

##(step3) Linearize the kernels/filters;
LinF1=as.vector(F1);LinF2=as.vector(F2);LinF3=as.vector(F3);LinF4=as.vector(F4);
LinF1;LinF2;LinF3;LinF4;

##(step4) set up a weight matrix W between lst and 2nd layers;W is regarded a 
## function W:nodes1stLayer X nodes2ndLaye-->R, i.e., W=bar40,000 X bar10,140->R
W=matrix(0,nodes1stLayer,nodes2ndLayer);W;dim(W);

k1=1;
for(j in 1:pf1_2)
{
  for(i in 1:pf1_1)
  {
  M=encodedIndex_f1[[i,j]];M;
  W[M,k1]=F1;
  k1=k1+1;
  }
} 

k2=pf1_2*pf1_1+1;
for(j in 1:pf2_2)
{
  for(i in 1:pf2_1)
  {
  M=encodedIndex_f2[[i,j]];M;
  W[M,k2]=F2;
  k2=k2+1;
  }
} 

k3=pf1_2*pf1_1+pf2_2*pf2_1+1;
for(j in 1:pf3_2)
{
  for(i in 1:pf3_1)
  {
  M=encodedIndex_f3[[i,j]];M;
  W[M,k3]=F3;
  k3=k3+1;
  }
} 

k4=pf1_2*pf1_1+pf2_2*pf2_1+pf3_2*pf3_1+1;
for(j in 1:pf4_2)
{
  for(i in 1:pf4_1)
  {
  M=encodedIndex_f4[[i,j]];M;
  W[M,k4]=F4;
  k4=k4+1;
  }
} 

##(step5) check the equivalence between standard cnn and our linearized model
sop=t(W)%*%vecInputTensor;length(sop);
MA1=matrix(sop[1:(pf1_2*pf1_1)],pf1_1,pf1_2);
imageShow(MA1);
MA2=matrix(sop[(pf1_2*pf1_1+1):(pf1_2*pf1_1+pf2_2*pf2_1)],pf2_1,pf2_2);
imageShow(MA2);
MA3=matrix(sop[(pf1_2*pf1_1+pf2_2*pf2_1+1):(pf1_2*pf1_1+pf2_2*pf2_1+pf3_2*pf3_1)],
pf3_1,pf3_2);
imageShow(MA3);
MA4=matrix(sop[(pf1_2*pf1_1+pf2_2*pf2_1+pf3_2*pf3_1+1):(pf1_2*pf1_1+pf2_2*pf2_1+pf3_2*pf3_1+pf4_2*pf4_1)],
pf4_1,pf4_2);
imageShow(MA4);


############ presentations  ###################
imageShow(Puppy);
imageShow(InputTensor);

imageShow(F1);
imageShow(F2);
imageShow(F3);
imageShow(F4);
image(F1);
image(F2);
image(F3);
image(F4);
imageShow(conv1);
imageShow(conv2);
imageShow(conv3);
imageShow(conv4);
imageShow(MA1);
imageShow(MA2);
imageShow(MA3);
imageShow(MA4);
F1;F2;F3;F4;
Puppy_G=Puppy[,,2];imageShow(Puppy_G);      
Puppy_B=Puppy[,,3];imageShow(Puppy_B);      

