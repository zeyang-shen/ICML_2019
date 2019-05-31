# ICML_2019
Codes and data for **"Interpreting Spacing Constraints of Transcription Factor Motifs from Convolutional Neural Networks"** published on Workshop on Computional Biology @ ICML 2019

Summary: We used CNNs to capture TF motifs and spacing between motifs at the same time. We designed our inputs to enable the interpretation directly related to TF motifs. The importance scores/saliency map incorporated the information of motif score as well as spacing and acted as a good indicator of TF binding.

Abstract: Transcription factors (TFs) bind to DNA regulatory sequences by recognizing their motifs. Convolutional neural networks (CNNs) have been applied in genomics to predict either TF binding or open chromatin induced by TF binding. Even though current models using CNNs can achieve high predictive performance, there is a lack of studies that relate the model interpretation beyond TF motifs to more complex regulatory codes such as spacing between TF motifs. We constructed a CNN model that incorporates spacing constraints of motifs, together with a training pipeline that connects model interpretation directly to motifs. The results showed that, by training with open chromatin regions, our CNN model can capture the co-localized motifs with spacing constraints, and the important motifs regarded by the model reflect actual TF binding. 

Motif files used in this project are originally from <https://github.com/jenhantao/abtba/tree/master/default_motifs> and also copied [here](https://github.com/zeyang-shen/ICML_2019/tree/master/default_motifs).

Additionally, an utility package ["szy"](https://github.com/zeyang-shen/ICML_2019/tree/master/szy) was made for this project to serve purposes of loading motif files, reading FASTA files, and converting DNA sequence to one-hot vectors.
