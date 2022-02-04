#Reports the depth of piRNA reads at each genome position within each piRNA clusters.

#USAGE: perl ReadCountsCoverage.pl [uniquely mapped piRNA read file] [piRNA cluster Annotation file]

perl ReadCountsCoverage.pl 21147-R1.piRNA.txt  piRNA_clusterAnnotation.txt
perl ReadCountsCoverage.pl 21147-R2.piRNA.txt  piRNA_clusterAnnotation.txt
perl ReadCountsCoverage.pl 21147-R3.piRNA.txt  piRNA_clusterAnnotation.txt

perl ReadCountsCoverage.pl 21183-R1.piRNA.txt  piRNA_clusterAnnotation.txt
perl ReadCountsCoverage.pl 21183-R2.piRNA.txt  piRNA_clusterAnnotation.txt
perl ReadCountsCoverage.pl 21183-R3.piRNA.txt  piRNA_clusterAnnotation.txt

perl ReadCountsCoverage.pl 21188-R1.piRNA.txt  piRNA_clusterAnnotation.txt
perl ReadCountsCoverage.pl 21188-R2.piRNA.txt  piRNA_clusterAnnotation.txt
perl ReadCountsCoverage.pl 21188-R3.piRNA.txt  piRNA_clusterAnnotation.txt

perl ReadCountsCoverage.pl 21213-R1.piRNA.txt  piRNA_clusterAnnotation.txt
perl ReadCountsCoverage.pl 21213-R2.piRNA.txt  piRNA_clusterAnnotation.txt
perl ReadCountsCoverage.pl 21213-R3.piRNA.txt  piRNA_clusterAnnotation.txt

perl ReadCountsCoverage.pl 21291-R1.piRNA.txt  piRNA_clusterAnnotation.txt
perl ReadCountsCoverage.pl 21291-R2.piRNA.txt  piRNA_clusterAnnotation.txt
perl ReadCountsCoverage.pl 21291-R3.piRNA.txt  piRNA_clusterAnnotation.txt

perl ReadCountsCoverage.pl 21346-R1.piRNA.txt  piRNA_clusterAnnotation.txt
perl ReadCountsCoverage.pl 21346-R2.piRNA.txt  piRNA_clusterAnnotation.txt
perl ReadCountsCoverage.pl 21346-R3.piRNA.txt  piRNA_clusterAnnotation.txt


#Reports the depth of piRNA reads at piRNA cluster intervals.

#USAGE: perl piRNA_clusterCoverage.pl [ReadCountCoverage read file] 

perl piRNA_clusterCoverage.pl 21147-R1_ReadCountCoverage.txt
perl piRNA_clusterCoverage.pl 21147-R2_ReadCountCoverage.txt
perl piRNA_clusterCoverage.pl 21147-R3_ReadCountCoverage.txt

perl piRNA_clusterCoverage.pl 21183-R1_ReadCountCoverage.txt
perl piRNA_clusterCoverage.pl 21183-R2_ReadCountCoverage.txt
perl piRNA_clusterCoverage.pl 21183-R3_ReadCountCoverage.txt

perl piRNA_clusterCoverage.pl 21188-R1_ReadCountCoverage.txt
perl piRNA_clusterCoverage.pl 21188-R2_ReadCountCoverage.txt
perl piRNA_clusterCoverage.pl 21188-R3_ReadCountCoverage.txt

perl piRNA_clusterCoverage.pl 21213-R1_ReadCountCoverage.txt
perl piRNA_clusterCoverage.pl 21213-R2_ReadCountCoverage.txt
perl piRNA_clusterCoverage.pl 21213-R3_ReadCountCoverage.txt

perl piRNA_clusterCoverage.pl 21291-R1_ReadCountCoverage.txt
perl piRNA_clusterCoverage.pl 21291-R2_ReadCountCoverage.txt
perl piRNA_clusterCoverage.pl 21291-R3_ReadCountCoverage.txt

perl piRNA_clusterCoverage.pl 21346-R1_ReadCountCoverage.txt
perl piRNA_clusterCoverage.pl 21346-R2_ReadCountCoverage.txt
perl piRNA_clusterCoverage.pl 21346-R3_ReadCountCoverage.txt
