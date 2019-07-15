"Mitra Menon" />


<title>INTROGRESS and Fold Enrichment analysis</title>

<!-- rnb-text-begin -->
<div id="background" class="section level1">
<h1>BACKGROUND</h1>
</div>
<div id="overall-goal-is-to-identify-adaptively-introgressed-loci-and-determine-the-major-environmental-axis-associated-with-adaptive-introgression" class="section level1">
<h1><em>Overall goal is to identify adaptively introgressed loci and determine the major environmental axis associated with adaptive introgression</em></h1>
<p>Here I utilize the regression approach implemented within <code>INTROGRESS</code> to first identify loci exhibiting exceptional patterns of introgression. These will then be subjected to a fold enrichment approach to determine if loci associated with a given environmental variable are more likely to be candidates for adaptive introgression than those associated with other environmental variables. We will focussing primarily on environmental variables that are strongly distinguished between the two hybridizing species, or those driving the niches of the two species (if you have prior information on that).</p>
<!-- rnb-text-end -->
<!-- rnb-chunk-begin -->
<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxubGlicmFyeShkYXRhLnRhYmxlKVxubGlicmFyeShpbnRyb2dyZXNzKVxuYGBgIn0= -->
<pre class="r"><code>library(data.table)
library(introgress)</code></pre>
<!-- rnb-source-end -->
<!-- rnb-chunk-end -->
<!-- rnb-text-begin -->
</div>
<div id="step-1-prep-data-and-run-introgress-for-large-number-of-loci-this-needs-to-be-run-on-a-computing-cluster" class="section level1">
<h1>STEP 1: Prep data and run introgress (<em>For large number of loci this needs to be run on a computing cluster</em>)</h1>
<p>Load data</p>
<!-- rnb-text-end -->
<!-- rnb-chunk-begin -->
<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuR2RhdGE8LWZyZWFkKFwiLi4vc25wRGF0YS9JbGx1bWluYU1peF8zNDg5L3Jhd0RhdGEvbWlub3IwMTIudHh0XCIsc2VwPVwiXFx0XCIsaGVhZGVyPVQsZGF0YS50YWJsZT1GKVxuUG9wczwtcmVhZC50YWJsZShcIi4uL3NucERhdGEvUG9wSW5kLnR4dFwiLGhlYWRlcj1ULHNlcD1cIlxcdFwiLHN0cmluZ3NBc0ZhY3RvcnMgPSBGQUxTRSlcbnBhcmVudGFsczwtcmVhZC50YWJsZShcIi4uL3NucERhdGEvc3VtbWFyeVN0YXRzXzM0ODkvTUFGX3BhcmVudGFscy50eHRcIixoZWFkZXI9VCxzZXA9XCJcXHRcIilcbmBgYCJ9 -->
<pre class="r"><code>Gdata&lt;-fread(&quot;../snpData/IlluminaMix_3489/rawData/minor012.txt&quot;,sep=&quot;\t&quot;,header=T,data.table=F)
Pops&lt;-read.table(&quot;../snpData/PopInd.txt&quot;,header=T,sep=&quot;\t&quot;,stringsAsFactors = FALSE)
parentals&lt;-read.table(&quot;../snpData/summaryStats_3489/MAF_parentals.txt&quot;,header=T,sep=&quot;\t&quot;)</code></pre>
<!-- rnb-source-end -->
<!-- rnb-chunk-end -->
<!-- rnb-text-begin -->
<p>Due to low overall differentiation between the two species I am using, I will use a very low allele frequency difference cutoff to obtain the input loci. As states in the paper, the approach should not be biased by this. Maximising the number of loci will also help the assesment of fold enrichment.</p>
<!-- rnb-text-end -->
<!-- rnb-chunk-begin -->
<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuYWxsZWxlRnJlcTwtY2JpbmQoZGlmZj1hYnMocGFyZW50YWxzJExQX3B1cmUtcGFyZW50YWxzJFNXV1BfcHVyZSkscGFyZW50YWxzKVxuYWxsZWxlRnJlcTwtYWxsZWxlRnJlcVthbGxlbGVGcmVxJGRpZmY+MCwgXSBcblxuR2RhdGE8LUdkYXRhWyAsY29sbmFtZXMoR2RhdGEpJWluJWFsbGVsZUZyZXEkbG9jaV1cbkdkYXRhPC1jYmluZChQb3BzLEdkYXRhKVxuR2RhdGE8LUdkYXRhWyEoR2RhdGEkUG9wPT1cIk1FU0MxTFwiKSwgXVxuUG9wczwtUG9wc1shKFBvcHMkUG9wPT1cIk1FU0MxTFwiKSwgXVxuXG5sb2NpPC1hcy5tYXRyaXgoR2RhdGFbICwtYygxOjUpXSlcbmxvY2lbaXMubmEobG9jaSldPC1cIk5BL05BXCJcbmxvY2lbbG9jaT09MV08LVwiQS9EXCJcbmxvY2lbbG9jaT09MF08LVwiQS9BXCJcbmxvY2lbbG9jaT09Ml08LVwiRC9EXCJcbmBgYCJ9 -->
<pre class="r"><code>alleleFreq&lt;-cbind(diff=abs(parentals$LP_pure-parentals$SWWP_pure),parentals)
alleleFreq&lt;-alleleFreq[alleleFreq$diff&gt;0, ] 

Gdata&lt;-Gdata[ ,colnames(Gdata)%in%alleleFreq$loci]
Gdata&lt;-cbind(Pops,Gdata)
Gdata&lt;-Gdata[!(Gdata$Pop==&quot;MESC1L&quot;), ]
Pops&lt;-Pops[!(Pops$Pop==&quot;MESC1L&quot;), ]

loci&lt;-as.matrix(Gdata[ ,-c(1:5)])
loci[is.na(loci)]&lt;-&quot;NA/NA&quot;
loci[loci==1]&lt;-&quot;A/D&quot;
loci[loci==0]&lt;-&quot;A/A&quot;
loci[loci==2]&lt;-&quot;D/D&quot;</code></pre>
<!-- rnb-source-end -->
<!-- rnb-chunk-end -->
<!-- rnb-text-begin -->
<p>Define groups</p>
<!-- rnb-text-end -->
<!-- rnb-chunk-begin -->
<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuTFA8LWMoXCJCQ0tcIixcIkJVXCIsXCJDSFwiLFwiQ1BcIixcIkhPXCIsXCJKRU5cIixcIkpWXCIsXCJMUFwiLFwiTUVMXCIsXCJNU0xcIixcIlZDXCIsXCJDT0xPMkxcIixcIkNPTE8zSFwiKVxuU1dXUDwtYyhcIkRBVjJMXCIsXCJBTFRcIixcIkhVQTJIXCIsXCJIVUEyTFwiLFwiU0FOUklUMUxcIixcIkJBTFwiLFwiQkFTXCIsXCJDSElcIixcIkVQQVwiLFwiRVJNXCIsXCJHVUFcIixcIkdVSVwiLFwiU0FNXCIsXCJUT01cIixcIlZFUlwiLFwiWUFRXCIsXCJDSEkxSFwiLFwiQ0hJMUxcIixcIkNISTNcIixcIkRBVjFMXCIsXCJEQVYxSFwiLFwiREFWMkhcIilcbnBvcDk4PC1yZWFkLnRhYmxlKFwiLi4vc25wRGF0YS9zdW1tYXJ5U3RhdHNfMzQ4OS9iYXllbnZPdXQvOThQb3BJRHNcIilcblxuR2RhdGE8LWNiaW5kLmRhdGEuZnJhbWUoUG9wcyxsb2NpKVxuXG5QMV9zd3dwPC1HZGF0YVtHZGF0YSRQb3AlaW4lU1dXUCwgXVxuUDJfbHA8LUdkYXRhW0dkYXRhJFBvcCVpbiVMUCwgXVxubWl4ZWQ8LUdkYXRhW0dkYXRhJFBvcCVpbiVwb3A5OCRWMSwgXVxuYGBgIn0= -->
<pre class="r"><code>LP&lt;-c(&quot;BCK&quot;,&quot;BU&quot;,&quot;CH&quot;,&quot;CP&quot;,&quot;HO&quot;,&quot;JEN&quot;,&quot;JV&quot;,&quot;LP&quot;,&quot;MEL&quot;,&quot;MSL&quot;,&quot;VC&quot;,&quot;COLO2L&quot;,&quot;COLO3H&quot;)
SWWP&lt;-c(&quot;DAV2L&quot;,&quot;ALT&quot;,&quot;HUA2H&quot;,&quot;HUA2L&quot;,&quot;SANRIT1L&quot;,&quot;BAL&quot;,&quot;BAS&quot;,&quot;CHI&quot;,&quot;EPA&quot;,&quot;ERM&quot;,&quot;GUA&quot;,&quot;GUI&quot;,&quot;SAM&quot;,&quot;TOM&quot;,&quot;VER&quot;,&quot;YAQ&quot;,&quot;CHI1H&quot;,&quot;CHI1L&quot;,&quot;CHI3&quot;,&quot;DAV1L&quot;,&quot;DAV1H&quot;,&quot;DAV2H&quot;)
pop98&lt;-read.table(&quot;../snpData/summaryStats_3489/bayenvOut/98PopIDs&quot;)

Gdata&lt;-cbind.data.frame(Pops,loci)

P1_swwp&lt;-Gdata[Gdata$Pop%in%SWWP, ]
P2_lp&lt;-Gdata[Gdata$Pop%in%LP, ]
mixed&lt;-Gdata[Gdata$Pop%in%pop98$V1, ]</code></pre>
<!-- rnb-source-end -->
<!-- rnb-chunk-end -->
<!-- rnb-text-begin -->
<p>Data prep for running INTROGRESS</p>
<!-- rnb-text-end -->
<!-- rnb-chunk-begin -->
<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuUDFfc3d3cDwtUDFfc3d3cFsgLC1jKDE6NSldXG5QMl9scDwtUDJfbHBbICwtYygxOjUpXVxubWl4ZWQ8LWNiaW5kKFBvcD1taXhlZCRQb3AsSW5kPW1peGVkJFBvcEluZCxtaXhlZFsgLC1jKDE6NSldKVxuXG5QMV90cmFuc3Bvc2U8LXQoUDFfc3d3cClcblAyX3RyYW5zcG9zZTwtdChQMl9scClcbm1peGVkX3RyYW5zcG9zZTwtdChtaXhlZClcblxuXG5Mb2N1czwtY2JpbmQobG9jdXM9YWxsZWxlRnJlcSRsb2NpLHR5cGU9cmVwKFwiQ1wiLG5yb3coYWxsZWxlRnJlcSkpKVxuXG5hZG1peENvdW50PC1wcmVwYXJlLmRhdGEoYWRtaXguZ2VuPW1peGVkX3RyYW5zcG9zZSxsb2NpLmRhdGE9TG9jdXMscGFyZW50YWwxPVAxX3RyYW5zcG9zZSxwYXJlbnRhbDI9UDJfdHJhbnNwb3NlLHBvcC5pZD1UUlVFLGluZC5pZD1UUlVFLGZpeGVkPUZBTFNFLHNlcC5yb3dzID0gRkFMU0Usc2VwLmNvbHVtbnMgPSBGQUxTRSlcbmBgYCJ9 -->
<pre class="r"><code>P1_swwp&lt;-P1_swwp[ ,-c(1:5)]
P2_lp&lt;-P2_lp[ ,-c(1:5)]
mixed&lt;-cbind(Pop=mixed$Pop,Ind=mixed$PopInd,mixed[ ,-c(1:5)])

P1_transpose&lt;-t(P1_swwp)
P2_transpose&lt;-t(P2_lp)
mixed_transpose&lt;-t(mixed)


Locus&lt;-cbind(locus=alleleFreq$loci,type=rep(&quot;C&quot;,nrow(alleleFreq)))

admixCount&lt;-prepare.data(admix.gen=mixed_transpose,loci.data=Locus,parental1=P1_transpose,parental2=P2_transpose,pop.id=TRUE,ind.id=TRUE,fixed=FALSE,sep.rows = FALSE,sep.columns = FALSE)</code></pre>
<!-- rnb-source-end -->
<!-- rnb-chunk-end -->
<!-- rnb-text-begin -->
<p>Estimate hybrid index and interspecfic het, this is important to determine outlier loci with respect to genome wide ancestry. (<em>takes a long time and so save the output after it is run</em>)</p>
<!-- rnb-text-end -->
<!-- rnb-chunk-begin -->
<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuSEluZGV4PC1lc3QuaChpbnRyb2dyZXNzLmRhdGEgPSBhZG1peENvdW50LGxvY2kuZGF0YSA9IExvY3VzLGluZC50b3VzZSA9IE5VTEwsZml4ZWQgPSBGQUxTRSlcbmhlYWQoSEluZGV4KVxuI3dyaXRlLnRhYmxlKGNiaW5kKEhJbmRleCxtaXhlZFsgLDE6Ml0pLGZpbGU9XCJIeWJyaWRJbmRleDk4UG9wcy50eHRcIixyb3cubmFtZXM9RixxdW90ZT1GLHNlcD1cIlxcdFwiKVxuXG5pbnQuaGV0PC1jYWxjLmludGVyc3AuaGV0KGludHJvZ3Jlc3MuZGF0YT1hZG1peENvdW50KVxuYGBgIn0= -->
<pre class="r"><code>HIndex&lt;-est.h(introgress.data = admixCount,loci.data = Locus,ind.touse = NULL,fixed = FALSE)
head(HIndex)
#write.table(cbind(HIndex,mixed[ ,1:2]),file=&quot;HybridIndex98Pops.txt&quot;,row.names=F,quote=F,sep=&quot;\t&quot;)

int.het&lt;-calc.intersp.het(introgress.data=admixCount)</code></pre>
<!-- rnb-source-end -->
<!-- rnb-chunk-end -->
<!-- rnb-text-begin -->
<p>Now, conduct genomic cline analysis using the parametric approach since the SNPs used don’t exhibit fixed differences (<em>Takes a long time and needs to be run on the cluster</em>)</p>
<!-- rnb-text-end -->
<!-- rnb-chunk-begin -->
<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuR2NsaW5lc19wYXJhPC1nZW5vbWljLmNsaW5lcyhpbnRyb2dyZXNzLmRhdGEgPSBhZG1peENvdW50LGhpLmluZGV4ID0gSEluZGV4LGxvY2kuZGF0YSA9IExvY3VzLFxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIG1ldGhvZD1cInBhcmFtZXRyaWNcIixzaWcudGVzdD1UUlVFLGxvY2kudG91c2U9TlVMTCxpbmQudG91c2U9TlVMTClcbmBgYCJ9 -->
<pre class="r"><code>Gclines_para&lt;-genomic.clines(introgress.data = admixCount,hi.index = HIndex,loci.data = Locus,
                               method=&quot;parametric&quot;,sig.test=TRUE,loci.touse=NULL,ind.touse=NULL)</code></pre>
<!-- rnb-source-end -->
<!-- rnb-chunk-end -->
<!-- rnb-text-begin -->
<p>Write output files</p>
<!-- rnb-text-end -->
<!-- rnb-chunk-begin -->
<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxud3JpdGUudGFibGUoR2NsaW5lc19wYXJhJFN1bW1hcnkuZGF0YSwgZmlsZT1cIlN1bW1fcGFyYTEwMDAudHh0XCIscXVvdGU9RkFMU0UsIHNlcD1cIlxcdFwiKVxud3JpdGUudGFibGUoR2NsaW5lc19wYXJhJEZpdHRlZC5BQSxmaWxlPVwiSG9tb1AxX3BhcmExMDAwLnR4dFwiLHNlcD1cIlxcdFwiLHF1b3RlPUYpXG53cml0ZS50YWJsZShHY2xpbmVzX3BhcmEkRml0dGVkLmFhLGZpbGU9XCJIb21vUDJfcGFyYTEwMDAudHh0XCIsc2VwPVwiXFx0XCIscXVvdGU9RilcbndyaXRlLnRhYmxlKEdjbGluZXNfcGFyYSROZXV0cmFsLkFBLGZpbGU9XCJIb21vUDFfQ0lwYXJhMTAwMC50eHRcIixzZXA9XCJcXHRcIixxdW90ZT1GKVxud3JpdGUudGFibGUoR2NsaW5lc19wYXJhJE5ldXRyYWwuYWEsZmlsZT1cIkhvbW9QMl9DSXBhcmExMDAwLnR4dFwiLHNlcD1cIlxcdFwiLHF1b3RlPUYpXG53cml0ZS50YWJsZShHY2xpbmVzX3BhcmEkUXVhbnRpbGVzLGZpbGU9XCJRdWFudGlsZV9wYXJhMTAwMC50eHRcIixzZXA9XCJcXHRcIixxdW90ZT1GKVxuYGBgIn0= -->
<pre class="r"><code>write.table(Gclines_para$Summary.data, file=&quot;Summ_para1000.txt&quot;,quote=FALSE, sep=&quot;\t&quot;)
write.table(Gclines_para$Fitted.AA,file=&quot;HomoP1_para1000.txt&quot;,sep=&quot;\t&quot;,quote=F)
write.table(Gclines_para$Fitted.aa,file=&quot;HomoP2_para1000.txt&quot;,sep=&quot;\t&quot;,quote=F)
write.table(Gclines_para$Neutral.AA,file=&quot;HomoP1_CIpara1000.txt&quot;,sep=&quot;\t&quot;,quote=F)
write.table(Gclines_para$Neutral.aa,file=&quot;HomoP2_CIpara1000.txt&quot;,sep=&quot;\t&quot;,quote=F)
write.table(Gclines_para$Quantiles,file=&quot;Quantile_para1000.txt&quot;,sep=&quot;\t&quot;,quote=F)</code></pre>
<!-- rnb-source-end -->
<!-- rnb-chunk-end -->
<!-- rnb-text-begin -->
<p>The parametric approach is prone to false positives, specifically due to the large number of tests done, we will conduct a P.val correction and use only loci that pass p val threshold.</p>
<!-- rnb-text-end -->
<!-- rnb-chunk-begin -->
<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuU3VtbTwtcmVhZC50YWJsZShcIi4uL3NucERhdGEvc3VtbWFyeVN0YXRzXzM0ODkvaW50cm9ncmVzcy9hbGwxMDAwcmVwcy9TdW1tX3BhcmExMDAwLnR4dFwiLGhlYWRlcj1ULHNlcD1cIlxcdFwiKVxuXG5udW08LTAuMDUvbnJvdyhTdW1tKSAjQm9uZmVyb25uaSBjb3JyZWN0aW9uXG5TdW1tQkg8LVN1bW1bU3VtbSRQLnZhbHVlPG51bSwgXVxuYGBgIn0= -->
<pre class="r"><code>Summ&lt;-read.table(&quot;../snpData/summaryStats_3489/introgress/all1000reps/Summ_para1000.txt&quot;,header=T,sep=&quot;\t&quot;)

num&lt;-0.05/nrow(Summ) #Bonferonni correction
SummBH&lt;-Summ[Summ$P.value&lt;num, ]</code></pre>
<!-- rnb-source-end -->
<!-- rnb-chunk-end -->
<!-- rnb-text-begin -->
<p>Across all individuals we will only retain loci that have passed the Bonferonni correction. Since we have estimates for each individual, we will then classify the loci within an individual as having excess ancestry from LP using the CI estimates.</p>
<!-- rnb-text-end -->
<!-- rnb-chunk-begin -->
<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuTFBjaTwtZnJlYWQoXCIuLi9zbnBEYXRhL3N1bW1hcnlTdGF0c18zNDg5L2ludHJvZ3Jlc3MvYWxsTG9jaS9Ib21vUDJfQ0lwYXJhMTAwMC50eHRcIixzZXA9XCJcXHRcIixkYXRhLnRhYmxlPUYpXG5MUG91dDwtZnJlYWQoXCIuLi9zbnBEYXRhL3N1bW1hcnlTdGF0c18zNDg5L2ludHJvZ3Jlc3MvYWxsTG9jaS9Ib21vUDJfcGFyYTEwMDAudHh0XCIsc2VwPVwiXFx0XCIsZGF0YS50YWJsZT1GKVxuXG5MUF91cHBlcjwtTFBjaVsgLDE6bmNvbChMUG91dCldXG5cblxuRXhjZXNzTFA8LW1hdHJpeChucm93ID0gbnJvdyhMUG91dCksbmNvbD05NTApXG5cbmZvciAoYyBpbiAyOm5jb2woTFBvdXQpKXtcbiAgbjwtYy0xXG4gIGZvciAociBpbiAxOm5yb3coTFBvdXQpKXtcbiAgICBFeGNlc3NMUFtyLG5dPC1pZmVsc2UoTFBvdXRbcixjXT5MUF91cHBlcltyLGNdLDEsMClcbiAgfVxufVxuXG5yb3duYW1lcyhFeGNlc3NMUCk8LUxQb3V0JFYxXG5jb2xuYW1lcyhFeGNlc3NMUCk8LWNvbG5hbWVzKExQb3V0KVsyOm5jb2woTFBvdXQpXVxuRXhjZXNzTFBfb3V0PC1FeGNlc3NMUFtyb3duYW1lcyhFeGNlc3NMUCklaW4lU3VtbUJIJGxvY3VzLCBdXG5cbmBgYCJ9 -->
<pre class="r"><code>LPci&lt;-fread(&quot;../snpData/summaryStats_3489/introgress/allLoci/HomoP2_CIpara1000.txt&quot;,sep=&quot;\t&quot;,data.table=F)
LPout&lt;-fread(&quot;../snpData/summaryStats_3489/introgress/allLoci/HomoP2_para1000.txt&quot;,sep=&quot;\t&quot;,data.table=F)

LP_upper&lt;-LPci[ ,1:ncol(LPout)]


ExcessLP&lt;-matrix(nrow = nrow(LPout),ncol=950)

for (c in 2:ncol(LPout)){
  n&lt;-c-1
  for (r in 1:nrow(LPout)){
    ExcessLP[r,n]&lt;-ifelse(LPout[r,c]&gt;LP_upper[r,c],1,0)
  }
}

rownames(ExcessLP)&lt;-LPout$V1
colnames(ExcessLP)&lt;-colnames(LPout)[2:ncol(LPout)]
ExcessLP_out&lt;-ExcessLP[rownames(ExcessLP)%in%SummBH$locus, ]
</code></pre>
<!-- rnb-source-end -->
<!-- rnb-chunk-end -->
<!-- rnb-text-begin -->
<p>After assigning introgressed loci as 1 and non-introgressed as 0, we can add them across all individuals to determine what proportion of times it was introgressed across the hybrid zone. We can then use a cutoff (mostly arbitary) to classify a loci as significantly introgressed. This arbitary cutoff for here is set as 0.2 (20% of individiuals) and was assessed through a series of cutoffs being subjected to downstream analysis, all of which gave similar results.</p>
<!-- rnb-text-end -->
<!-- rnb-chunk-begin -->
<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuRXhjZXNzTFBfb3V0PC1hcHBseShFeGNlc3NMUCwxLGZ1bmN0aW9uKFgpIHJldHVybihzdW0oWCkvbGVuZ3RoKFgpKSlcbkV4Y2Vzc0xQX291dDI8LUV4Y2Vzc0xQX291dFtFeGNlc3NMUF9vdXQ+MC4yMF1cbmBgYCJ9 -->
<pre class="r"><code>ExcessLP_out&lt;-apply(ExcessLP,1,function(X) return(sum(X)/length(X)))
ExcessLP_out2&lt;-ExcessLP_out[ExcessLP_out&gt;0.20]</code></pre>
<!-- rnb-source-end -->
<!-- rnb-chunk-end -->
<!-- rnb-text-begin -->
</div>
<div id="step-2-determining-whether-the-introgressed-loci-are-adaptive." class="section level1">
<h1>STEP 2: Determining whether the introgressed loci are adaptive.</h1>
<p>Here we will utilize the outliers identified through GEA (bayenv in my case) and intersect them with the introgressed loci. Since the number of outliers loci vary across environmental variables, we will utilize a permutation approach to assess the significance of fold enrichment.</p>
<p>Load GEA outliers SNPs</p>
<!-- rnb-text-end -->
<!-- rnb-chunk-begin -->
<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuY2xpbTwtbGlzdC5maWxlcyhcIi4uL3NucERhdGEvc3VtbWFyeVN0YXRzXzM0ODkvYmF5ZW52T3V0L2NsaW0vY29udmVyZ2VuY2UvXCIsZnVsbC5uYW1lcyA9IFQscGF0dGVybj1cIk92ZXJsYXBcIilcblNvaWw8LWxpc3QuZmlsZXMoXCIuLi9zbnBEYXRhL3N1bW1hcnlTdGF0c18zNDg5L2JheWVudk91dC9Tb2lsL2NvbnZlcmdlbmNlL1wiLGZ1bGwubmFtZXMgPSBULHBhdHRlcm49XCJPdmVybGFwXCIpXG5vdXRsaWVyczwtYyhjbGltLFNvaWwpXG5cblxuZmlsZXNCRjwtdmVjdG9yKFwibGlzdFwiLGxlbmd0aChvdXRsaWVycykpXG5cbmZvciAoZiBpbiAxOmxlbmd0aChmaWxlc0JGKSl7XG4gIFxuICBkZjwtcmVhZC50YWJsZShvdXRsaWVyc1tmXSxoZWFkZXI9VCxzZXA9XCJcXHRcIilcbiAgZGY8LWRmWyAsMTo1XVxuICBmaWxlc0JGW1tmXV08LWRmXG59XG5cbmZpbGVzPC1saXN0LmZpbGVzKFwiLi4vc25wRGF0YS9zdW1tYXJ5U3RhdHNfMzQ4OS9iYXllbnZPdXQvY2xpbS9jb252ZXJnZW5jZS9cIixwYXR0ZXJuPVwiT3ZlcmxhcFwiKVxuZmlsZXM8LWMoZmlsZXMsbGlzdC5maWxlcyhcIi4uL3NucERhdGEvc3VtbWFyeVN0YXRzXzM0ODkvYmF5ZW52T3V0L1NvaWwvY29udmVyZ2VuY2UvXCIscGF0dGVybj1cIk92ZXJsYXBcIikpXG52YXJJRDwtc2FwcGx5KHN0cnNwbGl0KGZpbGVzLFwiT3ZlcmxhcDFfM2NoYWluQ29udmdcIiksXCJbXCIsMilcbnZhcklEPC1zYXBwbHkoc3Ryc3BsaXQodmFySUQsXCIudHh0XCIpLFwiW1wiLDEpXG5cblxubmFtZXMoZmlsZXNCRik8LXZhcklEXG5yZW1vdmU8LWMoXCJFbGV2YXRpb25cIiAsXCJMYXRpdHVkZVwiLFwiTG9uZ2l0dWRlXCIpXG5vdXRsaWVyczwtZmlsZXNCRlshKG5hbWVzKGZpbGVzQkYpJWluJXJlbW92ZSldXG5gYGAifQ== -->
<pre class="r"><code>clim&lt;-list.files(&quot;../snpData/summaryStats_3489/bayenvOut/clim/convergence/&quot;,full.names = T,pattern=&quot;Overlap&quot;)
Soil&lt;-list.files(&quot;../snpData/summaryStats_3489/bayenvOut/Soil/convergence/&quot;,full.names = T,pattern=&quot;Overlap&quot;)
outliers&lt;-c(clim,Soil)


filesBF&lt;-vector(&quot;list&quot;,length(outliers))

for (f in 1:length(filesBF)){
  
  df&lt;-read.table(outliers[f],header=T,sep=&quot;\t&quot;)
  df&lt;-df[ ,1:5]
  filesBF[[f]]&lt;-df
}

files&lt;-list.files(&quot;../snpData/summaryStats_3489/bayenvOut/clim/convergence/&quot;,pattern=&quot;Overlap&quot;)
files&lt;-c(files,list.files(&quot;../snpData/summaryStats_3489/bayenvOut/Soil/convergence/&quot;,pattern=&quot;Overlap&quot;))
varID&lt;-sapply(strsplit(files,&quot;Overlap1_3chainConvg&quot;),&quot;[&quot;,2)
varID&lt;-sapply(strsplit(varID,&quot;.txt&quot;),&quot;[&quot;,1)


names(filesBF)&lt;-varID
remove&lt;-c(&quot;Elevation&quot; ,&quot;Latitude&quot;,&quot;Longitude&quot;)
outliers&lt;-filesBF[!(names(filesBF)%in%remove)]</code></pre>
<!-- rnb-source-end -->
<!-- rnb-chunk-end -->
<!-- rnb-text-begin -->
<p>Some simple manipulation to reatain only the SNPs used in <code>INTROGRESS</code></p>
<!-- rnb-text-end -->
<!-- rnb-chunk-begin -->
<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxub3V0bGllcnMyPC1sYXBwbHkob3V0bGllcnMsZnVuY3Rpb24oZGYpIGRmW2RmJGxvY2klaW4lTFBvdXQkVjEsIF0pICNOVUxMIFNFVFxuQkZfZ0NsaW5lczwtbGFwcGx5KG91dGxpZXJzMixmdW5jdGlvbihkZikgcmV0dXJuKGRmW2RmJGxvY2klaW4lbmFtZXMoRXhjZXNzTFBfb3V0MiksIF0pKSAjU0hBUkVEIFNFVFxuQkZfZ0NsaW5lc0lEPC11bmxpc3QobGFwcGx5KEJGX2dDbGluZXMsZnVuY3Rpb24oZGYpIG5yb3coZGYpKSlcbmBgYCJ9 -->
<pre class="r"><code>outliers2&lt;-lapply(outliers,function(df) df[df$loci%in%LPout$V1, ]) #NULL SET
BF_gClines&lt;-lapply(outliers2,function(df) return(df[df$loci%in%names(ExcessLP_out2), ])) #SHARED SET
BF_gClinesID&lt;-unlist(lapply(BF_gClines,function(df) nrow(df)))</code></pre>
<!-- rnb-source-end -->
<!-- rnb-chunk-end -->
<!-- rnb-text-begin -->
<p>Estimating fold change for each environmental variable</p>
<!-- rnb-text-end -->
<!-- rnb-chunk-begin -->
<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuTnI8LUJGX2dDbGluZXNJRC91bmxpc3QobGFwcGx5KG91dGxpZXJzMixmdW5jdGlvbihkZikgbnJvdyhkZikpKVxuRG48LWxlbmd0aChFeGNlc3NMUF9vdXQyKS9ucm93KExQb3V0KVxuXG5GQzwtTnIvRG5cbmBgYCJ9 -->
<pre class="r"><code>Nr&lt;-BF_gClinesID/unlist(lapply(outliers2,function(df) nrow(df)))
Dn&lt;-length(ExcessLP_out2)/nrow(LPout)

FC&lt;-Nr/Dn</code></pre>
<!-- rnb-source-end -->
<!-- rnb-chunk-end -->
<!-- rnb-text-begin -->
<p><em>Function for randomization to obtain null distribution of FC</em></p>
<p>Takes 4 inputs:</p>
<p><code>df</code> is a dataframe of bayenv outliers, where nrow = total outliers <code>loci</code> is a character vector of IDs of all snps used in introgress analysis <code>Int</code> total number of loci identified as significantly introgressed from LP in STEP 1 <code>R</code> number of bootstrap replicates to run</p>
<!-- rnb-text-end -->
<!-- rnb-chunk-begin -->
<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxucmFuZEZDPC1mdW5jdGlvbihkZixsb2NpLEludCxSKXtcblxuICBib290PC1OVUxMXG4gIGZvcihpIGluIDE6Uil7XG4gICAgXG4gICAgQkZyYW5kPC1zYW1wbGUobG9jaSxucm93KGRmKSxyZXBsYWNlID0gRilcbiAgICBQRnJhbmQ8LXNhbXBsZShsb2NpLEludCxyZXBsYWNlID0gRilcbiAgICBzaGFyZWQ8LWxlbmd0aChQRnJhbmRbUEZyYW5kJWluJUJGcmFuZF0pXG4gICAgXG4gICAgYm9vdFtpXTwtKHNoYXJlZC9sZW5ndGgoQkZyYW5kKSkvKGxlbmd0aChQRnJhbmQpL2xlbmd0aChsb2NpKSlcbiAgfVxuXG4gIHJldHVybihib290KVxufVxuXG5gYGAifQ== -->
<pre class="r"><code>randFC&lt;-function(df,loci,Int,R){

  boot&lt;-NULL
  for(i in 1:R){
    
    BFrand&lt;-sample(loci,nrow(df),replace = F)
    PFrand&lt;-sample(loci,Int,replace = F)
    shared&lt;-length(PFrand[PFrand%in%BFrand])
    
    boot[i]&lt;-(shared/length(BFrand))/(length(PFrand)/length(loci))
  }

  return(boot)
}
</code></pre>
<!-- rnb-source-end -->
<!-- rnb-chunk-end -->
<!-- rnb-text-begin -->
<p>The actual test</p>
<!-- rnb-text-end -->
<!-- rnb-chunk-begin -->
<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuSUQ8LUxQb3V0JFYxXG5JbnQ8LWxlbmd0aChFeGNlc3NMUF9vdXQyKVxuXG5cbkZDcmFuZE91dDwtbGFwcGx5KG91dGxpZXJzMixmdW5jdGlvbihkZikgcmV0dXJuKHJhbmRGQyhkZixJRCxJbnQsMTAwMDApKSlcbmBgYCJ9 -->
<pre class="r"><code>ID&lt;-LPout$V1
Int&lt;-length(ExcessLP_out2)


FCrandOut&lt;-lapply(outliers2,function(df) return(randFC(df,ID,Int,10000)))</code></pre>
<!-- rnb-source-end -->
<!-- rnb-chunk-end -->
<!-- rnb-text-begin -->
<p>Check which variables have observed FC outside the 99th percentile of null</p>
<!-- rnb-text-end -->
<!-- rnb-chunk-begin -->
<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuXG5mb3IgKGkgaW4gMTpsZW5ndGgoRkNyYW5kT3V0KSl7XG4gIEVudjwtbmFtZXMoRkNyYW5kT3V0KVtpXVxuICBlbXA8LUZDcmFuZE91dFtbRW52XV1cbiAgaWYgKHF1YW50aWxlKGVtcCxwcm9icyA9IGMoMC45OTkpKTw9RkNbW0Vudl1dKXtcbiAgICBjYXQoXCJ2YXJpYWJsZVwiLEVudixcImlzIHNpZ25pZmljYW50bHkgZW5yaWNoZWRcIixcIlxcblwiKVxuICB9XG4gIFxufVxuXG5gYGAifQ== -->
<pre class="r"><code>
for (i in 1:length(FCrandOut)){
  Env&lt;-names(FCrandOut)[i]
  emp&lt;-FCrandOut[[Env]]
  if (quantile(emp,probs = c(0.999))&lt;=FC[[Env]]){
    cat(&quot;variable&quot;,Env,&quot;is significantly enriched&quot;,&quot;\n&quot;)
  }
  
}
</code></pre>
<!-- rnb-source-end -->
<!-- rnb-chunk-end -->
<!-- rnb-text-begin -->
<!-- rnb-text-end -->
</div>

<div id="rmd-source-code">

</div>

<script>

// add bootstrap table styles to pandoc tables
function bootstrapStylePandocTables() {
  $('tr.header').parent('thead').parent('table').addClass('table table-condensed');
}
$(document).ready(function () {
  bootstrapStylePandocTables();
});

$(document).ready(function () {
  $('.knitsql-table').addClass('kable-table');
  var container = $('.kable-table');
  container.each(function() {

    // move the caption out of the table
    var table = $(this).children('table');
    var caption = table.children('caption').detach();
    caption.insertBefore($(this)).css('display', 'inherit');
  });
});

</script>

<!-- dynamically load mathjax for compatibility with self-contained -->
<script>
  (function () {
    var script = document.createElement("script");
    script.type = "text/javascript";
    script.src  = "https://mathjax.rstudio.com/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML";
    document.getElementsByTagName("head")[0].appendChild(script);
  })();
</script>

</body>
</html>
