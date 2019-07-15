"Mitra Menon" />


<title>INTROGRESS and Fold Enrichment analysis</title>

.kable-table table>thead>tr>th {
  border: none;
  border-bottom: 2px solid #dddddd;
}

.kable-table table>thead {
  background-color: #fff;
}
</style>


<div class="container-fluid main-container">

<!-- tabsets -->
<script>
$(document).ready(function () {
  window.buildTabsets("TOC");
});
</script>

<!-- code folding -->
<style type="text/css">
.code-folding-btn { margin-bottom: 4px; }
</style>
<script>
$(document).ready(function () {
  window.initializeSourceEmbed("AdapIntro_FE.Rmd");
  window.initializeCodeFolding("show" === "show");
});
</script>






<div class="fluid-row" id="header">

<div class="btn-group pull-right">
<button type="button" class="btn btn-default btn-xs dropdown-toggle" data-toggle="dropdown" aria-haspopup="true" aria-expanded="false"><span>Code</span> <span class="caret"></span></button>
<ul class="dropdown-menu" style="min-width: 50px;">
<li><a id="rmd-show-all-code" href="#">Show All Code</a></li>
<li><a id="rmd-hide-all-code" href="#">Hide All Code</a></li>
<li role="separator" class="divider"></li>
<li><a id="rmd-download-source" href="#">Download Rmd</a></li>
</ul>
</div>



<h1 class="title toc-ignore">INTROGRESS and Fold Enrichment analysis</h1>
<h4 class="author"><em>Mitra Menon</em></h4>

</div>


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

<div id="rmd-source-code">LS0tCnRpdGxlOiAiSU5UUk9HUkVTUyBhbmQgRm9sZCBFbnJpY2htZW50IGFuYWx5c2lzIgphdXRob3I6ICJNaXRyYSBNZW5vbiIKb3V0cHV0OiBodG1sX25vdGVib29rCi0tLQoKI0JBQ0tHUk9VTkQjCiMqT3ZlcmFsbCBnb2FsIGlzIHRvIGlkZW50aWZ5IGFkYXB0aXZlbHkgaW50cm9ncmVzc2VkIGxvY2kgYW5kIGRldGVybWluZSB0aGUgbWFqb3IgZW52aXJvbm1lbnRhbCBheGlzIGFzc29jaWF0ZWQgd2l0aCBhZGFwdGl2ZSBpbnRyb2dyZXNzaW9uKiMKCkhlcmUgSSB1dGlsaXplIHRoZSByZWdyZXNzaW9uIGFwcHJvYWNoIGltcGxlbWVudGVkIHdpdGhpbiBgSU5UUk9HUkVTU2AgdG8gZmlyc3QgaWRlbnRpZnkgbG9jaSBleGhpYml0aW5nIGV4Y2VwdGlvbmFsIHBhdHRlcm5zIG9mIGludHJvZ3Jlc3Npb24uIFRoZXNlIHdpbGwgdGhlbiBiZSBzdWJqZWN0ZWQgdG8gYSBmb2xkIGVucmljaG1lbnQgYXBwcm9hY2ggdG8gZGV0ZXJtaW5lIGlmIGxvY2kgYXNzb2NpYXRlZCB3aXRoIGEgZ2l2ZW4gZW52aXJvbm1lbnRhbCB2YXJpYWJsZSBhcmUgbW9yZSBsaWtlbHkgdG8gYmUgY2FuZGlkYXRlcyBmb3IgYWRhcHRpdmUgaW50cm9ncmVzc2lvbiB0aGFuIHRob3NlIGFzc29jaWF0ZWQgd2l0aCBvdGhlciBlbnZpcm9ubWVudGFsIHZhcmlhYmxlcy4gV2Ugd2lsbCBmb2N1c3NpbmcgcHJpbWFyaWx5IG9uIGVudmlyb25tZW50YWwgdmFyaWFibGVzIHRoYXQgYXJlIHN0cm9uZ2x5IGRpc3Rpbmd1aXNoZWQgYmV0d2VlbiB0aGUgdHdvIGh5YnJpZGl6aW5nIHNwZWNpZXMsIG9yIHRob3NlIGRyaXZpbmcgdGhlIG5pY2hlcyBvZiB0aGUgdHdvIHNwZWNpZXMgKGlmIHlvdSBoYXZlIHByaW9yIGluZm9ybWF0aW9uIG9uIHRoYXQpLiAKCmBgYHtyfQpsaWJyYXJ5KGRhdGEudGFibGUpCmxpYnJhcnkoaW50cm9ncmVzcykKYGBgCgojU1RFUCAxOiBQcmVwIGRhdGEgYW5kIHJ1biBpbnRyb2dyZXNzICgqRm9yIGxhcmdlIG51bWJlciBvZiBsb2NpIHRoaXMgbmVlZHMgdG8gYmUgcnVuIG9uIGEgY29tcHV0aW5nIGNsdXN0ZXIqKQoKTG9hZCBkYXRhCmBgYHtyfQpHZGF0YTwtZnJlYWQoIi4uL3NucERhdGEvSWxsdW1pbmFNaXhfMzQ4OS9yYXdEYXRhL21pbm9yMDEyLnR4dCIsc2VwPSJcdCIsaGVhZGVyPVQsZGF0YS50YWJsZT1GKQpQb3BzPC1yZWFkLnRhYmxlKCIuLi9zbnBEYXRhL1BvcEluZC50eHQiLGhlYWRlcj1ULHNlcD0iXHQiLHN0cmluZ3NBc0ZhY3RvcnMgPSBGQUxTRSkKcGFyZW50YWxzPC1yZWFkLnRhYmxlKCIuLi9zbnBEYXRhL3N1bW1hcnlTdGF0c18zNDg5L01BRl9wYXJlbnRhbHMudHh0IixoZWFkZXI9VCxzZXA9Ilx0IikKYGBgCgpEdWUgdG8gbG93IG92ZXJhbGwgZGlmZmVyZW50aWF0aW9uIGJldHdlZW4gdGhlIHR3byBzcGVjaWVzIEkgYW0gdXNpbmcsIEkgd2lsbCB1c2UgYSB2ZXJ5IGxvdyBhbGxlbGUgZnJlcXVlbmN5IGRpZmZlcmVuY2UgY3V0b2ZmIHRvIG9idGFpbiB0aGUgaW5wdXQgbG9jaS4gQXMgc3RhdGVzIGluIHRoZSBwYXBlciwgdGhlIGFwcHJvYWNoIHNob3VsZCBub3QgYmUgYmlhc2VkIGJ5IHRoaXMuIE1heGltaXNpbmcgdGhlIG51bWJlciBvZiBsb2NpIHdpbGwgYWxzbyBoZWxwIHRoZSBhc3Nlc21lbnQgb2YgZm9sZCBlbnJpY2htZW50LgpgYGB7cn0KYWxsZWxlRnJlcTwtY2JpbmQoZGlmZj1hYnMocGFyZW50YWxzJExQX3B1cmUtcGFyZW50YWxzJFNXV1BfcHVyZSkscGFyZW50YWxzKQphbGxlbGVGcmVxPC1hbGxlbGVGcmVxW2FsbGVsZUZyZXEkZGlmZj4wLCBdIAoKR2RhdGE8LUdkYXRhWyAsY29sbmFtZXMoR2RhdGEpJWluJWFsbGVsZUZyZXEkbG9jaV0KR2RhdGE8LWNiaW5kKFBvcHMsR2RhdGEpCkdkYXRhPC1HZGF0YVshKEdkYXRhJFBvcD09Ik1FU0MxTCIpLCBdClBvcHM8LVBvcHNbIShQb3BzJFBvcD09Ik1FU0MxTCIpLCBdCgpsb2NpPC1hcy5tYXRyaXgoR2RhdGFbICwtYygxOjUpXSkKbG9jaVtpcy5uYShsb2NpKV08LSJOQS9OQSIKbG9jaVtsb2NpPT0xXTwtIkEvRCIKbG9jaVtsb2NpPT0wXTwtIkEvQSIKbG9jaVtsb2NpPT0yXTwtIkQvRCIKYGBgCgpEZWZpbmUgZ3JvdXBzCmBgYHtyfQpMUDwtYygiQkNLIiwiQlUiLCJDSCIsIkNQIiwiSE8iLCJKRU4iLCJKViIsIkxQIiwiTUVMIiwiTVNMIiwiVkMiLCJDT0xPMkwiLCJDT0xPM0giKQpTV1dQPC1jKCJEQVYyTCIsIkFMVCIsIkhVQTJIIiwiSFVBMkwiLCJTQU5SSVQxTCIsIkJBTCIsIkJBUyIsIkNISSIsIkVQQSIsIkVSTSIsIkdVQSIsIkdVSSIsIlNBTSIsIlRPTSIsIlZFUiIsIllBUSIsIkNISTFIIiwiQ0hJMUwiLCJDSEkzIiwiREFWMUwiLCJEQVYxSCIsIkRBVjJIIikKcG9wOTg8LXJlYWQudGFibGUoIi4uL3NucERhdGEvc3VtbWFyeVN0YXRzXzM0ODkvYmF5ZW52T3V0Lzk4UG9wSURzIikKCkdkYXRhPC1jYmluZC5kYXRhLmZyYW1lKFBvcHMsbG9jaSkKClAxX3N3d3A8LUdkYXRhW0dkYXRhJFBvcCVpbiVTV1dQLCBdClAyX2xwPC1HZGF0YVtHZGF0YSRQb3AlaW4lTFAsIF0KbWl4ZWQ8LUdkYXRhW0dkYXRhJFBvcCVpbiVwb3A5OCRWMSwgXQpgYGAKCkRhdGEgcHJlcCBmb3IgcnVubmluZyBJTlRST0dSRVNTCmBgYHtyfQpQMV9zd3dwPC1QMV9zd3dwWyAsLWMoMTo1KV0KUDJfbHA8LVAyX2xwWyAsLWMoMTo1KV0KbWl4ZWQ8LWNiaW5kKFBvcD1taXhlZCRQb3AsSW5kPW1peGVkJFBvcEluZCxtaXhlZFsgLC1jKDE6NSldKQoKUDFfdHJhbnNwb3NlPC10KFAxX3N3d3ApClAyX3RyYW5zcG9zZTwtdChQMl9scCkKbWl4ZWRfdHJhbnNwb3NlPC10KG1peGVkKQoKCkxvY3VzPC1jYmluZChsb2N1cz1hbGxlbGVGcmVxJGxvY2ksdHlwZT1yZXAoIkMiLG5yb3coYWxsZWxlRnJlcSkpKQoKYWRtaXhDb3VudDwtcHJlcGFyZS5kYXRhKGFkbWl4Lmdlbj1taXhlZF90cmFuc3Bvc2UsbG9jaS5kYXRhPUxvY3VzLHBhcmVudGFsMT1QMV90cmFuc3Bvc2UscGFyZW50YWwyPVAyX3RyYW5zcG9zZSxwb3AuaWQ9VFJVRSxpbmQuaWQ9VFJVRSxmaXhlZD1GQUxTRSxzZXAucm93cyA9IEZBTFNFLHNlcC5jb2x1bW5zID0gRkFMU0UpCmBgYAoKRXN0aW1hdGUgaHlicmlkIGluZGV4IGFuZCBpbnRlcnNwZWNmaWMgaGV0LCB0aGlzIGlzIGltcG9ydGFudCB0byBkZXRlcm1pbmUgb3V0bGllciBsb2NpIHdpdGggcmVzcGVjdCB0byBnZW5vbWUgd2lkZSBhbmNlc3RyeS4gKCp0YWtlcyBhIGxvbmcgdGltZSBhbmQgc28gc2F2ZSB0aGUgb3V0cHV0IGFmdGVyIGl0IGlzIHJ1biopCmBgYHtyfQpISW5kZXg8LWVzdC5oKGludHJvZ3Jlc3MuZGF0YSA9IGFkbWl4Q291bnQsbG9jaS5kYXRhID0gTG9jdXMsaW5kLnRvdXNlID0gTlVMTCxmaXhlZCA9IEZBTFNFKQpoZWFkKEhJbmRleCkKI3dyaXRlLnRhYmxlKGNiaW5kKEhJbmRleCxtaXhlZFsgLDE6Ml0pLGZpbGU9Ikh5YnJpZEluZGV4OThQb3BzLnR4dCIscm93Lm5hbWVzPUYscXVvdGU9RixzZXA9Ilx0IikKCmludC5oZXQ8LWNhbGMuaW50ZXJzcC5oZXQoaW50cm9ncmVzcy5kYXRhPWFkbWl4Q291bnQpCmBgYAoKTm93LCBjb25kdWN0IGdlbm9taWMgY2xpbmUgYW5hbHlzaXMgdXNpbmcgdGhlIHBhcmFtZXRyaWMgYXBwcm9hY2ggc2luY2UgdGhlIFNOUHMgdXNlZCBkb24ndCBleGhpYml0IGZpeGVkIGRpZmZlcmVuY2VzICgqVGFrZXMgYSBsb25nIHRpbWUgYW5kIG5lZWRzIHRvIGJlIHJ1biBvbiB0aGUgY2x1c3RlciopCgpgYGB7cn0KR2NsaW5lc19wYXJhPC1nZW5vbWljLmNsaW5lcyhpbnRyb2dyZXNzLmRhdGEgPSBhZG1peENvdW50LGhpLmluZGV4ID0gSEluZGV4LGxvY2kuZGF0YSA9IExvY3VzLAogICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgbWV0aG9kPSJwYXJhbWV0cmljIixzaWcudGVzdD1UUlVFLGxvY2kudG91c2U9TlVMTCxpbmQudG91c2U9TlVMTCkKYGBgCgpXcml0ZSBvdXRwdXQgZmlsZXMgCmBgYHtyfQp3cml0ZS50YWJsZShHY2xpbmVzX3BhcmEkU3VtbWFyeS5kYXRhLCBmaWxlPSJTdW1tX3BhcmExMDAwLnR4dCIscXVvdGU9RkFMU0UsIHNlcD0iXHQiKQp3cml0ZS50YWJsZShHY2xpbmVzX3BhcmEkRml0dGVkLkFBLGZpbGU9IkhvbW9QMV9wYXJhMTAwMC50eHQiLHNlcD0iXHQiLHF1b3RlPUYpCndyaXRlLnRhYmxlKEdjbGluZXNfcGFyYSRGaXR0ZWQuYWEsZmlsZT0iSG9tb1AyX3BhcmExMDAwLnR4dCIsc2VwPSJcdCIscXVvdGU9RikKd3JpdGUudGFibGUoR2NsaW5lc19wYXJhJE5ldXRyYWwuQUEsZmlsZT0iSG9tb1AxX0NJcGFyYTEwMDAudHh0IixzZXA9Ilx0IixxdW90ZT1GKQp3cml0ZS50YWJsZShHY2xpbmVzX3BhcmEkTmV1dHJhbC5hYSxmaWxlPSJIb21vUDJfQ0lwYXJhMTAwMC50eHQiLHNlcD0iXHQiLHF1b3RlPUYpCndyaXRlLnRhYmxlKEdjbGluZXNfcGFyYSRRdWFudGlsZXMsZmlsZT0iUXVhbnRpbGVfcGFyYTEwMDAudHh0IixzZXA9Ilx0IixxdW90ZT1GKQpgYGAKClRoZSBwYXJhbWV0cmljIGFwcHJvYWNoIGlzIHByb25lIHRvIGZhbHNlIHBvc2l0aXZlcywgc3BlY2lmaWNhbGx5IGR1ZSB0byB0aGUgbGFyZ2UgbnVtYmVyIG9mIHRlc3RzIGRvbmUsIHdlIHdpbGwgY29uZHVjdCBhIFAudmFsIGNvcnJlY3Rpb24gYW5kIHVzZSBvbmx5IGxvY2kgdGhhdCBwYXNzIHAgdmFsIHRocmVzaG9sZC4KYGBge3J9ClN1bW08LXJlYWQudGFibGUoIi4uL3NucERhdGEvc3VtbWFyeVN0YXRzXzM0ODkvaW50cm9ncmVzcy9hbGwxMDAwcmVwcy9TdW1tX3BhcmExMDAwLnR4dCIsaGVhZGVyPVQsc2VwPSJcdCIpCgpudW08LTAuMDUvbnJvdyhTdW1tKSAjQm9uZmVyb25uaSBjb3JyZWN0aW9uClN1bW1CSDwtU3VtbVtTdW1tJFAudmFsdWU8bnVtLCBdCmBgYAoKQWNyb3NzIGFsbCBpbmRpdmlkdWFscyB3ZSB3aWxsIG9ubHkgcmV0YWluIGxvY2kgdGhhdCBoYXZlIHBhc3NlZCB0aGUgQm9uZmVyb25uaSBjb3JyZWN0aW9uLiBTaW5jZSB3ZSBoYXZlIGVzdGltYXRlcyBmb3IgZWFjaCBpbmRpdmlkdWFsLCB3ZSB3aWxsIHRoZW4gY2xhc3NpZnkgdGhlIGxvY2kgd2l0aGluIGFuIGluZGl2aWR1YWwgYXMgaGF2aW5nIGV4Y2VzcyBhbmNlc3RyeSBmcm9tIExQIHVzaW5nIHRoZSBDSSBlc3RpbWF0ZXMuCmBgYHtyfQpMUGNpPC1mcmVhZCgiLi4vc25wRGF0YS9zdW1tYXJ5U3RhdHNfMzQ4OS9pbnRyb2dyZXNzL2FsbExvY2kvSG9tb1AyX0NJcGFyYTEwMDAudHh0IixzZXA9Ilx0IixkYXRhLnRhYmxlPUYpCkxQb3V0PC1mcmVhZCgiLi4vc25wRGF0YS9zdW1tYXJ5U3RhdHNfMzQ4OS9pbnRyb2dyZXNzL2FsbExvY2kvSG9tb1AyX3BhcmExMDAwLnR4dCIsc2VwPSJcdCIsZGF0YS50YWJsZT1GKQoKTFBfdXBwZXI8LUxQY2lbICwxOm5jb2woTFBvdXQpXQoKCkV4Y2Vzc0xQPC1tYXRyaXgobnJvdyA9IG5yb3coTFBvdXQpLG5jb2w9OTUwKQoKZm9yIChjIGluIDI6bmNvbChMUG91dCkpewogIG48LWMtMQogIGZvciAociBpbiAxOm5yb3coTFBvdXQpKXsKICAgIEV4Y2Vzc0xQW3Isbl08LWlmZWxzZShMUG91dFtyLGNdPkxQX3VwcGVyW3IsY10sMSwwKQogIH0KfQoKcm93bmFtZXMoRXhjZXNzTFApPC1MUG91dCRWMQpjb2xuYW1lcyhFeGNlc3NMUCk8LWNvbG5hbWVzKExQb3V0KVsyOm5jb2woTFBvdXQpXQpFeGNlc3NMUF9vdXQ8LUV4Y2Vzc0xQW3Jvd25hbWVzKEV4Y2Vzc0xQKSVpbiVTdW1tQkgkbG9jdXMsIF0KCmBgYAoKQWZ0ZXIgYXNzaWduaW5nIGludHJvZ3Jlc3NlZCBsb2NpIGFzIDEgYW5kIG5vbi1pbnRyb2dyZXNzZWQgYXMgMCwgd2UgY2FuIGFkZCB0aGVtIGFjcm9zcyBhbGwgaW5kaXZpZHVhbHMgdG8gZGV0ZXJtaW5lIHdoYXQgcHJvcG9ydGlvbiBvZiB0aW1lcyBpdCB3YXMgaW50cm9ncmVzc2VkIGFjcm9zcyB0aGUgaHlicmlkIHpvbmUuIFdlIGNhbiB0aGVuIHVzZSBhIGN1dG9mZiAobW9zdGx5IGFyYml0YXJ5KSB0byBjbGFzc2lmeSBhIGxvY2kgYXMgc2lnbmlmaWNhbnRseSBpbnRyb2dyZXNzZWQuIFRoaXMgYXJiaXRhcnkgY3V0b2ZmIGZvciBoZXJlIGlzIHNldCBhcyAwLjIgKDIwJSBvZiBpbmRpdmlkaXVhbHMpIGFuZCB3YXMgYXNzZXNzZWQgdGhyb3VnaCBhIHNlcmllcyBvZiBjdXRvZmZzIGJlaW5nIHN1YmplY3RlZCB0byBkb3duc3RyZWFtIGFuYWx5c2lzLCBhbGwgb2Ygd2hpY2ggZ2F2ZSBzaW1pbGFyIHJlc3VsdHMuCmBgYHtyfQpFeGNlc3NMUF9vdXQ8LWFwcGx5KEV4Y2Vzc0xQLDEsZnVuY3Rpb24oWCkgcmV0dXJuKHN1bShYKS9sZW5ndGgoWCkpKQpFeGNlc3NMUF9vdXQyPC1FeGNlc3NMUF9vdXRbRXhjZXNzTFBfb3V0PjAuMjBdCmBgYAoKI1NURVAgMjogRGV0ZXJtaW5pbmcgd2hldGhlciB0aGUgaW50cm9ncmVzc2VkIGxvY2kgYXJlIGFkYXB0aXZlLgpIZXJlIHdlIHdpbGwgdXRpbGl6ZSB0aGUgb3V0bGllcnMgaWRlbnRpZmllZCB0aHJvdWdoIEdFQSAoYmF5ZW52IGluIG15IGNhc2UpIGFuZCBpbnRlcnNlY3QgdGhlbSB3aXRoIHRoZSBpbnRyb2dyZXNzZWQgbG9jaS4gU2luY2UgdGhlIG51bWJlciBvZiBvdXRsaWVycyBsb2NpIHZhcnkgYWNyb3NzIGVudmlyb25tZW50YWwgdmFyaWFibGVzLCB3ZSB3aWxsIHV0aWxpemUgYSBwZXJtdXRhdGlvbiBhcHByb2FjaCB0byBhc3Nlc3MgdGhlIHNpZ25pZmljYW5jZSBvZiBmb2xkIGVucmljaG1lbnQuCgpMb2FkIEdFQSBvdXRsaWVycyBTTlBzCmBgYHtyfQpjbGltPC1saXN0LmZpbGVzKCIuLi9zbnBEYXRhL3N1bW1hcnlTdGF0c18zNDg5L2JheWVudk91dC9jbGltL2NvbnZlcmdlbmNlLyIsZnVsbC5uYW1lcyA9IFQscGF0dGVybj0iT3ZlcmxhcCIpClNvaWw8LWxpc3QuZmlsZXMoIi4uL3NucERhdGEvc3VtbWFyeVN0YXRzXzM0ODkvYmF5ZW52T3V0L1NvaWwvY29udmVyZ2VuY2UvIixmdWxsLm5hbWVzID0gVCxwYXR0ZXJuPSJPdmVybGFwIikKb3V0bGllcnM8LWMoY2xpbSxTb2lsKQoKCmZpbGVzQkY8LXZlY3RvcigibGlzdCIsbGVuZ3RoKG91dGxpZXJzKSkKCmZvciAoZiBpbiAxOmxlbmd0aChmaWxlc0JGKSl7CiAgCiAgZGY8LXJlYWQudGFibGUob3V0bGllcnNbZl0saGVhZGVyPVQsc2VwPSJcdCIpCiAgZGY8LWRmWyAsMTo1XQogIGZpbGVzQkZbW2ZdXTwtZGYKfQoKZmlsZXM8LWxpc3QuZmlsZXMoIi4uL3NucERhdGEvc3VtbWFyeVN0YXRzXzM0ODkvYmF5ZW52T3V0L2NsaW0vY29udmVyZ2VuY2UvIixwYXR0ZXJuPSJPdmVybGFwIikKZmlsZXM8LWMoZmlsZXMsbGlzdC5maWxlcygiLi4vc25wRGF0YS9zdW1tYXJ5U3RhdHNfMzQ4OS9iYXllbnZPdXQvU29pbC9jb252ZXJnZW5jZS8iLHBhdHRlcm49Ik92ZXJsYXAiKSkKdmFySUQ8LXNhcHBseShzdHJzcGxpdChmaWxlcywiT3ZlcmxhcDFfM2NoYWluQ29udmciKSwiWyIsMikKdmFySUQ8LXNhcHBseShzdHJzcGxpdCh2YXJJRCwiLnR4dCIpLCJbIiwxKQoKCm5hbWVzKGZpbGVzQkYpPC12YXJJRApyZW1vdmU8LWMoIkVsZXZhdGlvbiIgLCJMYXRpdHVkZSIsIkxvbmdpdHVkZSIpCm91dGxpZXJzPC1maWxlc0JGWyEobmFtZXMoZmlsZXNCRiklaW4lcmVtb3ZlKV0KYGBgCgpTb21lIHNpbXBsZSBtYW5pcHVsYXRpb24gdG8gcmVhdGFpbiBvbmx5IHRoZSBTTlBzIHVzZWQgaW4gYElOVFJPR1JFU1NgCmBgYHtyfQpvdXRsaWVyczI8LWxhcHBseShvdXRsaWVycyxmdW5jdGlvbihkZikgZGZbZGYkbG9jaSVpbiVMUG91dCRWMSwgXSkgI05VTEwgU0VUCkJGX2dDbGluZXM8LWxhcHBseShvdXRsaWVyczIsZnVuY3Rpb24oZGYpIHJldHVybihkZltkZiRsb2NpJWluJW5hbWVzKEV4Y2Vzc0xQX291dDIpLCBdKSkgI1NIQVJFRCBTRVQKQkZfZ0NsaW5lc0lEPC11bmxpc3QobGFwcGx5KEJGX2dDbGluZXMsZnVuY3Rpb24oZGYpIG5yb3coZGYpKSkKYGBgCgpFc3RpbWF0aW5nIGZvbGQgY2hhbmdlIGZvciBlYWNoIGVudmlyb25tZW50YWwgdmFyaWFibGUKYGBge3J9Ck5yPC1CRl9nQ2xpbmVzSUQvdW5saXN0KGxhcHBseShvdXRsaWVyczIsZnVuY3Rpb24oZGYpIG5yb3coZGYpKSkKRG48LWxlbmd0aChFeGNlc3NMUF9vdXQyKS9ucm93KExQb3V0KQoKRkM8LU5yL0RuCmBgYAoKKkZ1bmN0aW9uIGZvciByYW5kb21pemF0aW9uIHRvIG9idGFpbiBudWxsIGRpc3RyaWJ1dGlvbiBvZiBGQyoKClRha2VzIDQgaW5wdXRzOgoKYGRmYCBpcyBhIGRhdGFmcmFtZSBvZiBiYXllbnYgb3V0bGllcnMsIHdoZXJlIG5yb3cgPSB0b3RhbCBvdXRsaWVycwpgbG9jaWAgaXMgYSBjaGFyYWN0ZXIgdmVjdG9yIG9mIElEcyBvZiBhbGwgc25wcyB1c2VkIGluIGludHJvZ3Jlc3MgYW5hbHlzaXMKYEludGAgdG90YWwgbnVtYmVyIG9mIGxvY2kgaWRlbnRpZmllZCBhcyBzaWduaWZpY2FudGx5IGludHJvZ3Jlc3NlZCBmcm9tIExQIGluIFNURVAgMQpgUmAgbnVtYmVyIG9mIGJvb3RzdHJhcCByZXBsaWNhdGVzIHRvIHJ1bgoKYGBge3J9CnJhbmRGQzwtZnVuY3Rpb24oZGYsbG9jaSxJbnQsUil7CgogIGJvb3Q8LU5VTEwKICBmb3IoaSBpbiAxOlIpewogICAgCiAgICBCRnJhbmQ8LXNhbXBsZShsb2NpLG5yb3coZGYpLHJlcGxhY2UgPSBGKQogICAgUEZyYW5kPC1zYW1wbGUobG9jaSxJbnQscmVwbGFjZSA9IEYpCiAgICBzaGFyZWQ8LWxlbmd0aChQRnJhbmRbUEZyYW5kJWluJUJGcmFuZF0pCiAgICAKICAgIGJvb3RbaV08LShzaGFyZWQvbGVuZ3RoKEJGcmFuZCkpLyhsZW5ndGgoUEZyYW5kKS9sZW5ndGgobG9jaSkpCiAgfQoKICByZXR1cm4oYm9vdCkKfQoKYGBgCgpUaGUgYWN0dWFsIHRlc3QKYGBge3J9CklEPC1MUG91dCRWMQpJbnQ8LWxlbmd0aChFeGNlc3NMUF9vdXQyKQoKCkZDcmFuZE91dDwtbGFwcGx5KG91dGxpZXJzMixmdW5jdGlvbihkZikgcmV0dXJuKHJhbmRGQyhkZixJRCxJbnQsMTAwMDApKSkKYGBgCgpDaGVjayB3aGljaCB2YXJpYWJsZXMgaGF2ZSBvYnNlcnZlZCBGQyBvdXRzaWRlIHRoZSA5OXRoIHBlcmNlbnRpbGUgb2YgbnVsbApgYGB7cn0KCmZvciAoaSBpbiAxOmxlbmd0aChGQ3JhbmRPdXQpKXsKICBFbnY8LW5hbWVzKEZDcmFuZE91dClbaV0KICBlbXA8LUZDcmFuZE91dFtbRW52XV0KICBpZiAocXVhbnRpbGUoZW1wLHByb2JzID0gYygwLjk5OSkpPD1GQ1tbRW52XV0pewogICAgY2F0KCJ2YXJpYWJsZSIsRW52LCJpcyBzaWduaWZpY2FudGx5IGVucmljaGVkIiwiXG4iKQogIH0KICAKfQoKYGBgCgo=</div>



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
