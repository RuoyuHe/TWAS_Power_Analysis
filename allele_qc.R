allele.qc = function(a1,a2,ref1,ref2) {
  ref = ref1
  flip = ref
  flip[ref == "A"] = "T"
  flip[ref == "T"] = "A"
  flip[ref == "G"] = "C"
  flip[ref == "C"] = "G"
  flip1 = flip
  
  ref = ref2
  flip = ref
  flip[ref == "A"] = "T"
  flip[ref == "T"] = "A"
  flip[ref == "G"] = "C"
  flip[ref == "C"] = "G"
  flip2 = flip;
  
  snp = list()
  snp[["keep"]] = (!((a1=="A" & a2=="T") | (a1=="T" & a2=="A") | (a1=="C" & a2=="G") | (a1=="G" & a2=="C"))) *
    ((a1 == ref2 & a2 == ref1) | (a1 == flip2 & a2 == flip1) | ((a1 == ref1 & a2 == ref2)) | (a1 == flip1 & a2 == flip2))
  snp[["keep"]] = as.logical(snp[["keep"]])
  snp[["flip"]] = (a1 == ref2 & a2 == ref1) | (a1 == flip2 & a2 == flip1)
  
  return(snp)
}
